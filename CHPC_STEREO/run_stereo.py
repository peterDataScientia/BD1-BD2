"""Matched conventional/chiral Morgan RBF-SVR benchmark. No screening or data changes."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
import rdkit, sklearn
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupKFold
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.svm import SVR

FILES={
 "BD1":("BD1/data/BRD4_BD1_QSAR_model_ready_2816.csv","BD1/results/heldout_predictions.csv","standardized_smiles",2816,541),
 "BD2":("BD2/data/BRD4_BD2_STRICT_MODELING_775.csv","BD2/validation/heldout_predictions_with_AD.csv","canonical_isomeric_smiles",775,187)}
GRID=list(itertools.product((.005,.01),(1.,10.,100.)))
EPS=.1
SEED=20261126

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def atomic_json(p,x):
    t=p.with_suffix(".tmp")
    t.write_text(json.dumps(x,indent=2,allow_nan=False))
    os.replace(t,p)

def groups_for(mols,smiles):
    # Reproduce the previous validation scaffold convention exactly.
    return np.array([MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)
                     or "ACYCLIC_"+s for m,s in zip(mols,smiles)])

def metric(y,p):
    e=p-y
    return dict(R2=float(1-np.sum(e**2)/np.sum((y-y.mean())**2)),
                RMSE=float(np.sqrt(np.mean(e**2))), MAE=float(np.mean(np.abs(e))),
                bias=float(np.mean(e)))

def tune(K,y,tr,splits):
    scores=[]
    for g,c in GRID:
        rms=[]
        for a,b in splits:
            m=SVR(kernel="precomputed",C=c,epsilon=EPS).fit(K[g][np.ix_(a,a)],y[a])
            p=m.predict(K[g][np.ix_(b,a)])
            rms.append(np.sqrt(np.mean((p-y[b])**2)))
        scores.append((float(np.mean(rms)),g,c))
    best=min(scores)
    return best,[dict(gamma=g,C=c,RMSE=r) for r,g,c in scores]

def folds(tr,groups,n):
    return [(tr[a],tr[b]) for a,b in GroupKFold(n).split(tr,groups=groups[tr])]

def stereo_info(m):
    # Include potential atom and double-bond stereo, not just tetrahedral centers.
    infos=Chem.FindPotentialStereo(m)
    specified=sum(str(x.specified)=="Specified" for x in infos)
    return len(infos),specified

def bitkey(x):
    return np.packbits(x.astype(np.uint8)).tobytes()

def audit(d,mols,Xs,partition,fold_id,out):
    conn=[Chem.MolToSmiles(m,isomericSmiles=False) for m in mols]
    iso=[Chem.MolToSmiles(m,isomericSmiles=True) for m in mols]
    st=[stereo_info(m) for m in mols]
    y=d.pIC50.to_numpy(float)
    rows=pd.DataFrame(dict(index=np.arange(len(d)),isomeric_SMILES=iso,
          connectivity_SMILES=conn,partition=partition,outer_fold=fold_id,
          pIC50=y,potential_stereo=[a for a,b in st],specified_stereo=[b for a,b in st]))
    rows["unspecified_stereo"]=rows.potential_stereo-rows.specified_stereo
    for label,X in Xs.items():
        rows[label+"_fingerprint_sha256"]=[hashlib.sha256(bitkey(x)).hexdigest() for x in X]
    rows.to_csv(out/"stereo_audit.csv",index=False)
    pairs=[]; groups=0
    for key,frame in rows.groupby("connectivity_SMILES",sort=False):
        if frame.isomeric_SMILES.nunique()<2: continue
        groups+=1
        for a,b in itertools.combinations(frame.index,2):
            if iso[a]==iso[b]: continue
            same_a=np.array_equal(Xs["achiral"][a],Xs["achiral"][b])
            same_c=np.array_equal(Xs["chiral"][a],Xs["chiral"][b])
            # Both fully specified: comparison is a stereochemical pair, not missing stereo.
            fully=all(st[i][0]>0 and st[i][0]==st[i][1] for i in (a,b))
            same_eval=partition[a]==partition[b] and (partition[a]=="holdout" or fold_id[a]==fold_id[b])
            pairs.append(dict(index_a=a,index_b=b,SMILES_a=iso[a],SMILES_b=iso[b],
                 connectivity_SMILES=key,partition_a=partition[a],partition_b=partition[b],
                 fold_a=int(fold_id[a]),fold_b=int(fold_id[b]),fully_specified_pair=fully,
                 same_evaluation_set=bool(same_eval),achiral_collision=bool(same_a),
                 chiral_collision=bool(same_c),observed_delta=float(y[a]-y[b]),
                 abs_observed_delta=float(abs(y[a]-y[b]))))
    cols=["index_a","index_b","SMILES_a","SMILES_b","connectivity_SMILES","partition_a","partition_b",
          "fold_a","fold_b","fully_specified_pair","same_evaluation_set","achiral_collision","chiral_collision",
          "observed_delta","abs_observed_delta"]
    p=pd.DataFrame(pairs,columns=cols)
    p.to_csv(out/"stereochemical_pairs.csv",index=False)
    report=dict(n=len(rows),n_with_potential_stereo=int((rows.potential_stereo>0).sum()),
        n_with_unspecified_stereo=int((rows.unspecified_stereo>0).sum()),
        connectivity_groups_with_distinct_stereo_records=groups,stereo_record_pairs=len(p),
        fully_specified_pairs=int(p.fully_specified_pair.sum()),
        achiral_colliding_pairs=int(p.achiral_collision.sum()),
        collisions_resolved_by_chirality=int((p.achiral_collision.astype(bool)&~p.chiral_collision.astype(bool)).sum()),
        pairs_crossing_holdout=int((p.partition_a!=p.partition_b).sum()),
        development_pairs_crossing_outer_folds=int(((p.partition_a=="development")&(p.partition_b=="development")&(p.fold_a!=p.fold_b)).sum()))
    # Any fingerprint-identical structures across partitions can inflate apparent generalization.
    for label in Xs:
        col=label+"_fingerprint_sha256"
        report[label+"_cross_partition_fingerprint_keys"]=len(set(rows.loc[partition=="development",col])&set(rows.loc[partition=="holdout",col]))
    atomic_json(out/"audit_summary.json",report)
    return p,report

def main():
    a=argparse.ArgumentParser()
    a.add_argument("domain",choices=FILES)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    data,hold,col,n,h=FILES[args.domain]
    d=pd.read_csv(root/data)
    test=pd.read_csv(root/hold,usecols=[col])
    assert len(d)==n and len(test)==h and d[col].is_unique and test[col].is_unique
    assert set(test[col])<=set(d[col])
    # Keep original development row order, then original holdout row order.
    is_test=d[col].isin(test[col])
    d=pd.concat([d.loc[~is_test],d.loc[is_test]],ignore_index=True)
    dev=np.arange(n-h); te=np.arange(n-h,n)
    partition=np.array(["development"]*(n-h)+["holdout"]*h)
    mols=[Chem.MolFromSmiles(s) for s in d[col]]
    assert all(m is not None for m in mols)
    y=d.pIC50.to_numpy(float)
    assert np.isfinite(y).all()
    groups=groups_for(mols,d[col])
    assert set(groups[dev]).isdisjoint(groups[te])
    outer=folds(dev,groups,5)
    fold_id=np.full(n,-1,dtype=int)
    for f,(tr,va) in enumerate(outer):
        fold_id[va]=f
        assert set(groups[tr]).isdisjoint(groups[va])
    assert (fold_id[dev]>=0).all()
    meta=dict(domain=args.domain,development_n=len(dev),holdout_n=len(te),
              inputs={data:sha(root/data),hold:sha(root/hold)},script_sha256=sha(Path(__file__)),
              versions=dict(rdkit=rdkit.__version__,sklearn=sklearn.__version__,numpy=np.__version__,pandas=pd.__version__),
              fold_sha256=hashlib.sha256(fold_id[dev].astype("<i8").tobytes()).hexdigest(),
              scaffold_rule="legacy Murcko achiral; acyclic unique by SMILES",
              grid=[dict(gamma=g,C=c,epsilon=EPS) for g,c in GRID],
              outer_folds=5,inner_folds=4,representation="Morgan radius 2; 2048 binary bits; includeChirality False vs True",
              holdout_note="Previously evaluated locked split; not a fresh external test. No holdout-guided tuning.")
    key=hashlib.sha256(json.dumps(meta,sort_keys=True).encode()).hexdigest()[:16]
    out=root/"CHPC_STEREO"/"outputs"/args.domain/key
    out.mkdir(parents=True,exist_ok=True)
    atomic_json(out/"protocol.json",meta)
    print("Output:",out,flush=True)
    Xs={}
    for label,chiral in [("achiral",False),("chiral",True)]:
        gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=chiral)
        Xs[label]=np.stack([gen.GetFingerprintAsNumPy(m) for m in mols]).astype(np.float64)
    pairs,report=audit(d,mols,Xs,partition,fold_id,out)
    print("Audit:",json.dumps(report),flush=True)
    predtable=pd.DataFrame(dict(index=np.arange(n),SMILES=d[col],partition=partition,
                               scaffold=groups,outer_fold=fold_id,observed_pIC50=y))
    results={}
    for label,X in Xs.items():
        dist=euclidean_distances(X,X,squared=True)
        np.maximum(dist,0,out=dist)
        K={g:np.exp(-g*dist) for g in (.005,.01)}
        del dist
        pred=np.full(n,np.nan)
        choices=[]
        for f,(tr,va) in enumerate(outer):
            (rmse,g,c),scores=tune(K,y,tr,folds(tr,groups,4))
            m=SVR(kernel="precomputed",C=c,epsilon=EPS).fit(K[g][np.ix_(tr,tr)],y[tr])
            pred[va]=m.predict(K[g][np.ix_(va,tr)])
            choices.append(dict(fold=f,C=c,gamma=g,inner_RMSE=rmse,grid_scores=scores))
        (rmse,g,c),scores=tune(K,y,dev,folds(dev,groups,4))
        m=SVR(kernel="precomputed",C=c,epsilon=EPS).fit(K[g][np.ix_(dev,dev)],y[dev])
        pred[te]=m.predict(K[g][np.ix_(te,dev)])
        assert np.isfinite(pred).all()
        predtable[label+"_predicted_pIC50"]=pred
        # AD calculated using the matching fingerprint and only the relevant training fold.
        ad=np.full(n,np.nan)
        for tr,va in outer+[(dev,te)]:
            dot=X[va]@X[tr].T
            den=X[va].sum(1)[:,None]+X[tr].sum(1)[None,:]-dot
            ad[va]=np.divide(dot,den,out=np.zeros_like(dot),where=den>0).max(1)
        predtable[label+"_AD"]=ad
        results[label]={"nested_scaffold_CV":metric(y[dev],pred[dev]),"holdout":metric(y[te],pred[te]),
                        "final_params":dict(C=c,gamma=g,epsilon=EPS),"development_grid_scores":scores,
                        "outer_fold_choices":choices}
        # Prediction requires an RBF cross-kernel to these fingerprints in this exact order.
        temp=out/(label+"_model.tmp")
        joblib.dump(dict(model=m,training_fingerprints=X[dev].astype(np.uint8),training_SMILES=d[col].iloc[dev].tolist(),
                         radius=2,fpSize=2048,includeChirality=(label=="chiral"),gamma=g,C=c,epsilon=EPS,
                         domain=args.domain,protocol=meta),temp,compress=3)
        os.replace(temp,out/(label+"_model.joblib"))
        del K
        print(label,json.dumps(results[label]),flush=True)
    predtable.to_csv(out/"matched_predictions.csv",index=False)
    if len(pairs):
        for label in Xs:
            p=predtable[label+"_predicted_pIC50"].to_numpy()
            pairs[label+"_predicted_delta"]=p[pairs.index_a.to_numpy(int)]-p[pairs.index_b.to_numpy(int)]
            pairs[label+"_delta_absolute_error"]=abs(pairs[label+"_predicted_delta"]-pairs.observed_delta)
        pairs.to_csv(out/"stereochemical_pairs.csv",index=False)
    # Evaluate pair deltas only for fully specified stereoisomers in the same held-out fold.
    eligible=pairs.loc[pairs.fully_specified_pair.astype(bool)&pairs.same_evaluation_set.astype(bool)].copy()
    pair_results={}
    for part in ("development","holdout"):
        p=eligible.loc[eligible.partition_a==part]
        pair_results[part]={"n_pairs":len(p),"n_connectivity_groups":int(p.connectivity_SMILES.nunique())}
        if len(p):
            for label in Xs:
                diff=p[label+"_predicted_delta"]-p.observed_delta
                nonzero=p.abs_observed_delta>=1.0
                # Equal-weight connectivity groups to avoid large stereochemical series dominating.
                tmp=p.assign(ae=abs(diff))
                pair_results[part][label]={"pair_delta_MAE":float(abs(diff).mean()),
                    "group_weighted_delta_MAE":float(tmp.groupby("connectivity_SMILES").ae.mean().mean()),
                    "n_pairs_observed_gap_ge_1":int(nonzero.sum()),
                    "direction_accuracy_gap_ge_1":float((np.sign(p.loc[nonzero,label+"_predicted_delta"])==np.sign(p.loc[nonzero,"observed_delta"])).mean()) if nonzero.any() else None}
    summary={"domain":args.domain,"audit":report,"models":results,"stereo_pair_evaluation":pair_results,
             "delta_chiral_minus_achiral":{part:{k:results["chiral"][part][k]-results["achiral"][part][k] for k in ("R2","RMSE","MAE")} for part in ("nested_scaffold_CV","holdout")}}
    atomic_json(out/"summary.json",summary)
    print("FINAL:",json.dumps(summary,indent=2),flush=True)

if __name__=="__main__": main()
