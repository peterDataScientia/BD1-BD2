"""Restartable natural-product screen of existing CDK input rows; no CDK scores used."""
import argparse, hashlib, json, os, time, shutil, platform
from pathlib import Path
import numpy as np
import pandas as pd
import rdkit, sklearn
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.svm import SVR
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from threadpoolctl import threadpool_limits

RDLogger.DisableLog("rdApp.error")
FILES={
 "BD1":("BD1/data/BRD4_BD1_QSAR_model_ready_2816.csv","BD1/results/heldout_predictions.csv","standardized_smiles",2816,541),
 "BD2":("BD2/data/BRD4_BD2_STRICT_MODELING_775.csv","BD2/validation/heldout_predictions_with_AD.csv","canonical_isomeric_smiles",775,187)}
INPUTS={"COCONUT":("CDK_COCONUT_full_SVR_ET.csv",735276),"LOTUS":("CDK_LOTUS_full_SVR_ET.csv",276401)}
GEN=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""): h.update(b)
    return h.hexdigest()

def write_json(p,obj):
    t=p.with_name(p.name+".tmp")
    t.write_text(json.dumps(obj,indent=2,allow_nan=False))
    os.replace(t,p)

def fp(m):
    return GEN.GetFingerprintAsNumPy(m).astype(np.float32)

def canonical(m):
    return Chem.MolToSmiles(m,isomericSmiles=True)

def bd1_parent(m):
    return rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(m))

def group(m,s):
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False) or "ACYCLIC_"+s

def cross(X,T):
    # binary dot products are exact in float32 (integers <=2048)
    dot=X@T.T
    count=X.sum(axis=1)[:,None]+T.sum(axis=1)[None,:]
    dist=np.maximum(count-2*dot,0)
    den=count-dot
    sim=np.divide(dot,den,out=np.zeros_like(dot),where=den>0)
    return dist,sim.max(axis=1)

def kernel(dist,g):
    return np.exp(-g*dist.astype(np.float64))

def metrics(y,p):
    return {"R2":float(r2_score(y,p)),"RMSE":float(np.sqrt(np.mean((p-y)**2))),
            "MAE":float(mean_absolute_error(y,p)),"bias":float(np.mean(p-y))}

def fit_model(root,domain):
    data,hold,col,n,h=FILES[domain]
    d=pd.read_csv(root/data)
    t=pd.read_csv(root/hold,usecols=[col])
    assert len(d)==n and len(t)==h and d[col].is_unique and t[col].is_unique
    assert set(t[col])<=set(d[col])
    tr=d.loc[~d[col].isin(t[col])].reset_index(drop=True)
    te=d.loc[d[col].isin(t[col])].reset_index(drop=True)
    assert len(tr)==n-h
    m=[Chem.MolFromSmiles(s) for s in tr[col]]
    mt=[Chem.MolFromSmiles(s) for s in te[col]]
    assert all(z is not None for z in m+mt)
    groups=np.array([group(z,s) for z,s in zip(m,tr[col])])
    assert set(groups).isdisjoint({group(z,s) for z,s in zip(mt,te[col])})
    X=np.stack([fp(z) for z in m]); Xt=np.stack([fp(z) for z in mt])
    y=tr.pIC50.to_numpy(float); yt=te.pIC50.to_numpy(float)
    assert np.isfinite(y).all() and np.isfinite(yt).all()
    dist,_=cross(X,X)
    scores=[]
    if domain=="BD1":
        C,g=10.,.01  # locked development parameters in BD1/results/summary.json
    else:
        # Historical BD2 hyperparameters not recovered: declare new development-only selection.
        folds=list(GroupKFold(4).split(X,y,groups))
        for g0 in (.005,.01):
            K=kernel(dist,g0)
            for c0 in (1.,10.,100.):
                rms=[]
                for a,b in folds:
                    assert set(groups[a]).isdisjoint(groups[b])
                    mm=SVR(kernel="precomputed",C=c0,epsilon=.1).fit(K[np.ix_(a,a)],y[a])
                    pr=mm.predict(K[np.ix_(b,a)])
                    rms.append(float(np.sqrt(np.mean((pr-y[b])**2))))
                scores.append((float(np.mean(rms)),g0,c0))
        _,g,C=min(scores)
    model=SVR(kernel="precomputed",C=C,epsilon=.1).fit(kernel(dist,g),y)
    td,ad=cross(Xt,X)
    pred=model.predict(kernel(td,g))
    met=metrics(yt,pred)
    meta={"domain":domain,"development_n":len(tr),"holdout_n":len(te),"C":C,"gamma":g,
          "epsilon":.1,"selection":"locked BD1" if domain=="BD1" else "new BD2 4-fold grouped development tuning",
          "BD2_grid_scores":[{"RMSE":s,"gamma":gg,"C":cc} for s,gg,cc in scores],
          "heldout_metrics":met,"input_sha256":{data:sha(root/data),hold:sha(root/hold)},
          "development_SMILES_sha256":hashlib.sha256("\n".join(tr[col]).encode()).hexdigest(),
          "standardization":"RDKit Cleanup + FragmentParent" if domain=="BD1" else "RDKit canonical isomeric SMILES; historical upstream BD2 standardization not fully recovered"}
    print(domain,json.dumps(meta),flush=True)
    # Known means matches any curated labeled structure, including untouched holdout.
    known={canonical(z) for z in m+mt}
    return (model,X,g,known,meta),pd.DataFrame({"SMILES":te[col],"observed_pIC50":yt,"predicted_pIC50":pred,"AD":ad})

def score_batch(mols,model):
    sv,T,g,known,_=model
    valid=[i for i,m in enumerate(mols) if m is not None]
    p=np.full(len(mols),np.nan); ad=p.copy()
    smiles=[""]*len(mols); k=np.zeros(len(mols),dtype=bool)
    if valid:
        X=np.stack([fp(mols[i]) for i in valid])
        dist,sim=cross(X,T)
        p[valid]=sv.predict(kernel(dist,g)); ad[valid]=sim
        for i in valid:
            smiles[i]=canonical(mols[i]); k[i]=smiles[i] in known
    return p,ad,smiles,k

def safe_mol(s):
    try: return Chem.MolFromSmiles(str(s)) if pd.notna(s) and str(s).strip() else None
    except Exception: return None

def safe_parent(m):
    try:
        z=bd1_parent(m) if m is not None else None
        return z if z is not None and z.GetNumAtoms()>0 else None
    except Exception: return None

def main():
    a=argparse.ArgumentParser()
    a.add_argument("source",choices=INPUTS)
    a.add_argument("--input-dir",type=Path,default=Path("/mnt/lustre/users/apeter1/CDK_NP_SCREENING"))
    a.add_argument("--output-dir",type=Path,default=Path("/mnt/lustre/users/apeter1/BRD4_NP_SCREENING"))
    a.add_argument("--chunk-size",type=int,default=512)
    a.add_argument("--threads",type=int,default=4)
    a.add_argument("--limit-chunks",type=int)
    args=a.parse_args()
    assert args.chunk_size>0 and args.threads>0 and (args.limit_chunks is None or args.limit_chunks>0)
    root=Path(__file__).resolve().parents[1]
    name,n=INPUTS[args.source]; inp=args.input_dir/name
    start=time.monotonic()
    with threadpool_limits(limits=args.threads):
        models={}; held={}
        for d in FILES: models[d],held[d]=fit_model(root,d)
        protocol={"source":args.source,"expected_rows":n,"input_sha256":sha(inp),
                  "chunk_size":args.chunk_size,"script_sha256":sha(Path(__file__)),
                  "models":{d:models[d][-1] for d in FILES},
                  "versions":{"python":platform.python_version(),"rdkit":rdkit.__version__,"sklearn":sklearn.__version__,"numpy":np.__version__,"pandas":pd.__version__},
                  "query":"original smiles column; BD1 Cleanup/FragmentParent; BD2 canonical isomeric; CDK SMILES retained for provenance",
                  "AD":"max ECFP4 binary Tanimoto to development compounds only"}
        # Version key excludes floating-point metrics, but includes chosen parameters, input hashes and versions.
        keymeta={k:v for k,v in protocol.items() if k!="models"}
        keymeta["models"]={d:{k:v for k,v in models[d][-1].items() if k not in ("heldout_metrics","BD2_grid_scores")} for d in FILES}
        key=hashlib.sha256(json.dumps(keymeta,sort_keys=True).encode()).hexdigest()[:16]
        out=args.output_dir/args.source/key
        parts=out/"parts"; parts.mkdir(parents=True,exist_ok=True)
        write_json(out/"protocol.json",protocol)
        for d in FILES: held[d].to_csv(out/(d+"_heldout_recheck.csv"),index=False)
        print("Output:",out,flush=True)
        done=0; total=0; expected=[]
        for i,raw in enumerate(pd.read_csv(inp,usecols=["identifier","smiles","SMILES"],chunksize=args.chunk_size,dtype=str,keep_default_na=False)):
            path=parts/("%08d.csv.gz"%total)
            expected.append(path)
            offset=total; total+=len(raw)
            if path.is_file(): continue
            if args.limit_chunks is not None and done>=args.limit_chunks:
                print("Pilot stopped after",done,"new chunks",flush=True); return
            # No input records discarded; invalid rows are explicitly flagged.
            mols=[safe_mol(s) for s in raw.smiles]
            parents=[safe_parent(m) for m in mols]
            r=raw.rename(columns={"SMILES":"CDK_standardized_SMILES"}).copy()
            r.insert(0,"source",args.source);r.insert(1,"source_row",np.arange(offset,offset+len(raw)))
            for d,m in [("BD1",parents),("BD2",mols)]:
                p,ad,s,k=score_batch(m,models[d])
                r[d+"_SMILES"]=s; r[d+"_pIC50"]=p; r[d+"_IC50_nM"]=10**(9-p)
                r[d+"_AD"]=ad; r[d+"_known"]=k; r[d+"_valid"]=np.isfinite(p)
            r["Delta_pIC50_BD1_minus_BD2"]=r.BD1_pIC50-r.BD2_pIC50
            r["min_AD"]=r[["BD1_AD","BD2_AD"]].min(axis=1,skipna=False)
            r["unseen_both"]=r.BD1_valid & r.BD2_valid & ~r.BD1_known & ~r.BD2_known
            tmp=path.with_name(path.name+".tmp")
            r.to_csv(tmp,index=False,compression="gzip"); os.replace(tmp,path)
            done+=1
            if done%20==0 or done==1: print("Rows reached",total,"of",n,"new chunks",done,"elapsed",round(time.monotonic()-start,1),flush=True)
        assert total==n,(total,n)
        assert set(parts.glob("*.csv.gz"))==set(expected)
        # Merge with strict sequential row audit, keep all input records.
        final=out/(args.source+"_BRD4_predictions.csv.gz")
        tmp=final.with_name(final.name+".tmp")
        import gzip
        counts={"rows":0,"valid_both":0,"unseen_both":0}
        for threshold in (.3,.4,.5): counts["unseen_both_AD_"+str(threshold)]=0
        with gzip.open(tmp,"wt",newline="") as f:
            for j,path in enumerate(expected):
                r=pd.read_csv(path,keep_default_na=True,dtype={"identifier":str})
                assert np.array_equal(r.source_row.to_numpy(),np.arange(counts["rows"],counts["rows"]+len(r)))
                counts["rows"]+=len(r);counts["valid_both"]+=int((r.BD1_valid&r.BD2_valid).sum())
                counts["unseen_both"]+=int(r.unseen_both.sum())
                for threshold in (.3,.4,.5): counts["unseen_both_AD_"+str(threshold)]+=int((r.unseen_both & (r.min_AD>=threshold)).sum())
                r.to_csv(f,index=False,header=j==0)
        assert counts["rows"]==n
        os.replace(tmp,final)
        counts["seconds_this_run"]=time.monotonic()-start
        write_json(out/"summary.json",counts)
        print(json.dumps(counts),flush=True)
if __name__=="__main__":main()
