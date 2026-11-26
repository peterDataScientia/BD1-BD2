"""Secondary chirality sensitivity screen of the fixed conventional AD>=0.30 subset."""
import hashlib
import json
import os
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
import rdkit, sklearn
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

ROOT=Path(__file__).resolve().parents[1]
BASE=Path("/mnt/lustre/users/apeter1/BRD4_NP_SCREENING")
INPUTS={"COCONUT":("c0c511883de106ba",735276,13076),"LOTUS":("0e67265c7593eb2a",276401,2546)}
GEN={c:rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=c) for c in (False,True)}
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""): h.update(b)
    return h.hexdigest()
def jsonwrite(p,d):
    tmp=p.with_suffix(".tmp")
    tmp.write_text(json.dumps(d,indent=2,allow_nan=False))
    os.replace(tmp,p)
def csvwrite(p,d):
    tmp=p.with_name(p.name+".tmp")
    d.to_csv(tmp,index=False,compression="gzip" if p.suffix==".gz" else None)
    os.replace(tmp,p)
def get_models():
    models={}; meta={}
    for domain in ("BD1","BD2"):
        dirs=[p.parent for p in (ROOT/"CHPC_STEREO/outputs"/domain).glob("*/summary.json")]
        if len(dirs)!=1: raise RuntimeError("Expected exactly one completed benchmark for "+domain+": "+str(dirs))
        folder=dirs[0]
        s=json.loads((folder/"summary.json").read_text())
        meta[domain]={"folder":str(folder),"summary_sha256":sha(folder/"summary.json")}
        for label,c in (("achiral",False),("chiral",True)):
            path=folder/(label+"_model.joblib")
            m=joblib.load(path) # locally generated trusted model
            assert m["domain"]==domain and m["includeChirality"]==c and m["radius"]==2 and m["fpSize"]==2048
            assert m["protocol"]["versions"]["rdkit"]==rdkit.__version__ and m["protocol"]["versions"]["sklearn"]==sklearn.__version__
            assert all(sha(ROOT/p)==h for p,h in m["protocol"]["inputs"].items())
            assert m["training_fingerprints"].shape==(m["protocol"]["development_n"],2048)
            assert m["model"].shape_fit_[0]==len(m["training_fingerprints"])
            assert all(m[k]==v for k,v in s["models"][label]["final_params"].items())
            models[domain,label]=m
            meta[domain][label]={"sha256":sha(path),"C":m["C"],"gamma":m["gamma"],"epsilon":m["epsilon"],
                                 "holdout":s["models"][label]["holdout"]}
        assert models[domain,"achiral"]["training_SMILES"]==models[domain,"chiral"]["training_SMILES"]
    return models,meta
def score(mols,m):
    T=m["training_fingerprints"].astype(np.float32)
    pred=[]; ad=[]
    for i in range(0,len(mols),256):
        X=np.stack([GEN[m["includeChirality"]].GetFingerprintAsNumPy(z) for z in mols[i:i+256]]).astype(np.float32)
        dot=X@T.T
        cnt=X.sum(1)[:,None]+T.sum(1)[None,:]
        dist=np.maximum(cnt-2*dot,0)
        den=cnt-dot
        pred.extend(m["model"].predict(np.exp(-m["gamma"]*dist.astype(np.float64))))
        ad.extend(np.divide(dot,den,out=np.zeros_like(dot),where=den>0).max(1))
    return np.array(pred),np.array(ad)
def stereo(m):
    info=Chem.FindPotentialStereo(m)
    return sum(str(x.specified)!="Specified" for x in info)
def connectivity(m):
    # Strip stereo but retain isotopes, unlike isomericSmiles=False.
    z=Chem.Mol(m); Chem.RemoveStereochemistry(z)
    return Chem.MolToSmiles(z,isomericSmiles=True)
def known_connectivity(domain,m):
    inputs=m["protocol"]["inputs"]
    path=next(p for p in inputs if "/data/" in p)
    col={"BD1":"standardized_smiles","BD2":"canonical_isomeric_smiles"}[domain]
    return {connectivity(Chem.MolFromSmiles(s)) for s in pd.read_csv(ROOT/path,usecols=[col])[col]}
def main():
    models,meta=get_models()
    chunks=[]; inpmeta={}
    for source,(key,n,want) in INPUTS.items():
        folder=BASE/source/key
        path=folder/(source+"_BRD4_predictions.csv.gz")
        protocol=json.loads((folder/"protocol.json").read_text())
        for domain in ("BD1","BD2"):
            old=protocol["models"][domain]
            m=models[domain,"achiral"]
            assert all(old[k]==m[k] for k in ("C","gamma","epsilon"))
            assert old["input_sha256"]==m["protocol"]["inputs"]
            assert old["development_SMILES_sha256"]==hashlib.sha256("\n".join(m["training_SMILES"]).encode()).hexdigest()
        count=0; selected=[]
        for x in pd.read_csv(path,chunksize=50000,dtype={"identifier":str}):
            assert np.array_equal(x.source_row.to_numpy(),np.arange(count,count+len(x)))
            assert x.source.eq(source).all()
            count+=len(x)
            selected.append(x.loc[x.unseen_both & (x.min_AD>=.3)].copy())
        x=pd.concat(selected,ignore_index=True)
        assert count==n and len(x)==want,(source,count,len(x))
        chunks.append(x)
        inpmeta[source]={"path":str(path),"sha256":sha(path),"rows":count,"eligible":len(x)}
    d=pd.concat(chunks,ignore_index=True)
    assert len(d)==15622
    protocol={"inputs":inpmeta,"models":meta,"script_sha256":sha(Path(__file__)),
              "selection":"original unseen_both and conventional min_AD>=0.30; not exhaustive chiral screen",
              "rank_rule":"agreement at both representations; conservative minimum potency and directional gap, then joint AD",
              "dedup_rule":"identical ordered pair of canonical isomeric BD1 and BD2 scoring structures; retain all source aliases",
              "limitations":"No experimentally calibrated uncertainty or natural origin verification; thresholds exploratory; original BD2 upstream standardization incomplete"}
    key=hashlib.sha256(json.dumps(protocol,sort_keys=True).encode()).hexdigest()[:16]
    out=BASE/"STEREO_SECONDARY"/key
    out.mkdir(parents=True,exist_ok=True)
    jsonwrite(out/"protocol.json",protocol)
    for dom in ("BD1","BD2"):
        mols=[Chem.MolFromSmiles(s) for s in d[dom+"_SMILES"]]
        assert all(m is not None for m in mols)
        a,ad=score(mols,models[dom,"achiral"])
        err=float(np.max(abs(a-d[dom+"_pIC50"].to_numpy())))
        aderr=float(np.max(abs(ad-d[dom+"_AD"].to_numpy())))
        print(dom,"baseline maximum absolute prediction / AD difference",err,aderr,flush=True)
        assert err<1e-5 and aderr<1e-6,"Original screening predictions not reproduced"
        p,sim=score(mols,models[dom,"chiral"])
        d[dom+"_chiral_pIC50"]=p
        d[dom+"_chiral_IC50_nM"]=10**(9-p)
        d[dom+"_chiral_AD"]=sim
        d[dom+"_chiral_minus_achiral"]=p-d[dom+"_pIC50"]
        d[dom+"_unspecified_stereo"]= [stereo(m) for m in mols]
        known=known_connectivity(dom,models[dom,"achiral"])
        d[dom+"_known_connectivity"]=[connectivity(m) in known for m in mols]
        d[dom+"_canonical_key"]=[Chem.MolToSmiles(m,isomericSmiles=True) for m in mols]
        protocol[dom+"_baseline_recheck"]={"max_pIC50_error":err,"max_AD_error":aderr}
    d["chiral_min_AD"]=d[["BD1_chiral_AD","BD2_chiral_AD"]].min(axis=1)
    d["joint_min_AD"]=d[["min_AD","chiral_min_AD"]].min(axis=1)
    d["chiral_delta"]=d.BD1_chiral_pIC50-d.BD2_chiral_pIC50
    d["achiral_delta"]=d.BD1_pIC50-d.BD2_pIC50
    d["max_representation_disagreement"]=d[["BD1_chiral_minus_achiral","BD2_chiral_minus_achiral"]].abs().max(axis=1)
    d["BD1_conservative_pIC50"]=d[["BD1_pIC50","BD1_chiral_pIC50"]].min(axis=1)
    d["BD2_conservative_pIC50"]=d[["BD2_pIC50","BD2_chiral_pIC50"]].min(axis=1)
    d["dual_conservative_pIC50"]=d[["BD1_conservative_pIC50","BD2_conservative_pIC50"]].min(axis=1)
    d["BD1_conservative_gap"]=d[["achiral_delta","chiral_delta"]].min(axis=1)
    d["BD2_conservative_gap"]=-d[["achiral_delta","chiral_delta"]].max(axis=1)
    d["consensus_BD1_selective"]=(d.BD1_conservative_pIC50>=6)&(d.BD1_conservative_gap>=1)
    d["consensus_BD2_selective"]=(d.BD2_conservative_pIC50>=6)&(d.BD2_conservative_gap>=1)
    d["consensus_dual"]=(d.dual_conservative_pIC50>=6)
    d["consensus_dual_high"]=(d.dual_conservative_pIC50>=7)
    d["known_connectivity_either"]=d.BD1_known_connectivity|d.BD2_known_connectivity
    d["unspecified_stereo_either"]=(d.BD1_unspecified_stereo>0)|(d.BD2_unspecified_stereo>0)
    d["manual_exclusion_reason"]=np.where(d.identifier.eq("CNP0577166.1"),"previously flagged synthetic BET-like compound; verify provenance","")
    d["chemical_key"]=[hashlib.sha256((a+"\n"+b).encode()).hexdigest() for a,b in zip(d.BD1_canonical_key,d.BD2_canonical_key)]
    # Propagate exclusion to aliases with identical scoring structures.
    flagged=set(d.loc[d.manual_exclusion_reason.ne(""),"chemical_key"])
    d.loc[d.chemical_key.isin(flagged),"manual_exclusion_reason"]="previously flagged synthetic BET-like structure; verify provenance"
    d["provenance_status"]="UNVERIFIED"
    csvwrite(out/"all_15622_records.csv.gz",d)
    csvwrite(out/"source_aliases.csv",d[["chemical_key","source","source_row","identifier","smiles"]])
    u=d.sort_values(["joint_min_AD","source","source_row"],ascending=[False,True,True]).drop_duplicates("chemical_key").copy()
    u["source_record_count"]=u.chemical_key.map(d.chemical_key.value_counts())
    csvwrite(out/"unique_scoring_structures.csv.gz",u)
    counts=[]
    criteria={"BD1_selective":("consensus_BD1_selective","BD1_conservative_gap","BD1_conservative_pIC50"),
              "BD2_selective":("consensus_BD2_selective","BD2_conservative_gap","BD2_conservative_pIC50"),
              "dual":("consensus_dual","dual_conservative_pIC50",None),
              "dual_high":("consensus_dual_high","dual_conservative_pIC50",None)}
    for t in (.3,.4,.5):
        eligible=u[u.joint_min_AD>=t]
        clean=eligible[(~eligible.known_connectivity_either)&eligible.manual_exclusion_reason.eq("")]
        counts.append(dict(AD=t,original_records=int((d.min_AD>=t).sum()),
            joint_AD_records=int((d.joint_min_AD>=t).sum()),joint_AD_unique=len(eligible),
            unique_excluding_known_connectivity_and_flagged=len(clean),
            **{k:int(eligible[v[0]].sum()) for k,v in criteria.items()},
            **{"priority_"+k:int(clean[v[0]].sum()) for k,v in criteria.items()}))
        for name,(flag,primary,secondary) in criteria.items():
            x=clean[clean[flag]].copy()
            sort=[primary]+([secondary] if secondary else [])+["joint_min_AD","max_representation_disagreement","chemical_key"]
            x=x.sort_values(sort,ascending=[False]*(len(sort)-2)+[True,True])
            x.insert(0,"rank",np.arange(1,len(x)+1))
            csvwrite(out/("%s_AD_%s.csv"%(name,str(t))),x)
    csvwrite(out/"counts.csv",pd.DataFrame(counts))
    jsonwrite(out/"protocol.json",protocol)
    summary={"input_records":len(d),"unique_scoring_structure_pairs":len(u),"known_connectivity_unique":int(u.known_connectivity_either.sum()),
             "unspecified_stereo_unique":int(u.unspecified_stereo_either.sum()),"flagged_unique":int(u.manual_exclusion_reason.ne("").sum()),
             "threshold_counts":counts,"note":"Overlapping categories; priority lists exclude known connectivity and flagged structure, not all synthetic molecules. Origin unverified. Agreement is not a confidence interval."}
    jsonwrite(out/"summary.json",summary)
    print("FINAL",out,json.dumps(summary,indent=2),flush=True)
if __name__=="__main__": main()
