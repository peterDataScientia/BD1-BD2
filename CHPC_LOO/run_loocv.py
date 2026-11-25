import argparse, os, json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.svm import SVR
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.metrics.pairwise import euclidean_distances
from joblib import Parallel, delayed
import sklearn, rdkit

def main():
    p=argparse.ArgumentParser()
    p.add_argument('domain', choices=['BD1','BD2'])
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--limit', type=int)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    files={'BD1':('BD1/data/BRD4_BD1_QSAR_model_ready_2816.csv','BD1/results/heldout_predictions.csv','standardized_smiles',2816,541),'BD2':('BD2/data/BRD4_BD2_STRICT_MODELING_775.csv','BD2/validation/heldout_predictions_with_AD.csv','canonical_isomeric_smiles',775,187)}
    data,hold,col,n,h=files[a.domain]
    d=pd.read_csv(root/data); test=pd.read_csv(root/hold)
    assert len(d)==n and len(test)==h and d[col].is_unique and test[col].is_unique
    assert set(test[col]) <= set(d[col])
    d=d.loc[~d[col].isin(test[col])].reset_index(drop=True)
    assert len(d)==n-h
    y=d.pIC50.to_numpy(float); assert np.isfinite(y).all()
    mols=[Chem.MolFromSmiles(s) for s in d[col]]
    assert all(m is not None for m in mols)
    scaffold=lambda m,s: MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False) or 'ACYCLIC_'+s
    groups=np.array([scaffold(m,s) for m,s in zip(mols,d[col])])
    assert set(groups).isdisjoint({scaffold(Chem.MolFromSmiles(s),s) for s in test[col]})
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    fps=[gen.GetFingerprint(m) for m in mols]
    X=np.array([gen.GetFingerprintAsNumPy(m) for m in mols],dtype=np.float64)
    dist=euclidean_distances(X,X,squared=True)
    kernels={g:np.exp(-g*dist) for g in (.005,.01)}
    # Nested grouped tuning avoids using each outer test label for parameter selection.
    meta={'domain':a.domain,'n':len(d),'grid':{'C':[1,10,100],'gamma':[.005,.01],'epsilon':[.1]},'inner_folds':4,'rdkit':rdkit.__version__,'sklearn':sklearn.__version__,'numpy':np.__version__,'script_sha':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':[hashlib.sha256((root/f).read_bytes()).hexdigest() for f in (data,hold)]}
    key=hashlib.sha256(json.dumps(meta,sort_keys=True).encode()).hexdigest()[:16]
    out=root/'CHPC_LOO'/'outputs'/a.domain/key; out.mkdir(parents=True,exist_ok=True)
    (out/'protocol.json').write_text(json.dumps(meta,indent=2))
    def fit(K,C,tr,va):
        model=SVR(kernel='precomputed',C=C,epsilon=.1)
        return model.fit(K[np.ix_(tr,tr)],y[tr]).predict(K[np.ix_(va,tr)])
    def one(i):
        target=out/f'fold_{i:05d}.json'
        if target.exists(): return
        tr=np.delete(np.arange(len(y)),i)
        splits=[(tr[t],tr[v]) for t,v in GroupKFold(4).split(tr,y[tr],groups[tr])]
        scores=[]
        for g,K in kernels.items():
            for C in (1,10,100):
                score=np.mean([np.sqrt(np.mean((fit(K,C,t,v)-y[v])**2)) for t,v in splits])
                scores.append((float(score),g,C))
        score,g,C=min(scores)
        pred=float(fit(kernels[g],C,tr,np.array([i]))[0])
        sim=max(DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in tr]))
        row={'index':i,'smiles':d[col].iloc[i],'observed_pIC50':float(y[i]),'predicted_pIC50':pred,'nearest_training_tanimoto':sim,'C':C,'gamma':g,'inner_RMSE':score}
        tmp=target.with_suffix('.tmp'); tmp.write_text(json.dumps(row)); os.replace(tmp,target)
        print(f'completed {a.domain} fold {i+1}/{len(y)}',flush=True)
    pending=[i for i in range(len(y)) if not (out/f'fold_{i:05d}.json').exists()]
    if a.limit is not None: pending=pending[:a.limit]
    Parallel(n_jobs=a.workers,prefer='processes',batch_size=1)(delayed(one)(i) for i in pending)
    paths=sorted(out.glob('fold_*.json'))
    print(f'Saved {len(paths)}/{len(y)} folds',flush=True)
    if len(paths)!=len(y): return
    r=pd.DataFrame([json.loads(f.read_text()) for f in paths]).sort_values('index')
    assert r['index'].tolist()==list(range(len(y)))
    obs=r.observed_pIC50; pred=r.predicted_pIC50
    r['residual_pIC50']=pred-obs
    r.to_csv(out/'LOO_predictions.csv',index=False)
    metrics={'n':len(r),'Q2_LOO':float(r2_score(obs,pred)),'RMSE':float(np.sqrt(np.mean((pred-obs)**2))),'MAE':float(mean_absolute_error(obs,pred)),'bias':float(np.mean(pred-obs))}
    (out/'metrics.json').write_text(json.dumps(metrics,indent=2))
    print(metrics)
if __name__=='__main__': main()
