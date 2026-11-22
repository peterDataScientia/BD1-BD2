"""BD2 SAR assay-stratified permutation: run after reproduce_sar.py."""
import csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
source=Path("data/BRD4_BD2_STRICT_MODELING_775.csv")
with source.open(newline="",encoding="utf-8-sig") as f: records=list(csv.DictReader(f))
p=np.array([float(x["pIC50"]) for x in records])
groups=defaultdict(list)
for i,r in enumerate(records):
    if int(r["n_assays"])==1:groups[r["assay_ids"]].append(i)
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fps=[fpgen.GetFingerprint(Chem.MolFromSmiles(r["canonical_isomeric_smiles"])) for r in records]
pairs=[]
for g,inds in groups.items():
    for k,i in enumerate(inds):
        for j in inds[k+1:]:
            t=DataStructs.TanimotoSimilarity(fps[i],fps[j])
            if t>=0.70:pairs.append((i,j,t))
rng=np.random.default_rng(20261122)
strong=[(i,j) for i,j,t in pairs if t>=.80]
i,j=np.array(strong).T
observed=int(np.sum(np.abs(p[i]-p[j])>=1))
counts=[]
for _ in range(2000):
    shuffled=p.copy()
    for inds in groups.values():
        shuffled[inds]=p[rng.permutation(inds)]
    counts.append(int(np.sum(np.abs(shuffled[i]-shuffled[j])>=1)))
a=np.asarray(counts)
stats=dict(observed=observed,n_pairs=len(strong),null_mean=float(a.mean()),null_sd=float(a.std(ddof=1)),null_interval_95=[float(x) for x in np.quantile(a,[.025,.975])],p_lower=float((1+np.sum(a<=observed))/(len(a)+1)))
stats["sensitivity"]=[dict(threshold=t,n_pairs=sum(sim>=t for _,_,sim in pairs),n_cliffs=sum(sim>=t and abs(p[i]-p[j])>=1 for i,j,sim in pairs)) for t in [.7,.75,.8,.85,.9]]
print(json.dumps(stats,indent=2));Path("BD2_permutation_results.json").write_text(json.dumps(stats,indent=2))
