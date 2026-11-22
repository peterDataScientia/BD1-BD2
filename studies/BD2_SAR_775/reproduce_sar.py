"""Reproduce the BRD4-BD2 SAR feasibility screen from the original 775-compound CSV."""
import csv, json, itertools
from pathlib import Path
from collections import defaultdict,Counter
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMMPA
from rdkit.Chem.Scaffolds import MurckoScaffold
inputfile=Path("data/BRD4_BD2_STRICT_MODELING_775.csv")
if not inputfile.exists(): inputfile=Path("BRD4_BD2_STRICT_MODELING_775.csv")
if not inputfile.exists(): raise FileNotFoundError("Put original 775-compound CSV in data/ or cwd")
with inputfile.open(newline="",encoding="utf-8-sig") as f: rows=list(csv.DictReader(f))
mols=[Chem.MolFromSmiles(x["canonical_isomeric_smiles"]) for x in rows]
assert len(rows)==775 and all(mols)
p=[float(x["pIC50"]) for x in rows]
assays=[set(x["assay_ids"].split(";")) for x in rows]
scaffolds=[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m)) for m in mols]
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fps=[fpgen.GetFingerprint(m) for m in mols]
sim_pairs=[]
for i in range(len(rows)):
    for k,t in enumerate(DataStructs.BulkTanimotoSimilarity(fps[i],fps[i+1:]),i+1):
        if t>=.8:
            sim_pairs.append(dict(i=i,j=k,similarity=t,delta=abs(p[i]-p[k]),
                shared=bool(assays[i]&assays[k])))
idx=defaultdict(list)
for i,m in enumerate(mols):
    seen=set()
    for _, fragments in rdMMPA.FragmentMol(m,maxCuts=1,resultsAsMols=False):
        if not fragments or fragments.count(".")!=1:continue
        a,b=fragments.split(".")
        ma,mb=Chem.MolFromSmiles(a),Chem.MolFromSmiles(b)
        if ma is None or mb is None:continue
        na,nb=ma.GetNumHeavyAtoms(),mb.GetNumHeavyAtoms()
        core,var=(a,b) if na>=nb else (b,a)
        nc,nv=(na,nb) if na>=nb else (nb,na)
        if not (1<=nv<=12 and nc>=2*nv):continue
        core=Chem.MolToSmiles(Chem.MolFromSmiles(core))
        var=Chem.MolToSmiles(Chem.MolFromSmiles(var))
        if (core,var) in seen:continue
        seen.add((core,var));idx[core].append((i,var))
mmp={}
for core,entries in idx.items():
    for (i,a),(j,b) in itertools.combinations(entries,2):
        if i==j or a==b:continue
        key=tuple(sorted((i,j)))
        if key not in mmp or len(core)>len(mmp[key]):
            mmp[key]=core
stringent=0
for i,j in mmp:
    if assays[i]&assays[j] and DataStructs.TanimotoSimilarity(fps[i],fps[j])>=.8 and abs(p[i]-p[j])>=1:
        stringent+=1
result=dict(n_compounds=len(rows),scaffolds=len(set(scaffolds)),
 singleton_scaffolds=sum(n==1 for n in Counter(scaffolds).values()),
 n_assays=len(set().union(*assays)),
 n_single_measurement=sum(int(x["n_measurements"])==1 for x in rows),
 similarity_pairs=len(sim_pairs),
 candidate_cliffs=sum(x["delta"]>=1 for x in sim_pairs),
 same_assay_candidate_cliffs=sum(x["delta"]>=1 and x["shared"] for x in sim_pairs),
 matched_pairs=len(mmp),strict_mmp_cliffs=stringent)
print(json.dumps(result,indent=2))
Path("BD2_SAR_summary.json").write_text(json.dumps(result,indent=2))
