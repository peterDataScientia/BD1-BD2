"""Standalone assay and compound dependency audit from original 775-compound source."""
import csv
from collections import Counter
from pathlib import Path
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
P=Path(__file__).resolve().parent
with (P/"data/BRD4_BD2_STRICT_MODELING_775.csv").open(newline="",encoding="utf-8-sig") as f: rows=list(csv.DictReader(f))
assert len(rows)==775
mol=[Chem.MolFromSmiles(r["canonical_isomeric_smiles"]) for r in rows]
assert all(mol)
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
chigen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
fps=[fpgen.GetFingerprint(m) for m in mol]
chfps=[chigen.GetFingerprint(m) for m in mol]
p=[float(r["pIC50"]) for r in rows]
candidates=[]
for i in range(len(rows)):
 for j,sim in enumerate(DataStructs.BulkTanimotoSimilarity(fps[i],fps[i+1:]),i+1):
  shared=set(rows[i]["assay_ids"].split(";"))&set(rows[j]["assay_ids"].split(";"))
  if sim<.80 or abs(p[i]-p[j])<1 or not shared:continue
  candidates.append(dict(compound_1=rows[i]["molecule_chembl_ids"],compound_2=rows[j]["molecule_chembl_ids"],assay=";".join(sorted(shared)),document_1=rows[i]["document_ids"],document_2=rows[j]["document_ids"],fold=10**abs(p[i]-p[j]),similarity=sim,chiral_similarity=DataStructs.TanimotoSimilarity(chfps[i],chfps[j])))
assert len(candidates)==63
out=P/"results";out.mkdir(exist_ok=True)
with (out/"standalone_candidate_cliffs.csv").open("w",newline="") as f:
 writer=csv.DictWriter(f,fieldnames=candidates[0].keys());writer.writeheader();writer.writerows(candidates)
counts=Counter(x["assay"] for x in candidates)
print("cliffs",len(candidates),"assays",len(counts),"top4",sum(v for _,v in counts.most_common(4)))
print("unique compounds",len({x["compound_1"] for x in candidates}|{x["compound_2"] for x in candidates}))
print(counts.most_common())
