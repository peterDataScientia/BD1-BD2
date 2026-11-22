"""Companion audit for the original curated BRD4-BD2 775 data."""
from pathlib import Path
from collections import Counter
import pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
P=Path(__file__).resolve().parent
d=pd.read_csv(P/"data/BRD4_BD2_STRICT_MODELING_775.csv")
c=pd.read_csv(P/"results/priority_cliff_provenance_audit.csv")
lookup=dict(zip(d.molecule_chembl_ids,d.canonical_isomeric_smiles))
fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
fp_no=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
def scores(r):
 a=Chem.MolFromSmiles(lookup[r.compound_1]);b=Chem.MolFromSmiles(lookup[r.compound_2])
 return pd.Series({"chirality_sensitive_Tanimoto":DataStructs.TanimotoSimilarity(fp.GetFingerprint(a),fp.GetFingerprint(b)),"chirality_agnostic_Tanimoto":DataStructs.TanimotoSimilarity(fp_no.GetFingerprint(a),fp_no.GetFingerprint(b))})
c=pd.concat([c,c.apply(scores,axis=1)],axis=1)
assays=c.groupby(["shared_assays","document_1"],as_index=False).agg(n_cliffs=("priority_rank","count"),n_mmp=("also_single_cut_MMP","sum"),max_fold=("fold_difference","max")).sort_values("n_cliffs",ascending=False)
assays.to_csv(P/"results/cliff_concentration_by_assay.csv",index=False)
c.to_csv(P/"results/priority_cliff_chirality_audit.csv",index=False)
c.sort_values("fold_difference",ascending=False).drop_duplicates("shared_assays").to_csv(P/"results/one_priority_pair_per_assay.csv",index=False)
print("unique_compounds",len(set(c.compound_1)|set(c.compound_2)))
print("assays",len(assays),"top_four_share",assays.head(4).n_cliffs.sum()/len(c))
