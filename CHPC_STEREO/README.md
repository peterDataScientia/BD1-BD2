# Chirality-aware BRD4 QSAR comparison

This benchmark retains both existing curated datasets and locked holdout membership. No existing output is overwritten. It compares conventional Morgan ECFP4 (includeChirality=False) with chiral Morgan ECFP4 (includeChirality=True), both radius 2, 2048 bits and RBF-SVR. Each representation is independently tuned with the same development-only grid (C=1,10,100; gamma=.005,.01; epsilon=.1). Nested CV uses the same five outer scaffold folds and four inner scaffold folds as CHPC_Y_RANDOM. Final development-only models use four-fold grouped tuning. The earlier holdout is a matched secondary comparison, not a fresh untouched external test.

Outputs: stereo_audit.csv (potential and specified stereo and hashed fingerprints), stereochemical_pairs.csv (same non-isomeric canonical SMILES, distinct isomeric SMILES, fingerprint collisions, measured activity differences and predicted pair differences), matched_predictions.csv (OOF/holdout predictions and matching AD), audit_summary.json, protocol.json, summary.json and both serialized development models with fingerprints in training order. Non-isomeric SMILES grouping can also collapse isotopic differences; inspect flagged structures before interpreting all pairs as stereoisomers. Pair scores use fully specified pairs in the same evaluation fold or holdout, with group-weighted absolute difference errors; aggregated activity differences may also reflect assay heterogeneity. Legacy acyclic scaffold grouping is retained for comparability; cross-partition and cross-fold stereo pair counts are audited, and nonzero counts require a stricter connectivity-grouped sensitivity analysis before generalization claims. Aggregate improvement alone does not establish stereoisomer discrimination; inspect pair counts and paired errors. If pair data are sparse, report the limitation.

## CHPC

From existing BD1-BD2 checkout, with the same verified Python module:

```bash
module load python/3.9.6
export LOO_PYTHON="$(command -v python)"
"$LOO_PYTHON" -m py_compile CHPC_STEREO/run_stereo.py
"$LOO_PYTHON" -c 'from rdkit import Chem; from rdkit.Chem import rdFingerprintGenerator as F; print(Chem.FindPotentialStereo(Chem.MolFromSmiles("C[C@H](O)F"))); F.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True); print("Stereo API OK")'
qsub -V -P CHEM1937 -N STEREO_BD1 -v STEREO_DOMAIN=BD1 CHPC_STEREO/STEREO.pbs
qsub -V -P CHEM1937 -N STEREO_BD2 -v STEREO_DOMAIN=BD2 CHPC_STEREO/STEREO.pbs
qstat -u "$USER"
```

Review outputs under CHPC_STEREO/outputs/{BD1,BD2}/protocol_hash before selecting any new screening model. No new million-record screen or new permutation run is launched by this benchmark. No CHPC execution has been performed yet. Copy results before Lengau shutdown.
