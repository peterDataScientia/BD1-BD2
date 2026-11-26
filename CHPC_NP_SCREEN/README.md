# BRD4 natural-product screening on CHPC

Reuse the two existing CDK screening CSV files on Lustre as a source of identifiers and structures only. No CDK predictions are used as features. Expected 735,276 COCONUT and 276,401 LOTUS rows; cross-database duplicates are retained and not counted as unique chemistry. Every source row is retained, including failed structures with missing predictions and explicit validity flags.

The two SVRs use Morgan ECFP4 radius 2, 2048 bits and an RBF kernel. Development-only training (BD1 2,275; BD2 588); original holdout remains excluded from fitting. BD1 uses previously locked C=10, gamma=.01, epsilon=.1. Historical BD2 hyperparameters were not recoverable from the repository evaluation notebook: select anew using 4-fold scaffold grouped development CV with C=1,10,100 and gamma=.005,.01; epsilon=.1. Recalculate holdout metrics for each model and report them in the protocol rather than assuming historical results apply. This is a new BD2 fitted model, not a serialized historical model. No holdout-guided selection.

BD1 query preprocessing uses Cleanup and FragmentParent, matching the available BD1 curation notebook. BD2 query preprocessing uses RDKit isomeric canonical SMILES on original source smiles; upstream historical BD2 standardization has not been completely established and should be reviewed before final scientific interpretation. CDK standardized SMILES is retained for audit. AD is maximum ECFP4 Tanimoto to development structures. AD >= .30 is exploratory, not calibrated reliability. Known means exact canonical isomeric match to any curated labeled structure in the domain; it is not evidence of novelty against the literature. Selectivity differences are exploratory predictions, not measured selectivity. Cross-database structure deduplication and uncertainty-based ranking should follow screening.

## Pilot

From the existing BD1-BD2 checkout, use the same working Python environment used for the Y-randomization. Do not overwrite prior datasets or results. Check input paths and imports first.

```bash
ls -lh /mnt/lustre/users/apeter1/CDK_NP_SCREENING/CDK_{COCONUT,LOTUS}_full_SVR_ET.csv
export LOO_PYTHON="$(command -v python)"
"$LOO_PYTHON" -c 'import rdkit, numpy, pandas, sklearn, threadpoolctl; print("Environment OK")'
qsub -V -P CHEM1937 -v NP_SOURCE=COCONUT,NP_LIMIT_CHUNKS=2,NP_THREADS=4 -N NP_COC_P CHPC_NP_SCREEN/BRD4_NP.pbs
qsub -V -P CHEM1937 -v NP_SOURCE=LOTUS,NP_LIMIT_CHUNKS=2,NP_THREADS=4 -N NP_LOT_P CHPC_NP_SCREEN/BRD4_NP.pbs
```

Two chunks = 1,024 source rows per pilot. Verify protocol, heldout recheck CSVs, and chunk outputs before full screening. Full jobs omit NP_LIMIT_CHUNKS and resume existing chunks when protocol hashes match. Adjust walltime and CPU/memory based on observed runtime and CHPC queue policies. Never overlap jobs for the same source/protocol. Outputs live under /mnt/lustre/users/apeter1/BRD4_NP_SCREENING/SOURCE/protocol_hash. Full merge checks row identity and counts; summary JSON is written only after completion. The workflow has not been executed on CHPC yet.
