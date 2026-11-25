# CHPC nested leave-one-out validation

Run from repository root. Locked holdout structures are excluded using existing heldout CSV files, with scaffold disjointness checked. BD1: 2275 development, 541 holdout; BD2: 588 development, 187 holdout. Morgan ECFP4 2048 bits, RBF SVR, four-fold inner scaffold-grouped tuning (C 1, 10, 100; gamma .005, .01; epsilon .1). This grid matches BD1 development notebook and is a new declared protocol for BD2, not a claim of recovering its historical hyperparameters. 25 fits per outer fold. This is supplementary to independent scaffold holdout, not a replacement.

Install in a compatible Python environment: numpy pandas scikit-learn rdkit joblib. Check CHPC queue/resource syntax and environment first. The PBS templates assume smp, 24 CPUs, 48 GB, 3 hours. Do not run simultaneous jobs for the same domain/configuration.

```bash
export LOO_PYTHON="$(command -v python)"
"$LOO_PYTHON" -c 'import rdkit, sklearn, pandas, joblib'
qsub -v LOO_PYTHON,LOO_WORKERS=2,LOO_EXTRA="--limit 4" CHPC_LOO/BD1_LOO.pbs
# After pilot completes successfully:
qsub -v LOO_PYTHON CHPC_LOO/BD1_LOO.pbs
qsub -v LOO_PYTHON CHPC_LOO/BD2_LOO.pbs
qstat -u "$USER"
```

Resubmit same script after walltime if incomplete. Each completed fold is saved atomically under CHPC_LOO/outputs/domain/protocol_hash; do not delete this directory. Inputs, code and dependency versions identify checkpoints. Final metrics and predictions appear only once all folds complete. Copy outputs to persistent storage before cluster shutdown. No CHPC execution or timing benchmark has yet been performed.
