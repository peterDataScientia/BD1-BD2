# BRD4 BD1/BD2 Y-randomization on CHPC

200 deterministic label permutations per target plus one unpermuted baseline (402 evaluations total). The original scaffold-disjoint holdouts are excluded from modelling; only their SMILES membership is read. The data files, fingerprint and grid follow CHPC_LOO. This is a different outer validation design from LOOCV: five fixed scaffold-grouped outer folds, each with four scaffold-grouped inner folds, tuning all six SVR configurations by mean inner RMSE. Compare the unpermuted and permuted pooled OOF Q2 from this exact protocol, **not** historical LOOCV Q2. For acyclic molecules use the same unique-SMILES convention as CHPC_LOO.

Each evaluation requires 5*(6*4+1)=125 SVR fits; total 50,250. Permutation p-value = (1 + number of randomized Q2 >= original Q2)/201. Global label shuffling tests the null of exchangeable activity labels, not scaffold-specific or assay-batch confounding.

## CHPC submission

Copy this new directory to your existing checkout containing BD1, BD2 and CHPC_LOO, preserving earlier outputs. From repository root with yesterday's Python environment activated:

```bash
export LOO_PYTHON="$(command -v python)"
"$LOO_PYTHON" -c 'import rdkit, sklearn, pandas, joblib; print("Environment OK")'
# Pilot: baseline (0) plus four shuffled runs for each target
qsub -V -P CHEM1937 -q serial -N YR_BD1_P -v YR_DOMAIN=BD1,YR_START=0,YR_END=4,YR_WORKERS=4 CHPC_Y_RANDOM/Y_RANDOM.pbs
qsub -V -P CHEM1937 -q serial -N YR_BD2_P -v YR_DOMAIN=BD2,YR_START=0,YR_END=4,YR_WORKERS=4 CHPC_Y_RANDOM/Y_RANDOM.pbs
qstat -u "$USER"
```

After confirming pilot completion and resources, submit both full runs. Existing completed permutations are skipped automatically.

```bash
qsub -V -P CHEM1937 -q serial -N YR_BD1 -v YR_DOMAIN=BD1,YR_WORKERS=24 CHPC_Y_RANDOM/Y_RANDOM.pbs
qsub -V -P CHEM1937 -q serial -N YR_BD2 -v YR_DOMAIN=BD2,YR_WORKERS=24 CHPC_Y_RANDOM/Y_RANDOM.pbs
qstat -u "$USER"
```

If PBS Pro arrays are supported on the queue, an optional alternative **instead of** the full jobs is `qsub -J 1-20%4 -V -P CHEM1937 -q serial -N YR_BD1 -v YR_DOMAIN=BD1,YR_WORKERS=10 CHPC_Y_RANDOM/Y_RANDOM.pbs` (and BD2 equivalent). Each array index handles 10 randomized runs; ensure baseline 0 has completed first. Do not run overlapping ranges or submit both alternatives. Queue limits and timing must be checked locally.

Atomic checkpoints appear in CHPC_Y_RANDOM/outputs/DOMAIN/protocol_hash/. Preserve the same script, datasets and dependency versions to resume. Do not run jobs with overlapping permutation ranges for the same domain. When all 201 evaluations finish:

```bash
"$LOO_PYTHON" CHPC_Y_RANDOM/run_y_randomization.py BD1 --summarize
"$LOO_PYTHON" CHPC_Y_RANDOM/run_y_randomization.py BD2 --summarize
cat CHPC_Y_RANDOM/outputs/*/*/summary.json
```

Outputs: protocol.json, original_OOF_predictions.csv, perm_000.json through perm_200.json, Y_randomization_results.csv and summary.json. An incomplete run cannot produce a final permutation p-value. Copy outputs off Lengau before shutdown. No CHPC job is submitted from GitHub and no actual permutation results are claimed until the jobs run.
