# CHPC Y-randomization results — 26 November 2026

Actual results from the nested scaffold-based CV permutation workflow in the parent directory. The source CHPC archive was BRD4_Y_RANDOMIZATION_RESULTS.tar.gz (BD1 protocol 005b5dd7e32ac81f; BD2 protocol 373d7f3fcfb6a574).

- BD1: 2,275 development compounds; 200 randomized evaluations plus original Q2 = 0.7119243194296495.
- BD2: 588 development compounds; 200 randomized evaluations plus original Q2 = 0.6922306239685689.
- Neither domain had a randomized score at or above its original Q2. One-sided permutation p = 1/201 = 0.004975124378109453 for each.

## Files

- BD1_Q2.csv and BD2_Q2.csv: actual per-permutation Q2 scores from the CHPC export, index 0 is the original unshuffled baseline, 1–200 are randomized labels. Counts, mean and sample SD were checked against CHPC summaries.
- Y_randomization_summary.csv: all reported CHPC summary metrics, including RMSE, MAE and permutation p-values.
- figures/: regenerated SVG visualizations from the committed Q2 values (not byte-identical copies of the original Matplotlib exports).
- Y_randomization_figures.ipynb: Figentra reproduction notebook; validates Q2 summaries and exports PNG/PDF/SVG.

The original CHPC tar.gz and original figure ZIP remain separate backups; their binary bytes were not uploaded through the available GitHub connector. Full individual RMSE/MAE/bias/timing, fold tuning choices and OOF predictions are in the original CHPC archive, not in the Q2-only CSV files here. Retain both original archives.

The permutation test uses nested scaffold CV (five outer folds and four inner folds), not LOOCV. The independent holdout was excluded. Results support rejection of the exchangeable-label null under this protocol, not independent experimental validation or absence of all assay confounding.
