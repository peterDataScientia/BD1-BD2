# BD2

BRD4 bromodomain 2 (BD2) curated data, QSAR models, applicability-domain analyses, and FDA screening results.

This folder organizes materials for the BRD4 BD1/BD2 research project.

## Notebooks
- [BD2 QSAR modeling](notebooks/01_BD2_QSAR_Modeling.ipynb)
- [Validation-derived errors](notebooks/02_BD2_Validation_Derived_Errors.ipynb)

## Updated executed QSAR figures and metrics

[01_BD2_QSAR_Modeling.ipynb](notebooks/01_BD2_QSAR_Modeling.ipynb) is an **executed model-evaluation notebook** (not full retraining). It reports held-out metrics, 20 repeated scaffold splits, 30 y-randomization trials, 2,520 fold-derived AD predictions, and the full 1,764-drug FDA screening snapshot. Ten distinct figures are embedded with consistent publication styling (bold 15-pt axis labels, bold 12-pt ticks, 2-pt spines, no titles or grids). The notebook can export ten figures in 600-dpi PNG, PDF and SVG when executed.

Inputs are in `data/`, `validation/` and `screening/`. Structural similarity ≥0.30 is a permissive **screening inclusion threshold**, not calibrated prediction reliability.
