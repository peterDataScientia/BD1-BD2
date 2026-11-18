# BRD4 BD1/BD2 research data index

The repository separates domain-specific files into `BD1/` and `BD2/`, with shared documentation in `Shared/`.

## Uploaded BRD4-BD2 materials

- `BD2/data/BRD4_BD2_STRICT_MODELING_775.csv`: 775 molecular structures with curated biochemical IC50-derived pIC50 and assay provenance. This is a conservatively filtered modeling subset and is not equivalent to full manual source-paper auditing.
- `BD2/validation/heldout_predictions_with_AD.csv`: predictions for 187 scaffold-separated held-out compounds, residuals and nearest training ECFP4 Tanimoto similarities.
- `BD2/validation/development_AD_threshold_performance.csv`: development fold coverage and accuracy at tested similarity thresholds.
- `BD2/validation/similarity_bin_validation_errors.csv`: bin-specific validation-derived errors.
- `BD2/screening/BD2_FDA_AD03_unseen.csv`: 33 FDA-library candidates satisfying ECFP4 nearest-neighbor Tanimoto similarity >= 0.30.
- `BD2/notebooks/01_BD2_QSAR_Modeling.ipynb`: model training and evaluation pipeline.
- `BD2/notebooks/02_BD2_Validation_Derived_Errors.ipynb`: similarity-dependent error analysis.

## Methodological cautions

Screening predictions are **computational estimates**, not measured BRD4-BD2 inhibition. The 0.30 ECFP4 similarity cutoff is a selected structural inclusion criterion, not an empirically established individual-prediction reliability guarantee. Validation errors should be interpreted with scaffold splits and dataset coverage in mind.

The complete 1,764-drug screening table, complete assay audit, trained model binary, and the BD1-specific datasets are not included in this initial upload. They should not be inferred from the files above.
