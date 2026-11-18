# BRD4-BD1 data and preliminary screening

This folder currently contains **pilot BD1-specific QSAR screening outputs**, not a verified final 2,908-compound BD1 modeling dataset.

## Files

- `screening/Inside_AD.csv` — 2 approved-list compounds meeting the pilot model's training-derived threshold (LOO 5th percentile, 0.5342).
- `screening/Near_AD_0p40.csv` — 13 approved-list compounds with Tanimoto similarity >= 0.40.
- `screening/Top_Potency_20.csv` — 20 numerically highest potency predictions, which should **not** be treated as in-domain inhibitor discoveries.
- `screening/screening_report_pilot.md` — methodology, model limitations and literature interpretation.
- `validation/BD1_DrugCentral_Similarity_Calibrated_Report.md` — similarity-dependent validation diagnostics for the earlier BD1-specific cohort.

## Provenance and important distinctions

The uploaded pilot materials correspond to a **provisionally BD1-annotated 1,939-compound** cohort, split into 1,499 development and 440 held-out compounds, and a 2024 historical DrugCentral approvals snapshot. They are not yet the later reported 2,908-compound biochemical BD1 curation. The latter raw model-ready dataset was **not available for verified upload** and is not represented by these pilot files.

No measured BD1 inhibitory potency can be inferred from predicted IC50 values. The pilot training-derived AD cutoff of 0.5342 should **not** be silently replaced with the BD2 cutoff of 0.30, nor treated as optimized for a later final BD1 model. Screening records from different model versions must remain distinct.

The full BD1 raw dataset, complete DrugCentral screening workbook, serialized model and final repeat-scaffold predictions are **not yet uploaded**.
