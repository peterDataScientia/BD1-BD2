# BRD4-BD1 — curated biochemical IC50 data

## Current model-ready dataset (2026-11-18)

- **[BD1 model-ready data](data/BRD4_BD1_QSAR_model_ready_2816.csv)**: 2,816 unique standardized molecular structures, pIC50, back-converted IC50 nM, and assay/document provenance.
- **[Compound aggregation audit](audit/BD1_compound_aggregation_audit.csv)**: 2,827 distinct structures before cross-assay disagreement exclusions, including 11 compounds marked for review.
- **[Assay audit](audit/BD1_assay_audit.csv)**: 685 assay-level domain/format classification decisions, of which 228 were automatically accepted for BD1 biochemical IC50.
- **[Summary](audit/summary.json)**: source records and exclusion counts.
- **[Reproducible preprocessing script](notebooks/01_BD1_ChEMBL_Preprocessing.ipynb)**: reads the ChEMBL-derived exact-IC50 export and writes all audit/curation CSVs.

### Source and filtering

The accessible source was `BRD4_IC50_curated_preliminary.csv` (10,919 exact IC50/nM rows from the earlier BRD4 ChEMBL extraction). The **original raw ZIP with 58,225 activity records was not available in this workflow**. This means the published files are **curated from the available ChEMBL-derived extract**, not a fresh re-extraction of the original ZIP.

Conservative rules: explicit BD1-only biochemical assay descriptions; review mixed BD1/BD2 descriptions, mutant assays, DSF/thermal-shift, cellular, or otherwise ambiguous formats; exact positive IC50 nM; RDKit Cleanup/FragmentParent SMILES; within-assay median pIC50 followed by across-assay median; flag >1 pIC50 unit disagreement between assays. Automated assay assignments remain subject to manual source/construct review.

The 2,816 figure is **not** asserted to be identical to the earlier 2,908-count cohort: source completeness and conservative assay classification differ. No generated or hypothetical potency values have been added.

### Earlier pilot (kept separate)

`screening/` and `validation/BD1_DrugCentral_Similarity_Calibrated_Report.md` refer to an **older pilot** with 1,939 provisionally annotated compounds, 1,499 development, and 440 held out. These outputs must not be attributed to the newly curated 2,816-compound dataset. New QSAR model training and validation have not yet been performed on the 2,816 structures.
