# Joint BRD4-BD1/BD2 comparative screening and uncertainty

**Executed:** 18 November 2026.

## Materials
- [Executed joint uncertainty notebook](BD1_BD2_Joint_Activity_Uncertainty_Executed.ipynb)
- [29 FDA-candidate uncertainty audit](joint_screen_uncertainty_audit_29.csv)
- [Held-out error by structural-similarity range](heldout_error_by_AD_similarity.csv)

## Main findings
Of 1,764 historical DrugCentral 2024 FDA-listed compounds screened by both independent models, 29 unseen compounds passed ECFP4 maximum training similarity ≥0.30 for both domains; only two passed ≥0.40 in both models; none passed ≥0.50. Five had predicted IC50 ≤1 µM in both domains. Three showed ≥10-fold **predicted** BD2 preference; none showed ≥10-fold **predicted** BD1 preference.

Low-similarity (<0.40) independent hold-out mean absolute errors were 0.964 pIC50 (BD1, n=24) and 1.587 pIC50 (BD2, n=12). There are too few low-similarity observations for well-calibrated selectivity probabilities. Resampled residual-difference percentiles in the CSV are **illustrative sensitivity calculations, not confidence intervals or validated probabilities of real binding**; they assume independent residuals from unmatched datasets.

### Literature verification note
Alprazolam is **already a published BRD4-BD1 binder** with calorimetry and co-crystal evidence (PDB [3U5J](https://www.rcsb.org/structure/3U5J); [Filippakopoulos et al.](https://pubmed.ncbi.nlm.nih.gov/22137933/)). It is not a previously unknown BRD4 hit. The predicted BRD2/BD1 activity ratio is not a verified selectivity measurement.

Biological context: [Gilan et al., Science (2020)](https://pubmed.ncbi.nlm.nih.gov/32193360/) and [Faivre et al., Nature (2020)](https://www.nature.com/articles/s41586-020-1930-8/).

**Status:** computational hypothesis prioritization only. Matched biochemical BD1/BD2 measurements and orthogonal binding tests remain necessary. Full source data for rerunning the notebook are included in the accompanying dissertation analysis package, not directly in the notebook file.
