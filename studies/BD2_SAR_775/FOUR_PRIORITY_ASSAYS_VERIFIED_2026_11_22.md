# Verified original ChEMBL assay sources — 22 November 2026

Four ChEMBL assays account for **44/63 (69.8%)** same-assay candidate cliffs. The assay and document API records were read directly.

| Assay | Candidate cliffs | Source | Original research focus | Method |
|---|---:|---|---|---|
| CHEMBL6086970 | 17 | Li et al. 2024; DOI 10.1021/acs.jmedchem.4c02516 | First BRD4-BD2-selective inhibitors (XY221) | BRD4 BD2 residues 348–464; GST/His tagged; TR-FRET, 1 h |
| CHEMBL5529441 | 10 | Li et al. 2024; DOI 10.1021/acs.jmedchem.3c01837 | BD1-selective BET inhibitors, BD2 measured as comparator | His6 BRD4 BD2; TR-FRET, 15 min |
| CHEMBL5376137 | 9 | Bradley et al. 2023; DOI 10.1021/acs.jmedchem.3c00906 | Domain-selective BET BD1 chemical probe, BD2 comparator | BRD4 BD2 TR-FRET; origin unspecified |
| CHEMBL4733660 | 8 | Chen et al. 2021; DOI 10.1016/j.ejmech.2020.112868 | Dual HDAC/BRD4 inhibitors | Histone H4 acetyl peptide, His-tagged BD2; fluorescence-based |

Documents respectively: CHEMBL6083288, CHEMBL5528036, CHEMBL5375945, CHEMBL4732071. A BD2 endpoint measured for BD1-selective development is **not** evidence that the compounds were optimized for BD2.

## Eight original records independently confirmed

| Assay | Compound | Original ChEMBL activity ID | IC50 nM |
|---|---|---:|---:|
| CHEMBL6086970 | CHEMBL6145355 | 29105548 | 148.01 |
| CHEMBL6086970 | CHEMBL6147019 | 29105531 | 3.15 |
| CHEMBL5529441 | CHEMBL5569832 | 25925847 | 15.7 |
| CHEMBL5529441 | CHEMBL5573825 | 25925846 | 714.4 |
| CHEMBL5376137 | CHEMBL5416784 | 25670556 | 501.19 |
| CHEMBL5376137 | CHEMBL5437203 | 25670572 | 19952.62 |
| CHEMBL4733660 | CHEMBL4777392 | 22908436 | 3.8 |
| CHEMBL4733660 | CHEMBL4778531 | 22908434 | 103.0 |

All are IC50 nM with '=' relation and match the values in the original curated 775-row dataset. These confirm **four representative candidate cliffs**, not all 63 pair-level values. Three are single-cut MMPs; CHEMBL5416784/CHEMBL5437203 is not.

Original API:
- Assay: https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL6086970.json (replace ID for other assays)
- Document: https://www.ebi.ac.uk/chembl/api/data/document/CHEMBL6083288.json (replace ID)
- Activity: https://www.ebi.ac.uk/chembl/api/data/activity.json?assay_chembl_id=CHEMBL6086970&molecule_chembl_id=CHEMBL6145355

## Conclusion and limits

Novelty of any proposed chemical transformation is **not confirmed**: some source publications expressly report SAR themselves. This is an assay- and publication-aware retrospective landscape, not proof of new SAR mechanisms, reproducibility, selectivity across BD1/BD2, or clinical efficacy. Most dataset molecules have only one measurement.

Research candidate claim: chemical-series and assay clustering explain why 63 candidate pairs cannot be presented as 63 independent findings. Prior publication and assay-aware filtration are essential to meaningful SAR conclusions.
