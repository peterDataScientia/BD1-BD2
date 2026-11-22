# Completed direct ChEMBL source audit — 22 November 2026

## What was checked

- **775** curated BRD4-BD2 chemical structures, **390** Murcko scaffolds, **79** assay IDs overall.
- **665** pairs at ECFP4 (radius 2, 2048 bits) Tanimoto >=0.80; **78** candidate >=10-fold potency cliffs, of which **63** share a ChEMBL assay ID.
- These 63 candidate pairs involve **85 unique molecules** and arise from **15** ChEMBL assay IDs / primary publication IDs. Four assays contribute **44/63 (69.8%)**. They cannot be treated as 63 independent biological findings.
- Independently retrieved **all 15 assay descriptions** and associated **15 primary-source document records** from the original ChEMBL API; the complete assay/document mapping is [ALL_15_SOURCE_PAPERS_VERIFIED.csv](ALL_15_SOURCE_PAPERS_VERIFIED.csv).
- Independently retrieved compound-level ChEMBL IC50 records for the largest-fold candidate pair in every assay: **15 pairs / 30 compound-assay requests / 34 activity records**. The activity identifiers and numerical values are in [ALL_15_PRIORITY_PAIRS_34_ORIGINAL_ACTIVITIES.csv](ALL_15_PRIORITY_PAIRS_34_ORIGINAL_ACTIVITIES.csv).
- All 15 priority pairs have source records equal to their selected curated IC50 values. **Fourteen** pairs remain >=10-fold across all returned original-record comparisons. **One pair is not robust to record selection**.
- **CHEMBL3584940**: CHEMBL3581650 original IC50 entries are **478 and 1817 nM**; CHEMBL3581652 entries **8322 and 73970 nM**. Depending on comparison, difference spans **4.58–154.75-fold**. The 10-fold cliff assignment therefore depends on which experimental record is chosen.
- **CHEMBL4185405** also contains multiple source entries: CHEMBL4209937 is **146 / 436 nM**, CHEMBL4210754 is **3.8 / 12.7 nM**. The differences span **11.50–114.74-fold**, remaining >=10-fold in all pairwise comparisons.
- The other 13 priority pairs returned one activity per queried compound in the results inspected.

## Why source context matters

Not all source studies optimized selective BRD4 BD2 inhibitors. Some source studies primarily designed BET-BD1-selective probes, bivalent inhibitors, pan-BET inhibitors, HDAC/BRD4 dual inhibitors, BET degraders, or compounds for non-oncology indications; BRD4 BD2 potency was a comparison endpoint in some of them. The original assay formats include TR-FRET, HTRF, AlphaScreen, fluorescence polarization, fluorescence anisotropy and other fluorescence assays. These are not one standardized experimental protocol.

## Scientific interpretation and scope

The result supports **assay-aware and publication-aware retrospective SAR mapping**, plus a concrete example in which a source-record choice changes the activity-cliff interpretation. It does **not** establish new experimentally confirmed inhibitors, BD1-vs-BD2 selectivity, independent measurement uncertainty, clinical activity, or novelty for a specific chemical substitution. A primary medicinal chemistry source, Humphreys et al. (2022, DOI 10.1021/acs.jmedchem.2c01102), already discusses a BET inhibitor activity cliff.

**Source-level audit complete for all 15 assays; exact activity checks complete for 15 representative pairs out of 63, not for every one of the 63 pairs.** ChEMBL source activity IDs are available in the CSV to reproduce the targeted checks. Further original-paper structure-number matching and complete pair-level evidence would be needed before claiming novel medicinal chemistry rules. The main 775-compound SAR methods and dataset are in this repository.

## Original source APIs

- Assay record: https://www.ebi.ac.uk/chembl/api/data/assay/CHEMBL6086970.json
- Publication record: https://www.ebi.ac.uk/chembl/api/data/document/CHEMBL6083288.json
- Specific activity: https://www.ebi.ac.uk/chembl/api/data/activity/29105548.json

Replace these identifiers with those listed in the attached source-audit CSV files.
