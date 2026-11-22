# BRD4-BD2: Assay-aware SAR and activity cliffs (775 compounds)

Exploratory independent study of **775 unique BRD4-BD2 molecular structures** from the curated ChEMBL-derived `BRD4_BD2_STRICT_MODELING_775.csv`. All reported statistics were calculated from the curated file on 22 November 2026.

## Results
- 775 valid unique structures; 390 Murcko scaffolds (284 singleton scaffolds).
- 79 distinct assay identifiers; 771 molecules have only one potency measurement.
- ECFP4 Tanimoto >=0.80: 665 similar pairs; 78 candidate >=10-fold activity cliffs.
- 63 candidate cliffs share an assay identifier.
- RDKit single-cut matched molecular pairs: 3,108 pairs; 2,934 share an assay ID.
- Stringent MMP subset (shared assay, Tanimoto >=0.80, >=10-fold difference): **46 pairs involving 66 compounds from 12 assay IDs**.

## Criteria
Morgan radius 2 / 2048 bits; threshold Tanimoto>=0.80; |delta pIC50|>=1; RDKit single-cut matched fragment, <=12 variable heavy atoms, constant fragment >=2 times variable fragment. The paired same-assay designation means *assay ID overlap*, not experimental replication.

## Reproduction
Use the paired script and notebook in this folder. Place the source CSV alongside the notebook or in `data/BRD4_BD2_STRICT_MODELING_775.csv`. The full reproducibility bundle, tables, plots and working manuscript were provided in the ChatGPT conversation as `BRD4_BD2_Complete_Study_Package.zip`.

## Scope
Descriptive SAR candidate discovery. **No docking, MD, BD1/BD2 selectivity calculation, experimental confirmation, disease model or new binding measurement.** Before publication, verify original ChEMBL assay descriptions and compare each proposed SAR rule against primary medicinal-chemistry papers. This study does not establish clinically useful or newly discovered inhibitors.

## Added statistical robustness results

Among 635 pairs in single-assay strata with ECFP4 Tanimoto >=0.80, the observed number of >=10-fold potency differences was 63. In 2,000 potency permutations within assay groups, the null mean was 157.18 (SD 11.22, 95% empirical interval 136–179; lower-tail empirical p=0.00050). Similarity thresholds 0.70, 0.75, 0.80, 0.85, 0.90 yielded respective cliff fractions 14.97%, 12.98%, 9.92%, 7.26%, 6.45%. These figures describe structural association, **not independent experimental validation**.

The complete downloadable study package includes `robustness.py`, sensitivity tables, permutation distribution, five figures and a manuscript draft. The GitHub notebook currently reproduces baseline descriptive findings from the CSV in this repository.
