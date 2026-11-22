# Primary paper and dependency audit — 22 November 2026

## Pair non-independence
From the original 775-compound ChEMBL-derived dataset, 63 ECFP4-similar (>=0.80), >=10-fold, shared-assay candidate activity cliffs comprise 85 unique compounds from **15 assay IDs / 15 source-document IDs**. Four assays account for **44/63 (69.8%)** pairs: CHEMBL6086970 (17), CHEMBL5529441 (10), CHEMBL5376137 (9), CHEMBL4733660 (8). Pair counts are **not independent replications**.

Forty-six pairs also meet a single-cut MMP operational definition. One separate stereochemistry-only pair (CHEMBL5219280 / CHEMBL5220623) has 63095.73 versus 630.96 nM recorded potencies (~100x). Chirality-insensitive ECFP4 Tanimoto=1.000; chirality-aware Morgan Tanimoto=0.939. Do not equate ECFP4 identity with stereochemical identity.

## Confirmed primary literature
Humphreys et al., 2022, *J Med Chem*, DOI: 10.1021/acs.jmedchem.2c01102, **I-BET432 medicinal chemistry**, has been read from the accepted author manuscript (https://strathprints.strath.ac.uk/83588/1/Humphreys_etal_JMC_Identification_and_optimization_of_a_ligand_efficient_benzoazepinone.pdf). Table 3 and accompanying discussion report an N-linked imidazole (rac-37) vs C-linked imidazole (34) **activity cliff**: BD1/BD2 pIC50 = 4.4/4.6 versus 8.0/8.0. This is ~2512-fold BD2 difference, not a novel observation by our study. The full-text article also provides original experimental methodological context.

CHEMBL5218148 is associated with the BRD4 BD2 TR-FRET assay in BindingDB (https://www.bindingdb.org/rwd/jsp/dbsearch/assay.jsp?assayid=1&entryid=50018092). Compound-level mapping between the 100x stereochemical pair and paper's numbered structures is **not verified**; do not conflate them.

## Manuscript implication
Our contribution is under investigation: assay-aware *systematic* structure–potency mapping, chemical-series concentration, and chirality-dependent activity-cliff interpretation. The work is exploratory and not submission-ready until high-priority experimental records and prior SAR from source papers are reconciled. We cannot infer independent experimental uncertainty from the 771/775 molecules with a single measurement.
