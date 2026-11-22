# Literature and assay-source verification — 22 November 2026

The curated 775-compound source was audited for shared-assay activity cliffs. Of 78 candidate cliffs (Tanimoto >=0.8, >=10-fold IC50 difference), 63 shared a ChEMBL assay ID. All 63 share a document ID within each pair; 46 were also single-cut matched molecular pairs under the operational fragmentation rule. These are **candidate comparisons**, not 63 independent novel SAR discoveries.

## First independently cross-referenced original source

Assay **CHEMBL5218148** is named ChEMBL_2263137 in BindingDB and is described as an **inhibition of BRD4 BD2 by TR-FRET assay** (origin unspecified). The corresponding source study is Humphreys et al., *Journal of Medicinal Chemistry* (2022), *Identification and Optimization of a Ligand-Efficient Benzoazepinone Bromodomain and Extra Terminal (BET) Family Acetyl-Lysine Mimetic into the Oral Candidate Quality Molecule I-BET432*.

Source: https://www.bindingdb.org/rwd/jsp/dbsearch/assay.jsp?assayid=1&entryid=50018092

From the curated records, **CHEMBL5219280** and **CHEMBL5220623** have same-assay IC50 values **63095.73 nM** and **630.96 nM** respectively (~100-fold difference). The isomeric SMILES have identical achiral connectivity but opposite stereochemical assignment at one center; chirality-agnostic ECFP4 Tanimoto=1.0. This is a **stereochemical candidate activity cliff** and NOT a single-cut substituent replacement. Original compound-level measurement and original paper's explanation are still to be independently checked.

The 2022 Humphreys article also explicitly reports an activity cliff for an N-linked imidazole (37) versus a C-linked counterpart (34) in its optimization discussion. Therefore, **activity-cliff observations in related BET chemistry predate this project**. Do not claim that this is the first such observation.

## Pending verification

The description and exact protocol for assays CHEMBL5369309, CHEMBL6086970 and other prioritized assays, plus pair-level original experimental results, are **not verified**. Existing original data are mostly single measurements (771/775), preventing independent replication-based variance estimates. Same-document and same-assay overlap do not certify reliability or novelty.

See the ChatGPT-delivered `BRD4_BD2_Verification_Updated_Package.zip` for the complete 63-pair source-prioritization CSV, figures, notebooks, and revised manuscript.

**Status:** exploratory manuscript; not submission-ready until source-paper pair mapping and independent novelty audit are complete.
