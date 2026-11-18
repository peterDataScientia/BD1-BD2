# BRD4-BD1 Morgan ECFP4 Screening of DrugCentral Approval-List Compounds

**Execution:** November 5, 2026. **Model:** Saved RBF-SVR fitted to the 1,499-development-compound fraction of the 1,939-compound *provisional BD1* cohort, 2,048-bit Morgan ECFP4 (radius 2). This model was chosen by 5-fold scaffold-grouped development CV (MAE 0.4314 pIC50); its independent 440-compound test MAE was 0.4091 pIC50, without subsequent tuning. This is an exploratory baseline, not experimentally validated and not independently revalidated on approved drugs.

## Sampling

- DrugCentral original structure dataset has **4,099** entries, of which the prior screening eligibility file excluded **8** with BRD4 training/cohort overlap (including nitroxoline and fedratinib). The current screen deliberately reuses the prior **4,091 unseen structures**, rather than counting these 8 a second time.
- **1,963** of the 4,091 match DrugCentral's FDA, EMA or PMDA public approval-list snapshots, including **1,853** in the FDA file. **1,873** meet the prespecified exploratory small-molecule cutoff of molecular weight ≤900 Da.
- Regulatory approvals are linked by exact DrugCentral integer **ID**, not by approximate drug-name matching. This is based on publicly available static lists with 2024-04-10 timestamps; **current marketing authorization or withdrawal has not been verified**.
- The full audit includes all 4,091 predictions, with separate files for the official-list subset, MW-eligible subset, training-AD subset, and sensitivity subset.

## Applicability domain

Nearest-neighbor Morgan Tanimoto similarity was computed against the **1,499 training compounds**, not the 440 untouched holdout compounds. Training-only leave-one-out nearest-neighbor similarity gives the fifth-percentile threshold of **0.5342**. Only **2** approved, eligible molecules reach that threshold; **13** reach a sensitivity threshold of 0.40 and **2** reach 0.50. The 0.40 threshold is only a sensitivity analysis; it is *not* a validated relaxed applicability domain.

## Two compounds within training-derived similarity domain

| Drug | Predicted pIC50 | Predicted IC50 (nM) | Max Tanimoto | Literature context |
|---|---:|---:|---:|---|
| Alprazolam | 6.058 | 874.5 | 0.739 | **Known BRD4-BD1 binder**: published ITC Kd ≈2.5 µM and BRD4-BD1 cocrystal; not novel. |
| Triazolam | 5.867 | 1359.6 | 0.653 | Published BET benzodiazepine SAR reports **lack of appreciable BET binding**; qualitative model–experiment disagreement. |

Alprazolam's binding **Kd** is not an IC50 and cannot be used as an exact numeric reference to score the predicted IC50. Likewise triazolam's reported negative BET result may have assay and concentration limitations; it is a counterexample requiring closer investigation, not a measured IC50 for regression scoring.

### Highest model-predicted potency among approved-list compounds

The five numerically highest model predictions are **venetoclax, rifaximin, rifabutin, telithromycin, and rifapentine** (estimated IC50 111–150 nM), but all five have **Morgan similarity <0.24**, far outside the threshold of 0.534. They are **not sufficiently supported to present as plausible confirmed BRD4 inhibitors or strong experimental candidates** on the basis of the present QSAR model.

## Scientific interpretation

With this first BD1-specific pilot model, the approved-drug library lies overwhelmingly outside the supported chemical domain. Its two inside-domain compounds are not validated novel repurposing leads: one is a literature-confirmed binder and the other has contrary experimental BET evidence. For drug prioritization, this implies the need to improve endpoint adjudication, representation/model selection, and chemical-space coverage before strong prospective claims. Predicted values have not been confirmed by BRD4 biochemical assays or target engagement.

## Key literature and data sources

1. DrugCentral official downloads: https://www.drugcentral.org/download
2. DrugCentral official FDA list: https://www.drugcentral.org/static/FDA_Approved.csv
3. DrugCentral official EMA list: https://www.drugcentral.org/static/EMA_Approved.csv
4. DrugCentral official PMDA list: https://www.drugcentral.org/static/PMDA_Approved.csv
5. Filippakopoulos et al. (2012), *Benzodiazepines and benzotriazepines as protein interaction inhibitors targeting bromodomains of the BET family*: https://pmc.ncbi.nlm.nih.gov/articles/PMC4748212/
6. RCSB BRD4-BD1–alprazolam complex (PDB 3U5J): https://www.rcsb.org/structure/3U5J
7. Privileged diazepines review describing triazolam and estazolam lacking BET binding: https://pmc.ncbi.nlm.nih.gov/articles/PMC4035449/

## Reproducibility and limitations

The script `screen_approved.py` regenerates scores and training-based similarity checks from the exact saved model and structure files. Approval-list membership is recorded using compact bitsets of DrugCentral official ID lists to permit offline screening. The approved-drug model uses fixed exploratory hyperparameters, and the 1,939-molecule BD1 annotations are still source-description-supported rather than all publication-verified. Approval snapshots may not reflect 2026 changes. **No experimental validation is claimed.**