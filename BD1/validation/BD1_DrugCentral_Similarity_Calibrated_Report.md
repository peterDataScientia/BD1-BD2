# BRD4-BD1 Morgan SVR: similarity-dependent DrugCentral screening

This is the recalibrated BD1-specific workflow corresponding to the prior SelleckChem L1300 strategy. It does **not** copy calibration errors from the earlier mixed-domain model. It uses the **440 held-out BD1 compounds** and **1,499 development compounds** to associate held-out absolute errors with nearest-neighbor Morgan Tanimoto similarity.

- DrugCentral entries screened: 4,091
- Entries matched to DrugCentral FDA/EMA/PMDA historical approval lists: 1,963
- Eligible for small-molecule screening (MW <=900 Da and no excluded exact model overlap): 1,873
- Historical approval-list status does not establish current authorization/marketing status.

## Calibration

Fixed similarity bands: [0,.2), [.2,.3), [.3,.4), [.4,.5), [.5,.6), [.6,.7), [.7,.8), [.8,1]. Categories based on held-out counts in each band: >=30 relatively supported; 10–29 limited; 1–9 low support; none not calibrated. These descriptive categories indicate the density of evidence supporting local error estimates, not proof of an individual drug's activity.

| similarity_band   |   n |        MAE |   p90_abs_error |   median_abs_error | support_category     |
|:------------------|----:|-----------:|----------------:|-------------------:|:---------------------|
| [0,0.2)           |   0 | nan        |      nan        |         nan        | not_calibrated       |
| [0.2,0.3)         |   2 |   0.559138 |        0.817283 |           0.559138 | low_support          |
| [0.3,0.4)         |   7 |   1.10877  |        2.21392  |           0.937498 | low_support          |
| [0.4,0.5)         |  26 |   0.553883 |        1.45343  |           0.272689 | limited_support      |
| [0.5,0.6)         |  69 |   0.36083  |        0.774758 |           0.261613 | relatively_supported |
| [0.6,0.7)         | 115 |   0.368011 |        0.908416 |           0.231691 | relatively_supported |
| [0.7,0.8)         | 125 |   0.441758 |        0.937906 |           0.317643 | relatively_supported |
| [0.8,1]           |  96 |   0.356998 |        0.703535 |           0.224945 | relatively_supported |

## DrugCentral eligible outcomes

{'low_support': 1166, 'not_calibrated': 694, 'limited_support': 11, 'relatively_supported': 2}

## Candidate shortlist (not experimentally validated)

| DrugCentral_name   |   predicted_BRD4_BD1_pIC50 |   predicted_BRD4_BD1_IC50_nM |   Morgan_nearest_training_Tanimoto |   heldout_bin_MAE |   heldout_bin_p90_abs_error | support_category     | literature_evidence_flag                                                                              |
|:-------------------|---------------------------:|-----------------------------:|-----------------------------------:|------------------:|----------------------------:|:---------------------|:------------------------------------------------------------------------------------------------------|
| alprazolam         |                    6.05825 |                      874.487 |                            0.73913 |          0.441758 |                    0.937906 | relatively_supported | Published BRD4-BD1 binding evidence; not novel; Kd not IC50                                           |
| triazolam          |                    5.86658 |                     1359.63  |                            0.65306 |          0.368011 |                    0.908416 | relatively_supported | Published comparative BET studies report weak/negligible binding; prediction needs assay confirmation |
| midazolam          |                    6.1743  |                      669.418 |                            0.4     |          0.553883 |                    1.45343  | limited_support      | Previously documented BRD4/BET binding; not novel                                                     |
| brigatinib         |                    6.09011 |                      812.629 |                            0.49398 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| estazolam          |                    5.72398 |                     1888.06  |                            0.43636 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| diazepam           |                    5.32979 |                     4679.64  |                            0.4     |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| donepezil          |                    5.25976 |                     5498.5   |                            0.41429 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| oxyquinoline       |                    5.17797 |                     6637.94  |                            0.41667 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| edaravone          |                    4.90381 |                    12479.4   |                            0.46512 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| metofoline         |                    4.70892 |                    19546.9   |                            0.41509 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| phensuximide       |                    4.45058 |                    35433.7   |                            0.4     |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| methyl salicylate  |                    4.26373 |                    54484.1   |                            0.44444 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |
| tetrabenazine      |                    4.17702 |                    66524.4   |                            0.43396 |          0.553883 |                    1.45343  | limited_support      | Independent BRD4 literature review and regulatory check required                                      |

The p90 absolute errors are retrospective held-out prediction error summaries and **not** formal 90% coverage intervals. BRD4 literature validation and compound regulatory status are pending for most shortlisted compounds. The BD1 1,939-compound cohort remains provisionally annotated. Low-support results should not be used for confirmatory potency claims.