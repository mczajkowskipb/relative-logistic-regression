# RLR vs binary ENET/MCP: 11 cohorts x 5 outer folds — COMPLETE (2026-10-30)

**Science status: EXPLORATORY, NOT MANUSCRIPT-READY CLAIM OF SUPERIORITY TO ORIGINAL PURIST/DECAF.**

## Reproducible computation
Real [GitHub Actions run #41064335985](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41064335985) 12/12 PASS (11 cohort jobs + independent complete patient-level Python audit). Source commit for trained runs: 01d0c42cd7e1a2eba0ff74e4386918be30cd9556. Independent audit replayed 165 method/fold outputs from heldout patient probabilities; each patient appears exactly once per model within a cohort's five outer folds.

**Source limitation:** input is the public Schliep reconstructed CompCancer collection. Exact original paper processed ARFF SHA256 parity is **NOT VERIFIED**. This is a single 5-fold repetition (55 outer splits), not the paper's 470-fold schedule, and these are the same 11 cohort identities that have been explored previously, not independent prospective datasets.

## Controlled comparison
Identical training-only switchBox top1000 candidate ranking, strict gene-disjointness d2, K20 pairs in all three models:
- Original RLR20: frozen original RLR20 m10 ternary rank-gap encoding, glmnet alpha .5.
- Binary ENET20, PurIST-inspired: binary pair comparison and glmnet alpha .5.
- Binary MCP20, DeCAF-inspired: binary pair comparison and ncvreg MCP alpha .05.
These binary models are methodological adaptations **not full original PurIST/DeCAF pipelines or their disease-specific fitted models**. Original DeCAF also incorporates a cross-cohort gene filtering/selection step.

## 11 dataset-equal mean results (five folds per cohort)
| Arm | F1 ↑ | AUC ↑ | Logloss ↓ | Brier ↓ | ECE ↓ | Mean active coefficients |
|---|---:|---:|---:|---:|---:|---:|
| **RLR20_ORIGINAL** | **0.876180** | **0.961673** | **0.237944** | **0.062638** | **0.092092** | 16.818 |
| BINARY_ENET20_PURIST_INSPIRED | 0.823985 | 0.948487 | 0.290442 | 0.077384 | 0.110474 | 16.745 |
| BINARY_MCP20_DECAF_INSPIRED | 0.813641 | 0.947843 | 0.291233 | 0.085455 | 0.134420 | 11.527 |

## Exact paired inference (2^11 sign assignments, dataset units)
RLR minus comparator. Negative logloss/Brier and positive AUC/F1 favor RLR. Exploratory BH over the eight listed tests:
| Comparator | Metric | Mean Δ | Exact two-sided p | BH exploratory q |
|---|---|---:|---:|---:|
| ENET | Logloss | -0.052498 | 0.034180 | **0.039062** |
| ENET | Brier | -0.014746 | 0.008789 | **0.023438** |
| ENET | AUC | +0.013187 | 0.078125 | 0.078125 |
| ENET | F1 | +0.052194 | 0.031250 | **0.039062** |
| MCP | Logloss | -0.053290 | 0.026367 | **0.039062** |
| MCP | Brier | -0.022817 | 0.001953 | **0.015625** |
| MCP | AUC | +0.013830 | 0.031250 | **0.039062** |
| MCP | F1 | +0.062539 | 0.007812 | **0.023438** |

Seven of eight exploratory FDR tests meet q<0.05. Differences are not confirmatory because the general data worlds and earlier one-split pilot have been analyzed; the eight-test family and model adaptations are method-specific and are not a test against the full clinical original algorithms. Correlated metrics and small cohort count warrant caution.

## All 11 dataset-level logloss estimates
| Cohort | RLR20 | Binary ENET20 | Binary MCP20 |
|---|---:|---:|---:|
| Alizadeh | 0.2365 | **0.2347** | 0.3786 |
| Armstrong | **0.0118** | 0.0179 | 0.0785 |
| Bittner | 0.4567 | **0.4509** | 0.6129 |
| Chen | **0.1699** | 0.1949 | 0.1896 |
| Golub | 0.1872 | 0.1581 | **0.1557** |
| Gordon | **0.0321** | 0.0344 | 0.0893 |
| Pomeroy | 0.5075 | 0.6437 | **0.4424** |
| Shipp | **0.1774** | 0.2683 | 0.2497 |
| Singh | **0.2580** | 0.3483 | 0.3236 |
| West | **0.4153** | 0.5901 | 0.4308 |
| Yeoh | **0.1650** | 0.2536 | 0.2524 |

## Science gate and further action
Report all positive and negative results. DO NOT modify the frozen manuscript or claim full superiority to PurIST or DeCAF from this exploratory benchmark. The owner offered an Ubuntu server containing the publication datasets. Verify the **exact original processed ARFF SHA256** with the read-only script at scripts/rlr_original_compcancer_sha_preflight_v2.py, then run an exact-original-data counterpart with a separately declared analysis scope.

## Retained artifact evidence
Actions run #41064335985 audit artifact: rlr-binary-outer5-all11-audit-41064335985 (all11_summary.csv, all11_paired_exact_inference.csv, all11_dataset_equal_metrics.csv, all11_165_outer_fold_metrics.csv, all11_outer_patient_predictions.csv, audit.log). Artifact retention 30 days. Frozen protocol: docs/RLR_BINARY_OUTER5_V2_LOCKED_PROTOCOL.md. Source runner: R/experimental/rlr_binary_outer5_v2.R; auditor: scripts/audit_rlr_binary_outer5_v2.py.
