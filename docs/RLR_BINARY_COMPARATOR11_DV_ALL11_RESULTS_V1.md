# RLR20 vs binary ENET and DeCAF-inspired MCP: first 11 D/V pilot — completed, exploratory

Date: 2026-10-30. **GitHub Actions [#41062249487](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41062249487): 12/12 jobs SUCCESS**, comprising 11 real R cohort runs plus independent Python audit. The auditor replayed logloss and Brier from the heldout per-patient probabilities and enumerated all 2^11 sign assignments.

Original frozen RLR helper, original manuscript lock, and original outer E remain unchanged. **These are public Schliep reconstructed CompCancer matrices, NOT verified SHA256-identical to original paper's processed ARFF.** This is an exploratory pilot on first outer T only split D/V, **not** a full 470-fold replication and **not** an independently held-out confirmatory comparison. The two binary models are **PurIST-/DeCAF-inspired controlled adaptations**, not the respective full published clinical algorithms.

## Design frozen BEFORE this pilot outcome
11 original-paper-named public cohorts. Same D-only switchBox supervised pair ranking top1000, strict gene-disjoint shortlist d2, K20, identical pairs in all arms. D-only CV. Original RLR uses rank margin m10, glmnet elastic net alpha=.5; binary ENET uses hard binary order, same glmnet alpha=.5; binary MCP uses hard binary order, ncvreg MCP alpha=.05 reflecting original DeCAF training logic.

## Dataset-equal means

| Arm | Positive F1 ↑ | AUC ↑ | Logloss ↓ | Brier ↓ | ECE ↓ | Nonzero pairs |
|---|---:|---:|---:|---:|---:|---:|
| **RLR20_ORIGINAL** | **0.908225** | **0.944545** | **0.260502** | **0.063266** | **0.105673** | 17.545 |
| BINARY_ENET20_PURIST_INSPIRED | 0.880313 | 0.937228 | 0.293891 | 0.071901 | 0.108822 | 17.273 |
| BINARY_MCP20_DECAF_INSPIRED | 0.761775 | 0.900056 | 0.367416 | 0.111287 | 0.146178 | 8.000 |

F1 and logloss per dataset (unweighted dataset units):

| Dataset | RLR LL | Binary ENET LL | Binary MCP LL | RLR F1 | ENET F1 | MCP F1 |
|---|---:|---:|---:|---:|---:|---:|
| Alizadeh | 0.0617 | 0.0535 | 0.5616 | 1.0000 | 1.0000 | 0.8889 |
| Armstrong | 0.0015 | 0.0032 | 0.1133 | 1.0000 | 1.0000 | 0.9474 |
| Bittner | 0.2258 | 0.2425 | 0.9766 | 1.0000 | 1.0000 | 0.6667 |
| Chen | 0.1405 | 0.1599 | 0.1354 | 0.9333 | 0.9032 | 0.9032 |
| Golub | 0.0258 | 0.0122 | 0.1469 | 1.0000 | 1.0000 | 0.9091 |
| Gordon | 0.0167 | 0.0246 | 0.0779 | 1.0000 | 0.9231 | 0.8000 |
| Pomeroy | 1.0076 | 1.0726 | 0.7782 | 0.4000 | 0.4000 | 0.0000 |
| Shipp | 0.1489 | 0.0843 | 0.2246 | 0.8571 | 0.8571 | 0.7500 |
| Singh | 0.1254 | 0.1512 | 0.1832 | 1.0000 | 1.0000 | 1.0000 |
| West | 1.1000 | 1.1821 | 0.6589 | 0.8000 | 0.8000 | 0.8000 |
| Yeoh | 0.0116 | 0.2468 | 0.1849 | 1.0000 | 0.8000 | 0.7143 |

## Exact paired 11-dataset two-sided sign-flip analysis

RLR minus binary comparator; negative is favorable for logloss/Brier, positive for F1/AUC.

| Comparator | Metric | Mean difference | Exact raw p | BH q over all eight exploratory tests |
|---|---|---:|---:|---:|
| Binary ENET | Logloss | -0.033389 | 0.154297 | 0.308594 |
| Binary ENET | Brier | -0.008635 | 0.283203 | 0.333705 |
| Binary ENET | AUC | +0.007317 | 0.500000 | 0.500000 |
| Binary ENET | F1 | +0.027912 | 0.250000 | 0.333705 |
| Binary MCP | Logloss | -0.106914 | 0.291992 | 0.333705 |
| Binary MCP | Brier | -0.048021 | 0.007812 | **0.031250** |
| Binary MCP | AUC | +0.044489 | 0.062500 | 0.166667 |
| Binary MCP | F1 | +0.146450 | 0.003906 | **0.031250** |

The exploratory BH adjustment is for eight pilot comparisons (two comparators × four metrics). Both MCP Brier and F1 remain under q=.05 within this pilot family. Correlated metrics and 11 small cohorts impose limitations. Do not conflate these p/q values with original publication Holm/BH comparisons against kTSP.

## Interpretation / science gate

- **RLR better numerically in all four primary descriptive metrics vs both binary controls**. RLR beats binary ENET in logloss on 8/11 cohorts, MCP in logloss on 8/11.
- Against the more relevant binary elastic-net comparator, logloss advantage **not statistically supported** (two-sided p=.1543); therefore **the added contribution of rank-margin encoding has NOT been established**.
- MCP result looks substantially weaker on this particular shared K20 selection, and the differences in Brier/F1 are exploratory. MCP hyperparameters, penalty and selection protocol differ from full original DeCAF, and this **cannot** be advertised as RLR superior to DeCAF.
- No changes to submission manuscript justified yet. Before any submission claim, verify original processed ARFF checksums and run a genuinely independent confirmatory comparison with correctly implemented previously published approaches. No selective exclusion of adverse cohorts.

## Provenance

- Frozen protocol: [RLR_BINARY_COMPARATOR11_DV_PROTOCOL_V1.md](RLR_BINARY_COMPARATOR11_DV_PROTOCOL_V1.md)
- R training runner: [rlr_binary_comparator_first11_DV_v1.R](../R/experimental/rlr_binary_comparator_first11_DV_v1.R)
- Independent auditor: [audit_rlr_binary_comparator11_v1.py](../scripts/audit_rlr_binary_comparator11_v1.py)
- Complete patient-level + cohort CSV audit is a GitHub Actions artifact in run [#41062249487](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41062249487), with 30-day retention.
