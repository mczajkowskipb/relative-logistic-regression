# RLR-GLVT v5 — complete 11-cohort D/V NEGATIVE science gate

**2026-10-29. Final verdict: TECHNICAL PASS, SCIENTIFIC NEGATIVE. No outer E or full470.**

## Execution and data integrity

- [GitHub Actions #40919146672](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40919146672): **12/12 PASS**, 11 source cohorts + independent audit of **55 model–cohort V patient score summaries**, on fixed commit **4e92695707c891624d85e6d69769732f08600c07**.
- [Independent v5 synthetic test #40918731712](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40918731712): PASS; exactly optimized finite tree verified against independently enumerated candidates, <=2 rank comparisons per patient, teacher-free prediction, strictly increasing patientwise rank invariance.
- Original publication outer first TRAIN T divided into disjoint D/V. All source selection, parent RLR20 and cross-fitted OOF predictions train on D only. E not read in study. Full D parent is frozen original RLR20, real switchBox kTSP5, unmodified RLR and manuscript.
- Public Schliep reconstructed CompCancer arrays with original cohort names; original ARFF SHA and probe mapping **NOT verified**; these datasets were seen repeatedly in prior research, therefore this is **exploratory**, not independent external validation.
- Three explicitly marked ridge surrogate teacher OOF folds for **Pomeroy only**; no surrogate replaces original full-D RLR20. All other cohort OOF teachers required zero fallbacks.

## Mechanism

A self-contained, tiny binary rank-pair tree. Each node compares one pair of genes; max depth 2, max internal nodes 3, maximum two actually tested pair relations per patient. Given identical D-only original switchBox TSP top20 candidate pairs, exact global optimization over eligible 0–3 internal-node structures minimizes dimensionless normalized genuine-label cross-entropy, plus the optional (w=0 or 0.25) **within-leaf OOF logit variance of original RLR20**, plus 0.005 per node and 0.005 mean D path length. Leaves use true D labels ONLY, smoothed by same rule in both arms. RLR never runs at prediction. It is not C-RET and has no exact teacher fidelity certificate.

[Preregistered objective and controls](RLR_GLVT_V5_GLOBAL_LOGIT_VARIANCE_TREE_PREREGISTRATION.md). Scientific novelty **not presumed**: globally optimal decision trees, surrogates and expected path complexity are prior art.

## Results — 11 equally weighted datasets, heldout D/V only

| Arm | Positive F1 ↑ | Macro F1 ↑ | AUC ↑ | Logloss ↓ | Brier ↓ | Evaluated gene-pair tests/patient ↓ | Internal nodes (mean) |
|---|---:|---:|---:|---:|---:|---:|---:|
| GLOBAL_LABEL_ONLY (identical global optimizer, w=0) | .834660 | .855530 | .880495 | .358958 | .106152 | 1.360140 | 1.7273 |
| **GLOBAL_RLR_LOGIT_OOF (w=.25)** | **.826367** | .850205 | .878139 | **.364068** | .108338 | **1.308192** | 1.6364 |
| GREEDY_LABEL_TSP_TOP20 (proxy, NOT real TSPDT) | .840313 | .864467 | .881017 | .349497 | .103147 | **1.298701** | 1.6364 |
| Original RLR20 | .908225 | .925915 | .944545 | .260502 | .063266 | 17.545455 | N/A |
| Real switchBox kTSP5 | .863363 | .881064 | .937638 | .289605 | .092827 | 5.000000 | N/A |

Mean distinct tree relations and different genes: v5 RLR 1.6364 / 3.0000; identical exact-global label-only 1.7273 / 3.1818. Full expression universe still needs measuring for within-patient ranks. The short path is NOT a claim of 2-pair wet-lab assays.

## Independent dataset-level tests, 11 units (2^11 exact sign flips, 8k seeded dataset bootstrap)

| A−B | Metric | Mean delta | Bootstrap 95% | Two-sided exact p |
|---|---|---:|---:|---:|
| GLVT_RLR − GLOBAL_LABEL | F1 | **−.008293** | [−.024880, .000000] | **1.000000** |
| GLVT_RLR − GLOBAL_LABEL | AUC | −.002356 | [−.007068, .000000] | 1.000000 |
| GLVT_RLR − GLOBAL_LABEL | Logloss | **+.005111** | [.000000, +.014428] | .500000 |
| GLVT_RLR − GLOBAL_LABEL | Pair checks | −.051948 | [−.155844, .000000] | 1.000000 |
| GLVT_RLR − GREEDY_LABEL | F1 | **−.013946** | [−.033757, .000000] | .500000 |
| GLVT_RLR − GREEDY_LABEL | Logloss | +.014572 | [+ .001181, +.035457] | .074219 |
| GLVT_RLR − original RLR20 | F1 | −.081858 | [−.157509, −.007871] | .070312 |
| GLVT_RLR − real kTSP5 | F1 | −.036996 | [−.085184, +.013663] | .195312 |

All p-values exploratory and not interpreted as evidence of equivalence. Cohorts, not folds or patients, are the primary inference units.

**Where RLR changes the outcome:** in 9/11 datasets both exact-global trees have exactly the same F1, AUC, logloss and path burden. West: F1/path identical, logloss worsened +.004971. Yeoh: RLR-guided F1 dropped **.933333 → .842105** (−.091228), AUC −.025915, logloss +.051246; pair paths shortened by .571429. There is **no dataset where RLR guidance improves F1**. The lower path burden is confined to Yeoh and obtained with a substantial F1 loss.

## Final scientific decision

**NEGATIVE, no scientific promotion.** This is a valid very small interpretable tree and an algorithmically correct globally optimized surrogate, yet RLR's OOF logit variance fails to add stable value over the *identical* global label-only search, and the simpler greedy TSP-like tree is at least as accurate with slightly fewer checks. V5 is not a defensible habilitation-specific advance based on current evidence.

Do not tune teacher weight, objective normalizations, leaf priors, depth or candidate panel on these repeatedly exposed V/E data. Do not run v5 on first E or full470. Do not call the greedy proxy the historical CUDA TSPDT. The original unpublished RLR remains preserved and offers the more mature paper path with actual kTSP10/15/20 and train-derived matched K.

Future *tree* research only if there is a genuinely different RLR-specific inferential principle and faithful comparisons with author's TSPDT/GTSPDT/EMTTree+RX, rather than optimization alone. Original processed ARFF/probe parity and an actually unseen independent cohort remain key publication gates.

**Audit archive:** [full run #40919146672](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40919146672), artifact rlr-glvt-v5-all11-independent-science-40919146672 (patient-score replay, cohort CSV and exact tests). This report preserves the main evidence beyond artifact expiration.
