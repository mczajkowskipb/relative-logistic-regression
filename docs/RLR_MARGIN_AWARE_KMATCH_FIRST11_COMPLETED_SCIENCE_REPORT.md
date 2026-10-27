# First outer E fold complete — margin-aware RLR vs actual capacity-matched kTSP (2026-10-27)

**Provenance:** all 11 original manuscript-named CompCancer public Schliep reconstructed expression matrices, not original processed paper ARFF SHA verified; all these datasets previously examined in method development, so this is an **exploratory SAME-SPLIT benchmark**, NOT an independent confirmation. Full original locked manuscript **UNPUBLISHED and untouched**.

**Real CI:** [GitHub run #40588678830](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40588678830), **12/12 successful** (11 cohorts + independent Python patient-level audit). Exactly 11 named cohort first original outer stratified seed42 folds, 12 arms per cohort = **132 independent replayed model-cohort metrics** from archived E patient probabilities and native kTSP decisions. Training-only K budgets checked against archived T-only RLR global coefficients and train per-patient activity. Original exact RLR training code frozen unchanged. [Algorithm and matchup protocol](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-margin-aware-budget-match-v1/docs/RLR_MARGIN_AWARE_CAPACITY_MATCHED_KTSP_PROTOCOL_V1.md).

## Dataset-equal-weighted FIRSTFOLD metrics

| Arm | Positive F1 ↑ | Macro-F1 ↑ | AUC ↑ | LL ↓ | Mean global nonzero | Mean E patient active |
|---|---:|---:|---:|---:|---:|---:|
| Original RLR20 | .868613 | .904867 | .954134 | .274682 | 16.818 | 15.366 |
| **New margin-aware RLR20** | **.877746** | **.908629** | **.976181** | **.264319** | 16.364 | 15.136 |
| Original RLR50 | .867924 | .897200 | .967657 | **.225754** | 29.091 | 26.516 |
| New margin-aware RLR50 | **.874418** | **.906810** | **.973909** | .230255 | 27.818 | 25.610 |
| Actual kTSP10 | **.901340** | **.918630** | .963302 | .301085 | 10 | 10 |
| Actual kTSP15 | .899870 | .918182 | .963331 | .328138 | 15 | 15 |
| Actual kTSP20 | .855341 | .887385 | .958391 | .340390 | 20 | 20 |
| Actual kTSP25 | .848381 | .876587 | .959063 | .354925 | 25 | 25 |
| Train matched kTSP to RLR20 global | .865401 | .893007 | .961035 | .326712 | 16.818 | 16.818 |
| Train matched kTSP to RLR20 active | .875275 | .898043 | .962729 | .321320 | 15.545 | 15.545 |
| Train matched kTSP to RLR50 global | .894134 | .908339 | .963575 | .349390 | 29.091 | 29.091 |
| Train matched kTSP to RLR50 active | .896743 | .913129 | .963796 | .350066 | 26.818 | 26.818 |

**Do not pool 11 outer folds as 11 independent replications of a tuned model:** each differs in patients/dataset/class support and the same cohorts were exposed earlier.

## Dataset-level exact uncertainty on FIRSTFOLD only

Computed directly from true 11 per-dataset E metrics in complete CI artifact, paired exact 2^11 two-sided random-sign test, and 15k seeded dataset-level bootstrap resamples. Exploratory and multiplicity-unadjusted:

| Paired outcome | Mean Δ | 95% dataset bootstrap | Exact signflip p |
|---|---:|---:|---:|
| MarginRLR20 − original RLR20 positive F1 | **+.009133** | [−.01636,+.04870] | **1.000** |
| MarginRLR20 − original RLR20 AUC | **+.022047** | [+.00117,+.04659] | **.1406** |
| MarginRLR20 − original RLR20 LL | **−.010363** | [−.10458,+.07964] | **.8369** |
| MarginRLR20 − actual kTSP15 F1 | −.022124 | [−.10604,+.05313] | .6875 |
| MarginRLR20 − actual kTSP20 F1 | +.022405 | [−.02620,+.07477] | .4375 |
| MarginRLR20 − same-T-budget kTSP20-global F1 | +.012345 | [−.04377,+.06962] | .7031 |

**F1 gains do NOT pass any 11-cohort-level significance check.** AUC mean improvement is promising but exact sign-flip p=.1406 is not persuasive for declaring a winner. With only 11 datasets and many ties, bootstrap percentiles and exact signflip tests may differ; report both, do not cherry-pick the one excluding zero. The dramatic 11-fold average partly reflects hard datasets, but kTSP10/15 remain stronger than new RLR20 in firstfold F1. Pomeroy 7 E patients is too small for an isolated categorical conclusion.

### What is actually new versus original complete RLR

New extension is not "RLR K=50 with few real pairs"; original **already** explicitly has that result. New test asks whether selecting the SAME original top1000 switchBox candidate pairs by class-conditional separation of the SIGNED-MARGIN **actual original ternary feature** `{-1,0,+1}` with finite inactive-support penalty improves predictive performance. Original disjointness, active coefficients, logistic fit and probabilities are unchanged. This is a modest *RLR feature-selection refinement*, and must not be called a newly learned tree. No absolute expression normalization/range, ratios, or between-sample thresholds.

Historical references: Tan et al. original kTSP k≤10; switchBox default k2–10; Schliep et al. k10/25/50. **Position kTSP50 as supplementary**, PRIMARY comparison predeclared kTSP10/15/20 and train-matched effective K, plus original K50 results for completeness. Do not delete old benchmark or alter original locked manuscript in place.

## Repeated 470 fold validation actually launched

[Actual full 470 workflow #40589807042](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40589807042): **94 (cohort, repetition) jobs × five original 5fold E splits = 470**, 12 arms including the train-only matched-K controls, and independent Python auditor to replay **5640 method-fold patient score summaries** and check every subject held out once per repetition; sample-size-weight independent **11-dataset averaging**, plus paired 11-dataset exact signflip. The first attempt #40589517311 failed before model fitting because a CI Rscript filename was misspelled; fixed *only the path*, not scientific code, and fresh full run #40589807042 started. **Do not report full470 outcomes as complete until all CI+auditor SUCCESS.** Still public reconstructed inputs and repeated E already exposed. No significance/clinical validation claim.
