# QRLR V1 — class-blind per-pair rank-quantile abstention on the 11 original cohorts: completed D/V science gate

**Date:** 2026-10-27. Original RLR unpublished and immutable. [Full actual GitHub Actions #40594140376](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40594140376): **12/12 jobs SUCCESS**, eleven original named CompCancer source reconstructions plus fully independent patientwise Python recalculation of **66 method×cohort inner validation metrics**, exact patient D/V split auditing, q0 exact original RLR predictions, monotone nonlinear per-patient rank invariance. **Original outer E UNREAD** by QRLR. Original processed paper ARFF SHA parity remains NOT VERIFIED; cohorts had been used in other development tests, so these are exploratory.

## Proposed algorithm

`r_{ij}(x;q)=sign(rank_x(j)-rank_x(i)) * I(|rank_x(j)-rank_x(i)| >= m_{ij}(q))`, `m_{ij}(q)=max(10,Q_q^{D}(|rank_x(j)-rank_x(i)|))`, q fixed in {0,.25,.5}, type-1 quantile from class-blind D only. Original switchBox top1000 and disjoint K=10/20, original alpha=.5 glmnet CV lambda.min unchanged. q=0 returns original m=10 with exact literal original RLR probability parity; no numeric raw expression/ratio cutoffs or rescaling, only within-patient gene ranks. Quantile thresholds are pair-specific and train-only (held-out patient labels unused). This is **NOT A LEARNED TREE** and must not be called a tree-based classifier.

## Audit-confirmed 11 equal dataset means on first outer TRAIN D/V, not the 470 E test folds

| K | Within-D pair quantile q | F1 positive | Macro-F1 | AUC | Logloss | Nonzero locally per V patient | Globally nonzero coefficients |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10 | 0 (original) | **.902543** | **.911109** | .918283 | .363898 | 8.839 | 9.636 |
| 10 | .25 | .841937 | .874745 | .937778 | **.323172** | 6.987 | 9.545 |
| 10 | .50 | .834079 | .861959 | **.942349** | .417188 | **4.436** | 9.091 |
| 20 | 0 (original) | **.908225** | **.925915** | .944545 | **.260502** | 16.048 | 17.545 |
| 20 | .25 | .853680 | .889394 | **.953581** | .261439 | 11.848 | 16.091 |
| 20 | .50 | .808596 | .854366 | .941011 | .355078 | **7.818** | 15.818 |

**Quantify genuine sparsity:** K20 q25 reduces V per-patient active signed ternary relations by **4.200** (16.048→11.848) and K20 q50 by **8.230** (16.048→7.818); exact 11-dataset sign-flip p=.00293 and .00098 on activity, respectively, expected construction effect. But both lower F1 by **−.054545** (q25, p=.5 with dataset ties and sparse small validation) and **−.099629** (q50, p=.02344 unadjusted, multiple candidate tests). K10 q25 similarly reduces active 8.839→6.987 with F1 .9025→.8419. AUC of K20 q25 increases .00904 but exact 11-dataset p=.625; K10 q25 increases .01949 but p=.25. This is an INTERESTING **uncertainty and abstention vs discrimination** tradeoff, not a proven better classifier.

**Critical distinction:** global nonzero coefficients are still ~16.1 for q25, so QRLR must still CHECK them to determine whether they are inactive in this patient. A lower active rank-term count is a smaller effective explanatory score, **not** a matching reduction in actual gene-pair comparisons or number of assay genes measured.

## Scientific decision

Quantile-adaptive ternary rank bands remain strictly patientwise monotone invariant and mathematically faithful to original RLR, but aggressive pair-specific abstention **did not preserve F1**. Do not tune q by outer E outcomes, do not claim AUC gain from this 11 small inner V snapshot, and do not push q25/q50 as a publication winner. Maintain q0 as original RLR baseline. Original RLR20 with K15/20 and train-matched K kTSP as PRIMARY article framing remains stronger and ready for work independent of exploratory tree and quantile extensions.

Another research route if pursuing improved classification is to compare **properly train-only OOF class-decision threshold calibration and original RLR** under a method-symmetrical protocol; q25 AUC-F1 disagreement might be a calibration/threshold phenomenon, but this is HYPOTHESIS ONLY. A threshold adjustment alone is prior art and not a new learned relational tree. For habilitation-specific novel tree extensions require genuine conditional predictor value beyond global no-tree comparator; previous SAVS and certified contextual C-RET were not convincing.

Evidence: [scientific full470 RLR–kTSP conclusion](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-margin-budget-full470-v1/docs/RLR_MARGIN_AWARE_BUDGET_FULL470_FINAL_SCIENCE_GATE.md), [QRLR fixed protocol](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-quantile-ordinal-margin-v1/docs/RLR_QRLR_PAIR_SPECIFIC_ORDINAL_ABSTENTION_V1_PROTOCOL.md), and exact immutable GitHub artifacts [run #40594140376](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40594140376).
