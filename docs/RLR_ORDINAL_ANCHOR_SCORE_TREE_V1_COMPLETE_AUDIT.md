# Ordinal Anchor-Score Decision Tree V1 — COMPLETE real-data evidence, negative science gate

**Execution:** 2026-10-27, [GitHub Actions #40467789138](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40467789138), **ALL JOBS SUCCESS**, fitted commit `c13fc3da469ddd41519dfac1ae32f8fb9f67c20a`. This document only reports that frozen run; no changes to the original RLR source or locked unpublished paper. Predeclared [protocol](RLR_ORDINAL_ANCHOR_SCORE_TREE_V1_PROTOCOL.md).

## Original user hypothesis

Compare one potentially changing target gene against several relatively stable reference genes in each **individual patient**. Instead of restricting the result to majority voting, define ordinal anchor score `S(g,A,x) = sum_{a in A} 1[x_g > x_a]`, with A containing three or five distinct reference genes. Let a node learn a threshold `1 <= t <= |A|`, splitting on `S>=t`. This is an integer score purely from **within-sample ordering**, with no raw expression values or ratios. A fixed threshold corresponds to a chosen order statistic of the references, so it can still be effectively redundant with a single pair in some populations.

## Algorithm and controls

- RLR baseline is unchanged: K20, supervised switchBox top1000, strict d2, signed ternary m10, elastic-net alpha0.5 and inner lambda.min.
- Target genes are top supervised pair hubs (up to 14) and reference stars are constructed from original pool neighbors and **train-only rank** stability/diversity. For each fixed reference panel of 3/5, threshold t searched on D using Laplace/empirical-prior-smoothed logloss and complexity penalty 0.0025*(k-1). Minimum 10 patients and 2 per class per root branch.
- Actual top-down child induction is attempted separately in each D root branch: n>=22, class counts>=5, leaf size>=6, within-child D gain>=.004, new local supervised pairs and a new reference panel/threshold. The system does NOT force local children when unsuitable.
- Inner D/V selects depth2 only if it improves on its own stump by >=.005 V logloss; tree acceptance requires >=.005 V logloss improvement vs RLR, >=.002 vs same-root best single-anchor pair and BA loss no greater than .02. Otherwise **exact RLR fallback**.
- Seven arms: ORIGINAL_RLR, ANCHOR_SCORE_STUMP, ANCHOR_SCORE_TREE_D2, SAME_ROOT_BEST_SINGLE_PAIR, SAME_COUNT_RIDGE, V_SELECTED_ANCHOR_SCORE_TREE, RLR_V_GATED_ANCHOR_SCORE_TREE. Global ridge is a control using the same node count features without routing (fixed lambda=.10, inert all-zero column only for glmnet's >=2 predictors restriction); it is not perfectly matched in degrees of freedom to the leaf models.

## Independent experiment verification

Public Schliep source-reconstructed Chowdary-2006 (104 patients/182 genes) and Singh-2002 (102 patients/339 genes) are two **reused original-paper development cohort identities**, not fresh validation or confirmed exact processed ARFF equivalents. Fixed 3-fold E splits (seed42) on each cohort. Six splits, 42 model-fold metrics, **1442 archived heldout patient-method probability rows**. Separate Python auditor reproduced for *every patient and arm* heldout raw Brier, clipped logloss, balanced accuracy, Macro-F1, delta loss, exact positive-class mapping, patient coverage, V-depth selection and bit-exact V-rejected original RLR fallback: `OASDT_1442_HELDOUT_PATIENT_SCORES_REPLAYED=PASS`, `OASDT_42_METHOD_FOLD_MACROF1_BA_LOGLOSS_BRIER_AUDITED=PASS`, `OASDT_ALL_PATIENTS_COVERED_AND_FALLBACK_EXACT=PASS`.

Mechanical tests PASS for 130 independently monotone-transformed patient vectors, rank-count and branch parity, integer threshold, score majority-equivalent to anchor median at t=3,k=5, one-anchor perturbation bound, fixed gene schema and a planted manually specified *true* depth2 tree with distinct routes (this is only a software functional test, not proof of discovery by top-down induction).

The first CI run was stopped by `glmnet` requiring two predictor columns for the one-score no-tree comparator. Adding a permanently-zero second column changed **no predictive mathematical model**. The next run fitted all six folds but inherited the wrong historical auditor script in the workflow; this was a **CI invocation error**, not a model-fitting failure. Final run #40467789138 has both full R fit and independent proper OASDT auditor SUCCESS. Frozen final fitted SHA above.

## Fold-level outcomes (delta E logloss: new minus original RLR, negative better)

| Dataset | E fold | Root anchors k | Count threshold t | Actual local child nodes | V decision | Inner V gain vs RLR | Ungated E stump deltaLL |
|---|---:|---:|---:|---:|---|---:|---:|
| Chowdary | 0 | 5 | 1 | 0 | reject | -.271 | +.158254 |
| Chowdary | 1 | 3 | 2 | 0 | reject | -.272 | +.289936 |
| Chowdary | 2 | 3 | 1 | 0 | reject | -.562 | +.264521 |
| Singh | 0 | 5 | 3 | 0 | reject | -.0169 | +.274368 |
| Singh | 1 | 3 | 1 | 0 | reject | -.523 | -.560942 |
| Singh | 2 | 5 | 2 | 0 | reject | -.483 | +.320216 |

The ordinal score threshold is genuinely **not always majority**: 4/6 selected root thresholds differ from majority of the chosen k. Nevertheless additional capacity did not generalize. The root was learned in 6/6 but **no local child passed the frozen feasibility/gain criteria in ANY fold**. ANCHOR_SCORE_TREE_D2 predictions equal ANCHOR_SCORE_STUMP, not evidence of successful deeper tree induction. No V gate accepted a tree, so final V-gated model falls back to RLR for all E patients.

## Six-fold unweighted mean heldout E relative to original RLR (negative better)

| Method | Delta logloss | E folds improved | Delta balanced accuracy |
|---|---:|---:|---:|
| Original RLR | 0 | — | 0 |
| ANCHOR_SCORE_STUMP | +.1243922 | 1/6 | -.0936411 |
| ANCHOR_SCORE_TREE_D2 | +.1243922 | 1/6 | -.0936411 |
| SAME_ROOT_BEST_SINGLE_PAIR | +.0853547 | 1/6 | -.0682938 |
| SAME_COUNT_RIDGE | +.0596607 | 1/6 | -.0787018 |
| V_SELECTED_ANCHOR_SCORE_TREE | +.1243922 | 1/6 | -.0936411 |
| RLR_V_GATED_ANCHOR_SCORE_TREE | 0 | 0/6 | 0 |

## Scientific decision and next gate

**NEGATIVE for this formulation on the two exposed development cohorts.** We have proven only that the ordinal-count induction is *implemented correctly* and samplewise rank invariant, not that it improves original RLR or yields learned depth2 trees. The only large E fold improvement (Singh fold1) was **correctly rejected by the inner V rule**, consistent with unstable small-sample generalization. Do not retrospectively tune k, root thresholds, depth criteria, smoothing or V acceptance using the observed E results.

A reasonable distinct next question is whether there is any reliable branch-conditional information at all after original RLR, rather than repeatedly fitting tiny sample-leaf models. Within-training D-only cross-fitted RLR residual analysis, reference stability stratified by class, local interaction permutation checks and a no-tree interaction control could test this without another E-search. If this mechanism has no stability, additional EA is unlikely to help and risks exploiting noise.

**Original completed RLR and locked manuscript are untouched and remain the scientifically strong fallback**; original CUDA TSPDT not run here and cannot be replaced with greedy proxy in a final paper. No clinical interpretation, novelty proof, external replication, or improvement claim arises from these results.
