# Anchor-majority gene top-down tree V1 — complete real evidence

**Date:** 2026-10-27 (Europe/Warsaw; GitHub workflow run timestamp 2026-10-26 UTC). **Technical gate: PASS. Scientific gate: NO E superiority.** Parent original RLR source and locked unpublished article unchanged.

## Scientific hypothesis

An active gene g that appears in many highly ranked gene pairs can be compared with a train-selected group of relatively stable reference (anchor) genes. A node votes 1 when g exceeds a strict majority of anchor genes (3/5 or 5/5). Prediction depends only on ordering within the *individual* patient, never on original expression magnitude or ratios.

Formally, v(g,A,x) = 1[sum_{a in A} 1(x_g > x_a) > |A|/2]. For odd |A| and no problematic ties, this equals 1[x_g > median_{a in A}(x_a)]. This is an ordinal median-comparison test, NOT 3 or 5 statistically independent votes.

**Important degeneracy:** If the anchor genes preserve a fixed relative ordering across all patients, their median is the same single gene in every sample. The group majority then reduces identically to a single TSP comparison with that median anchor. Nontrivial group behavior requires anchor ranks to re-order across samples; the balance between rank stability and reference diversity must be audited, not assumed. A common hub gene is not necessarily the causal driver; data alone cannot demonstrate causality.

## Frozen v1 algorithm

- Parent RLR is the existing K20, m10 signed ternary, top1000 switchBox supervised candidate pool, gene-disjoint d2, glmnet alpha .5 lambda.min.
- Candidate targets are the top-14 gene hubs by frequency of appearance among top-1000 supervised switchBox pairs, **before the original disjointness restriction**.
- Reference candidates are same-hub TSP neighbors and up to 70 genes with D-only stable within-patient normalized ranks and proximity to the target's class-rank midpoint.
- Greedy panel sizes 3 and 5; anchor score rewards TSP pair discrimination and penalizes normalized rank drift between classes, rank spread and distance from target rank. A near-duplicate comparison filter is applied.
- The root is selected for smoothed logloss on inner D, requiring at least 8 training samples in each branch.
- Separate node searches on the D subset of each root branch can add one level with its own anchor hub / majority rule. A child needs >=24 instances, both classes >=5, two child leaves >=7 instances and D smoothed logloss gain >.025. These conditions are part of the original frozen attempt, not posthoc.
- A single inner V chooses STUMP versus DEPTH2 when depth2 improves V logloss >.005. Outer E never influences structure or acceptance.
- Original RLR fallback selected if inner V logloss advantage over D-trained RLR <.005 or inner V BA loss >.02.
- Control uses ridge logistic regression on precisely the chosen anchor pair-order indicators, but **without tree hierarchy**. This is not an RLR-matched same-capacity control, so inference about unique tree benefit remains limited.

## Executed trial

**Real GitHub Actions:** https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40423066025 . Fitted SHA: `872f6d85be165bb8908bb9321d72f3ff0d68c2a4`. COMPLETE SUCCESS.

The two cohorts, Chowdary-2006 (104 patients, 182 features) and Singh-2002 (102 patients, 339 features), are **reconstructed public source matrices corresponding to names already in the original RLR manuscript**. Hash-equivalence to the original processed ARFF data is not established. Each was already inspected in previous RLR/RLMT development. Six frozen 3-fold outer E splits: 36 method-fold scores, 1236 E patient-method predictions, 100% independent Python replay of patient ID/class/probability, raw Brier, clipped logloss and BA, V-gate exact fallback PASS. Mechanical unit tests PASS for independent patientwise strictly increasing affine transforms, fixed gene schema, reference vote perturbation, median equivalence, and rank-invariant anchor selection. These are not fresh validation data.

## Per-fold diagnostic

| Cohort | Outer fold | V gate | Chosen depth | Child nodes | V gain vs RLR | Ungated anchor E delta logloss |
|---|---:|---|---|---:|---:|---:|
| Chowdary | 0 | reject | STUMP | 0 | -.002504 | -.082689 |
| Chowdary | 1 | accept | STUMP | 0 | +.009618 | +.119758 |
| Chowdary | 2 | reject | STUMP | 0 | -.031840 | +.181121 |
| Singh | 0 | reject | STUMP | 0 | +.162100 | +.078772 |
| Singh | 1 | reject | STUMP | 0 | -.033830 | -.301302 |
| Singh | 2 | accept | STUMP | 0 | +.173800 | +.085557 |

A positive V gain alone did NOT determine acceptance: BA can reject (Singh fold0); both accepted folds have worse E logloss than original RLR.

**Critical:** Zero children were added in all six outer folds. ANCHOR_TREE_D2 is numerically identical to ANCHOR_STUMP in this trial. This provides evidence about majority-anchor stumps and the feasibility of the fixed depth2 induction conditions, **not a substantive evaluation of deeper anchor trees**.

## Six-fold mean outer E logloss difference vs original RLR (negative better)

| Arm | Mean delta LL | E wins |
|---|---:|---:|
| ORIGINAL_RLR | 0 | — |
| ANCHOR_STUMP | +.0135362 | 2/6 |
| ANCHOR_TREE_D2 | +.0135362 | 2/6 |
| SAME_COMPARISONS_RIDGE | +.0518816 | 1/6 |
| V_SELECTED_ANCHOR_TREE | +.0135362 | 2/6 |
| RLR_GATED_ANCHOR_TREE | +.0342192 | 0/6 |

## Interpretation and next science gate

- **NO consistent improvement over original RLR.** Two spectacular-looking wins and several large losses emphasize small-sample selection instability; tuning against these exposed E folds would be invalid evidence.
- The original PI hypothesis ("top-pair hub often reflects one changing marker plus relatively stable reference genes") is **not yet directly verified** by hub counts or the tree experiment. V1 reference selection penalizes rank drift; it does not enforce strict stationarity, nor can rank-only data prove raw-expression stability.
- Strong immediate diagnostic question: is the identity of the median reference gene nearly constant across samples, causing the majority group to collapse to one ordinary gene-pair test? Evaluate this using TRAINING DATA ONLY without accessing new E labels, along with differences in target vs reference rank drift.
- It is legitimate to continue studying this **distinct reference-gene mechanism**, but new models must be predeclared, use E-unseen cohorts for eventual confirmatory claims and control ordinary TSP, grouped-gene TSP variants and a model with the same pair comparisons but no hierarchy.
- Prior work includes k-TSP and less restrictive shared-gene kTSP variants, stable-reference-gene single-sample signatures and top-scoring-pair-of-groups ideas; anchor majority alone is not obviously algorithmically original.
- The original complete RLR article remains untouched as a publication fallback. No results about actual original CUDA TSPDT were produced here.
