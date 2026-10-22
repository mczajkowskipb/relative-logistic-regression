# RLMT-v3B — one-relational-split ridge offsets with global-calibration control

Declared **after**, not before, review of the EXPOSED v3 pilot.
This is research development, NOT a second independent validation.

## Real previous outcomes (2026-10-22)

The operator completed RLMT-v3 pilot 2 actual patient cohorts
(\`chen-2002\`, \`gordon-2002\`) ×3 original outer fold
positions ×3 algorithms, 18 real fits. All model RDS synthetic
replays and independent held-out patient-level metrics passed:
\`RLMT_V3_PILOT_BATCH_EXIT=0\`,
\`RLMT_V3_PILOT_INDEPENDENT_AUDIT=PASS\`.

| Variant | Trees / 6 | Logloss win/loss/tie vs aligned root | Average ΔLL vs archived original RLR |
|---|---:|---:|---:|
| STUMP_TSP_BETA | 6 | 0/6/0 | +0.21759099 |
| STUMP_RIDGE_OFFSET | 6 | 5/1/0 | **−0.00247027** |
| TREE_D2_PRUNED_BETA | 6 | 0/6/0 | +0.20396658 |

Crucial interpretation: plain trees with Beta-smoothed probabilities
are substantially inferior to well-trained global probabilistic
RLR on the examined folds. The slight positive offset-stump
signal may not come from a **useful relation**; one global
intercept correction alone may produce a similar gain.
The improvement has NOT been validated independently and
is not proof the offset stump is better than RLR.

## New science question: does the route matter beyond calibration?

Run four methods on the EXACT same outer training/test samples:
1. \`GLOBAL_RIDGE_CAL\`: calibrated original aligned parent
   \`logit(p_aligned) + b\`, scalar \`b\` penalized by fixed
   ridge λ=1. No route, 0 tree nodes.
2. \`STUMP_TSP_RIDGE\`: **bit-identical** v3 STUMP_RIDGE_OFFSET,
   top training-label TSP class contrast relation;
   route false and route true have two λ=1 ridge log-odds
   corrections to the SAME aligned parent.
3. \`STUMP_GINI_RIDGE\`: top train-only Gini reduction
   relation with same min leaf support and two-coefficient
   ridge correction. No separate Gini-based full RLR leaves.
4. \`STUMP_RESID_RIDGE\`: top parent **training OOF residual
   contrast** candidate from same bounded 4096-pair inventory,
   with same two penalized corrections. This is an ALTERNATIVE
   hypothesis, not a validated improvement.

Fixed constants: seed=42, deterministic 4096 sampled unordered
expression-feature pairs (fewer when fewer possible);
route is \`x_i > x_j\`; each side >=8 training samples;
ridge λ=1, alpha=0, no glmnet intercept, no standardization,
no per-node K=10 reranking; one independently trained parent,
one strict 3-fold honest parent OOF (train-only supervised
pair selection) **shared** by all methods per outer fold.

Per parent, all 4 methods consume the same OOF vector and
identical full outer-training records. The reference RLR
model is read from the original frozen 285-model archive
without recomputation or alteration. Original 19 cohort
source files verified by SHA256; original fold training
and heldout sample IDs verified byte-for-byte.
On a dataset where the outer-training minority class has
<9 observations, strict parent OOF is mathematically
infeasible under the existing source protocol. All 4
methods MUST transparently return the unchanged aligned
parent (no route), rather than censoring the dataset
or leaking in-sample probabilities to fit a calibrator.
This declared identical fallback still counts in the
denominator.

The one-parameter intercept calibration uses the convex
ridge logistic objective optimized by Newton to tolerance
<1e-10. The two-route coefficients use the frozen glmnet
ridge fit, never use outer test labels at prediction time.
An independent patient-level Python audit scores all saved
predictions, compares **against the global calibration
control** (primary distinction), aligned parent and
archived RLR, and reports equal dataset averages, 19
dataset win/loss/tie and actual number of tree rules.
The four method outputs share the same parent fit and
its warnings; warnings should NOT be added four times.

Important: Calibrating \`RLR_ALIGNED\` is not recalibrating
the archived published RLR model. Include both in tables.

## Fixed computational gate

- First real synthetic R/glmnet fit/replay and preserved
  exact v3 TSP-RIDGE outcome check.
- Then 19 original datasets ×5 original outer folds =
  95 common training bundles/380 method-fold prediction
  grids in a new outcome root.
- Then independent sample-level audit must PASS. Failure
  is preserved, not overwritten, never counted as success.
- No change in method constants or selecting favorable
  cohorts following observed outcomes.
- If route-specific models do NOT consistently improve
  over \`GLOBAL_RIDGE_CAL\`, report that and stop claiming
  DT benefits. Do not launch the full 470-fold paper
  benchmark just to seek a favorable result.

These 19 cohorts and fold outcomes have been repeatedly
inspected and influenced method development: positive
results remain EXPOSED RETROSPECTIVE DEVELOPMENT.
A genuinely new external cohort with prospectively
frozen variants is necessary for a trustworthy
out-of-development assessment.


## FINAL RLMT-v3B 19×5 retrospective control result (2026-10-22)

**Real operator-provided final logs, independently audited:**
\`RLMT_V3B_19X5_BATCH_EXIT=0\`;
\`RLMT_V3B_19X5_REAL_TRAINING_BUNDLES=95\`;
\`RLMT_V3B_19X5_METHOD_FOLD_RESULTS=380\`;
\`RLMT_V3B_19X5_INDEPENDENT_PATIENT_AUDIT=PASS\`;
\`RLMT_V3B_19X5_AUDITED_FOLD_METHODS=380\`;
\`SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY\`.

**Equal-dataset mean paired held-out logloss delta, 19 datasets:**

| Method | Real split folds / 95 | Dataset W/L/T vs GLOBAL_RIDGE_CAL | ΔLL vs archived ORIGINAL RLR | ΔLL vs GLOBAL_RIDGE_CAL |
|---|---:|---:|---:|---:|
| GLOBAL_RIDGE_CAL | 0 | 0/0/19 | +0.000606427 | +0.000000000 |
| STUMP_TSP_RIDGE | 90 | 11/7/1 | +0.000963708 | +0.000357282 |
| STUMP_GINI_RIDGE | 90 | 9/9/1 | +0.001073118 | +0.000466691 |
| STUMP_RESID_RIDGE | 90 | 10/8/1 | +0.001535369 | +0.000928942 |

Negative delta would be better. All four models are **worse on
overall dataset-equal mean logloss than ORIGINAL RLR**. All three
relational stumps are **worse on mean logloss than the no-tree
GLOBAL_RIDGE_CAL control**. Route creation frequency 90/95
does not establish predictive merit. TSP 11/7/1 dataset W/L/T
relative to calibration does not overcome a positive average delta.

Previous v3 2×3 EXPOSED pilot looked promising for TSP ridge:
mean ΔLL vs original −0.00247027, but the broader full
19×5 audit contradicted the favorable small-sample conclusion.
The plain Beta-frequency stumps and depth2-pruned Beta
trees were worse than original by ~+0.218 and ~+0.204
in that earlier 6-fold pilot and should not be promoted.

### SCIENTIFIC NO-GO AND NEXT APPROPRIATE RESEARCH QUESTION

**No-go** for superiority claims about v3B TSP/Gini/residual
stumps. Do NOT edit original paper to include a claimed DT
improvement or launch the 470-position CompCancer experiment
to mine a favorable result. No retroactive selection of 11
winning datasets or adjustment of λ/split score based on these
95 exposed outer fold outcomes.

**Identified comparability caveat to analyze, not conceal:**
The calibrated no-split model penalizes ONE global offset
coefficient by ridge λ=1; the two-branch stump penalizes
TWO separate route coefficients with the same λ. A common
offset duplicated on both branches incurs a different
penalty than a single global intercept. Consequently the
current comparison evaluates two entire procedures but
does not cleanly isolate the incremental conditional effect
of the relational route while keeping calibration fixed.

**Potential NEW predeclared study (not yet implemented):**
a hierarchical conditional calibration model
\`logit(p(x)) = logit(p_RLR(x)) + b_global + γ(z(x)-π)\`,
where \`z(x) = 1{x_i>x_j}\`, \`π\` is the TRAIN-only route
frequency, \`b_global\` is the SAME global calibration
coefficient in control/tree, and \`γ\` alone receives
ridge shrinkage toward zero. γ=0 is an EXACT nested
global-calibrated model. The pair selection and γ penalty
must be determined wholly within training using honest
cross-fitted predictions with per-fold supervised pair
re-selection. Evaluate also training-only soft-pruning
(shrink γ toward zero), Brier/calibration and no-split
rate. If no meaningful incremental value, stop; do not
force nonzero trees just to have DT in a paper.

For a publishable **new DT method**, use a genuinely new
independent validation source/holdout frozen BEFORE reviewing
new outcomes, along with reproduced historical TSPDT.
The original RLR manuscript is already complete and
should stay unchanged. This plan is hypothesis generation,
not evidence of new superiority.
