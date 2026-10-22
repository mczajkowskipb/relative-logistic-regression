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
