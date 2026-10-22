# RLMT-v6 — original RLR pair-conditioned rank-margin reconstruction

**Status 2026-10-23: CODE PREPARED, NOT EXECUTED ON UBUNTU.
NO CLINICAL OR BIOMEDICAL GAIN CLAIMED.**

## Evidence motivating direction (already EXPOSED)

The full v4 19×5 ×5 output has 475 independently audited
patient-level fold-method predictions and exact original RLR
parent parity. Equal-dataset mean ΔLL vs original published RLR:

- GLOBAL_CAL +0.001301927
- CONF_ROUTE +0.001615592
- ACTIVITY_ROUTE +0.001233117
- SOFT_GAP_ROUTE +0.050673367

The 2×3 v5 honest train-only threefold router produced
0/6 activated corrections, 0.0 difference from original RLR,
and full independent audit PASS. A read-only examination of
the already saved router files found:
- chen fold00 best raw gain −0.00052108;
- chen fold01 best raw gain −0.00088067;
- chen fold02 best raw gain +0.00036690;
- gordon fold00 best raw gain +0.00074018;
- gordon fold01 best raw gain +0.00081053;
- gordon fold02 best raw gain +0.00052579.
Thus two of six had no positive internal evidence, four had
small positive internal gains rejected by the conservative
SE-complexity gate. Do not retrospectively lower the v5
acceptance threshold to mine these exposed outer tests.

## Different v6 hypothesis

The original RLR uses up to 10 supervised gene rank-order
pairs, each quantized to {-1, 0, +1} using the published
m=10 rank separation threshold. This discards how far
apart the two chosen genes are within the same ternary
category. The v6 direction is to restore ONLY this
disregarded continuous within-category rank separation,
rather than adding thousands of randomly screened gene
pairs or fitting full leafwise RLR models.

For each original-RLR selected pair (a,b) compute:

    d = [rank(b) - rank(a)] / max(1, n_features - 1)
    t = ternary_RLR(a,b; published_m=10)
    u = d - mean_training(d | t)
    m = clamp(u / max(0.02,RMS_training(u)), -3, 3)

Train-only conditional means per pair and ternary state,
RMS and original gene order are saved with the model and
applied unchanged to outer heldout patients. All original
K10 pairs are from the exact original-paper root training
of each outer fold; no candidate features from outside the
original model are searched. This is within-patient
monotone-transformation invariant, with biological gene
names and rank distances directly inspectable.

Evaluate FIVE paired fixed candidate models:
1. ORIGINAL_RLR — literally the original manuscript
   switchBox / cv.glmnet RLR.
2. GLOBAL_CAL — original RLR + penalized intercept.
3. ONEPAIR_MARGIN — original RLR + intercept and only ONE
   selected rank-margin among the original K10 pairs,
   selected by training OOF residual contrast.
4. CONSENSUS_MARGIN — original RLR + intercept and ONE
   rank-margin summary weighted by the ORIGINAL_RLR
   learned signed logistic coefficients.
5. ALLPAIR_MARGIN — original RLR + intercept and <=10
   original-rank-margin predictors with ridge shrinkage.

Penalized binomial offset loss on original-RLR 3-fold
train-only OOF probabilities, lambda intercept .02,
lambda additional margin gamma .15, log-odds correction
capped to [-2,2] at prediction. Single and multi
coefficient corrections all share the SAME original root;
if gamma is zero and intercept is held fixed they reduce
to the same global calibration (subject to identical cap).

Scientific distinction: v6 is NOT a hard-split binary
decision tree. It is a *transparent RLR with a low-dimensional
continuous rank-margin correction*; a one-pair variant can be
represented as a soft rule/decision gate.

## Mandatory honesty / limitations

- All 19 prior CompCancer cohorts and especially the six
  development pilot folds are previously exposed. Even a
  positive v6 score is not an independent scientific
  validation. Novel cohort/laboratory/batch replication is
  required before any improvement claim.
- The original candidate K10 genes are selected using ALL
  outer-training labels. While the crossfit original-parent
  probabilities come from genuinely independent 3-fold
  internal fits, features derived from the complete outer
  train's supervised original K10 shortlist are not themselves
  reselected independently inside each OOF holdout fold.
  Do not mislabel this as a fully nested original-pair
  feature-selection evaluation. This can overfit in development,
  and must be tested on true outer holdouts and independent data.
- One selected margin by correlation with y−pOOF can also
  overfit (though searches <=10, not 4096); compare
  CONSENSUS_MARGIN which does not inspect OOF labels again
  for pair choice. ALLPAIR_MARGIN can overfit 10 coefficients.
- Do not prechoose a winner by sorting the six outer-fold
  scores and then call a larger run confirmatory.
- A stable improvement beyond BOTH original RLR and
  GLOBAL_CAL, across multiple independent cohorts with
  unchanged parameters, is a requirement for a novel
  predictive method claim.

## Frozen first gate

Real 175/52 patient-like *synthetic* data with real switchBox,
cv.glmnet and independent original-wrapper parity, fixed
3-fold original RLR OOF, five predictive modes, strict RDS
replay and rank-preserving affine invariance. Then exposed
chen-2002 and gordon-2002 3 original folds each,
six shared training bundles × five methods=30
method-fold pred vectors, with source SHA256, original
frozen train/test sample IDs and exact original heldout
RLR probability parity. Finally an independent Python
patient-level saved prediction audit of logloss/Brier and
original alignment. NO modifications to original archive.
Do not start full 19x5 until this synthetic and pilot
gate has passed and inspection warrants it.
