# RLMT-v4 — original-RLR-anchored rank-repair trees (2026-10-22)

**STATUS: PROTOTYPE CODE PREPARED; NO R EXPERIMENT RUN BY ASSISTANT.**
Do not claim a gain until synthetic parity, original patient-paired
parity, and independent heldout evaluation all complete on Ubuntu.

## Motivation: previous negative experiments

- v2 19×5×4 audited, TSPDT_GINI 7 accepted/95 with equal-dataset
  mean Δ logloss vs original RLR −0.001087952 but 7/10/2
  dataset win/loss/tie and LODO flips 4/19.
- v3B 19×5×4 audited:
  * GLOBAL_RIDGE_CAL ΔLL original +0.000606427;
  * STUMP_TSP_RIDGE +0.000963708 (90 real splits/95),
    Δ vs global calibrator +0.000357282;
  * STUMP_GINI_RIDGE +0.001073118 (90/95),
    Δ vs calibrated +0.000466691;
  * STUMP_RESID_RIDGE +0.001535369 (90/95),
    Δ vs calibrated +0.000928942.
  No superiority. This is NOT merely an inability to split; forcing
  actual stump predictions often harms a strong global model.

## New question — correct locally only where RLR is unreliable

Original RLR has K=10 sparse/ternary selected signed rank pairs.
Their per-patient activity (how many encode ±1 rather than 0)
and the trained model's probability provide two distinct sources
of reliability information. Rank-gap gives a third previously
discarded signal. Tree route is a GENE ORDER, never a raw value
threshold. Interpretability: 1 global original RLR, at most
1 named gene-pair route and 2 coefficient calibration parameters.

**CRUCIAL CHANGE:** Start from the *literal frozen original-paper
RLR fitting procedure* (switchBox top1000 + disjoint K10,
m10,d2, cv.glmnet alpha0.5 lambda.min nfolds=safe_nfolds,
set.seed(42)). Not the re-fit stratifed RLR_ALIGNED_V2 root
that failed as the comparator in v2/v3B.

We do not simply say the parent is original: the run FAILS CLOSED
if its heldout probabilities deviate by >1e-9 from archived
published 2011/CompCancer 19x5 RLR fold probabilities.

## Four corrections + exactly nested no-correction root

- ORIGINAL_RLR: immutable baseline function output.
- GLOBAL_CAL: logit(p_orig) + β.
- CONF_ROUTE: logit(p_orig) + β + γ
  [4p_orig(1-p_orig)] [I(g_a>g_b) − π_train].
- ACTIVITY_ROUTE: logit(p_orig) + β + γ
  [1 − fraction active RLR signed pairs] [I(g_a>g_b)−π_train].
- SOFT_GAP_ROUTE: logit(p_orig) + β + γ
  [tanh((rank(g_a)-rank(g_b))/10) − μ_train].

For each pair and gate π_train/μ_train are selected only on
outer-training data. Original parent and all corrections use
strict **training-only** 3-fold OOF with supervised pair ranking
RECOMPUTED inside each OOF training fold. Correction β,γ are
fit by the convex penalized logistic objective with fixed
λβ=.02 and λγ=.10; gamma→0 recovers exactly the no-route
calibration within that model.

At most 4096 deterministic label-agnostic unordered feature
pairs plus ten original RLR shortlist pairs are screened via
a standardized weighted covariance of OOF residual and
route-feature. Each candidate must have at least 8 samples
per route side but no 6-per-class constraint (no local
glmnet leaf training). After selecting a gene pair, fit a
low-dimensional coefficient, NOT two full local RLR models.
This increases interpretability and avoids class<8 warnings
in nested leaf CV.

The candidate selection itself is supervised, and 4096 choices
can overfit even with OOF predictions. OOF predictions avoid
training the parent on the same labels, but are NOT a substitute
for nested/independent feature-selection validation.
These models must be treated as exploratory and can be
improved only using prospective train-only nested selection
or new independent validation—not by viewing outer fold
outcomes and tuning λ/weight/candidates post hoc.

## Guarded pilot

1. Real synthetic 160/54 patient-like data: 5 model modes,
   genuine supervised training and 3-fold OOF, bitwise RDS
   replay, exact parity against fit_rlr_train_test wrapper,
   affine rank-order invariance.
2. chen-2002 and gordon-2002, first 3 original seed42 folds:
   six bundles ×5 methods = **30 real heldout prediction vectors**.
   Check immutable old input SHA256, frozen patient IDs,
   class labels and exact original root probability parity.
3. Python independent audit recomputes logloss/Brier
   from the six heldout patient sets, compares all modes
   to archived original RLR; includes negative cases;
   stop on mismatched p_original or warning/error.

Result directory NEW only:
\`results/rlmt_v4_rank_repair_2x3_v1\`.
The source v4 branch is separate from main, old 285 RLR
models, v2/v3/v3B completed results and paper.
No clinical data, external validation or demonstrated
method gain claimed before real run.

## The scientifically meaningful next decision

If one candidate is promising in this exposed six-fold screen:
pre-register **ONE** method identity and ONE endpoint,
run 19x5 as retrospective development (not independent),
and then test on a new dataset/batch/independent cohort
frozen *before* any label examination or hyperparameter
adaptation. Include full 470-position CompCancer if needed
for publication comparability, explicitly labeled
development-exposed. Check Brier, logloss, balanced accuracy,
macro-F1, calibration, subgroup support and gene-pair
stability. If no method exceeds original RLR in reliable,
paired external evaluation, report that honestly and
publish the already-complete original RLR manuscript.
