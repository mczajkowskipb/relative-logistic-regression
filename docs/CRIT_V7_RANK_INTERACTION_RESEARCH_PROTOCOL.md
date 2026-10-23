# CRIT-v7 — Conditional Rank Interaction Trees (research hypothesis)

**Status 2026-10-23:** independent SYNTHETIC probe executed locally;
code saved on dedicated branch. No biological cohort fitted; no
ORIGINAL_RLR comparator fitted; no "better than published RLR" claim.

## Why NOT another global correction

The operator completed full, patient-level independent audits for
19 previously EXPOSED biomedical cohorts ×5 original heldout folds:
- v4 GLOBAL_CAL mean ΔLL +0.001301927 vs original RLR;
  CONF_ROUTE +0.001615592; ACTIVITY_ROUTE +0.001233117;
  SOFT_GAP_ROUTE +0.050673367.
- v5 2-cohort 2×3, honest 3-fold selection, produced 0/6
  activated correction; original-RLR exact predictions.
- v6 19×5: GLOBAL_CAL +0.001384147;
  ONEPAIR_MARGIN +0.005470103 (8 cohort wins / 11 losses);
  CONSENSUS_MARGIN +0.001894404;
  ALLPAIR_MARGIN +0.010726693.
- The v6 2×3 pilot produced 6/6 ONEPAIR logloss wins, but
  those wins did NOT generalize to the other exposed cohorts.
This is explicitly exploratory, not a prospective comparison.

The problem with all these approaches is their largely ADDITIVE
correction to a strong global RLR. We need a falsifiable
nonadditive data-generating mechanism instead of another
posthoc tuning search over the same 19 datasets.

## Specific CRIT hypothesis: relation signs reverse with context

Let g_ab(x) = sign(x_a - x_b) and r_cd(x) =
sign(x_c - x_d), interpreted as binary relation values ±1.

A planted conditional interaction:

  logit P(Y=1|x) = 2 * g_ab(x) * r_cd(x)

has zero marginal signal for each individual pair when the
four gene measurements are independent continuous noises and
the two relational signs are balanced. Any ADDITIVE
\`intercept + beta_g*g_ab + beta_r*r_cd\` cannot express the
conditional slope reversal, while a genuine two-leaf
model tree with separate leaf slopes can.

Candidate interpretable real method:

  z(x) = I(expression[A]>expression[B])
  logit p_CRIT(x) =
      logit(p_ORIGINAL_RLR(x)) +
      beta_{z(x)} +
      sum_{j in S_{z(x)}} gamma_{z(x),j} * r_j(x)

with 2-5 rank pairs per leaf as a DESIRED complexity range,
but one per leaf for the first technical feasibility gate.
This models *contextual effect modification*, unlike
GLOBAL_CAL or a global rank-margin feature. The split must
be an actual interpretable rank comparison. If no robust
interaction: no split, EXACT original RLR.

**Critical new insight:** pairs that participate ONLY in
interactions can have NO useful marginal TSP score. Restricting
all candidate gates and leaf features to the original
supervised TSP/RLR K10 shortlist would systematically
miss pure sign-reversal XOR contexts. A defensible
candidate mechanism must include some labels-independent
balanced/stable relation pool or a pathway-prior pool, with
the full shortlist produced only within each outer train.

## Synthetic feasibility study performed here

Code: \`scripts/crit_v7_synthetic_regime_switch_probe.py\`

10 independent Gaussian gene measurements; 45 comparisons,
±1 gene-order encoding. **Baseline is a strong ADDITIVE
ridge logistic using ALL 45 comparisons: NOT published RLR.**
True model regimes:
- SWITCH: true logit=2*g_(0,1)*r_(2,3).
- GLOBAL: true logit=1.4*r_(2,3), no interaction.
- NOISE: true logit=0, no prediction signal.

In each synthetic replication, 300 TRAIN and separate
1400 TEST patients. Only TRAIN labels and samples are used
to choose 1 real split and 1 leaf predictor pair.
Discovery D=68% of TRAIN; validation V=32% of TRAIN:
- Fit additive pair-logistic baseline on D.
- Screen all DISJOINT candidate gate/model pair products
  via the score |sum_{D}(y_i-p_global_i) g_i r_i|.
  This is coefficient instability, not scalar residual
  class contrast.
- Only 8 top D-scored pair combinations are fitted as
  2-leaf penalized logistic-offset trees with ridge=.08.
  Each leaf has an intercept and ONE relation coefficient.
- Select the best V conservative logloss improvement,
  mean paired gain minus 1.5 SE(gains); accept only if
  >0.015. With no convincing benefit, predict exact
  additive global baseline (zero branches).
- Refit the chosen model on all TRAIN, inspect only
  then the independent TEST labels.
- Logistic baseline C=.23. Exact gene-pair direction,
  train/valid split seeds, penalties and thresholds
  are fixed in the stored script.

Early exploratory 50 repeats seed 1200..1249 WITHOUT
the 1.5 SE penalty: SWITCH 50/50 selected, 50/50
heldout wins, mean ΔLL=-0.2565833; GLOBAL 3/50
selected, two heldout losses and one win. Adding
an independent noise negative control found undesirable
false splits ~9/50 in that early code. Therefore the
SE penalty was introduced during **synthetic development**,
BEFORE the next seed set. Do not disguise this tuning.

Fresh disjoint 100-replication CHECK seed 8100..8199,
with 1.5 SE acceptance criterion now frozen, executed
locally (numpy/scipy/sklearn):
- SWITCH 100/100 true interaction pair sets discovered,
  100/100 splits accepted, 100/100 heldout wins,
  equal-rep mean ΔLL **−0.2579684** vs additive control;
- GLOBAL 0/100 accepted, 0/100 improvement or
  deterioration (all exact additive fallback);
- NOISE 0/100 accepted, 0/100 improvement or
  deterioration (all exact additive fallback).

These results show **representational feasibility
under deliberately easy, strong planted synthetic
nonadditivity**, NOT efficacy on biology, not original
RLR outperforming, and no novelty proof.

Additional exploratory weak-strength inspection (NOT
independent validation) with gate×relation coefficient:
0.35: 1/70 activated, mean ΔLL -0.000219;
0.60: 25/70 activated, mean ΔLL -0.014099;
1.00: 66/70 activated, mean ΔLL -0.08483;
1.40: 70/70 activated, mean ΔLL -0.161742.
Thus conservative abstention has limited power for
weak interactions at N=300: another reason not to
expect universal biomedical benefit.

## Next technical gates — not yet implemented

**Gate 1.** Replace synthetic 45-pair additive baseline
with the *LITERAL published original RLR* wrapper:
switchBox top1000, K10, m10, d2, original cv.glmnet
training and archived heldout probability patient-ID
parity. No RLR_ALIGNED_V2 surrogate.

**Gate 2.** Redesign candidate generation to scale to
10k+ genes without exhaustive O(G^4) products or
data-snooping:
  - compare label-free balance/rank-repeatability pool;
  - separately predeclare a biological-pathway-constrained
    pool; avoid choosing the winner using exposed cohorts;
  - choose a modest discovery candidate budget;
  - preserve candidate selection ONLY inside outer train.
The shortlist MUST be able to include pure contextual
signal pairs that lack marginal TSP association.

**Gate 3.** Actual small trees of depth=1 first. Prefer
hierarchical shrinkage of the 2-5 leaf coefficients to
the global parent rather than fitting independent
high-variance local full K10 RLRs on tiny nodes.
Compare 1,2,3,5 per-leaf rank relations as a
predeclared model complexity study only after
baseline and synthetic null validation are stable.

**Gate 4.** True train-only D/V or nested CV acceptance,
null-noise and global-main-effect negative controls
at multiple sample sizes and realistic class imbalance;
false split frequency and benefit vs loss assessed on
unseen synthetic replications.

**Gate 5.** Compare with exact original RLR, global
calibration, pure TSPDT, and a non-tree **single global
rank interaction** logistic control. Without the latter,
performance might come from interaction features
rather than the *tree* structure itself.
Include model size, split / selected gene stability,
heldout LL, Brier, Macro-F1 and balanced accuracy.

**Gate 6.** Validate in a newly acquired untouched
external cohort/laboratory batch with freeze BEFORE
seeing outcomes. All 19 previous cohorts and their 95
splits have been repeatedly exposed to method
development, so another 19×5 study is NOT confirmatory.

## Literature / novelty check

Model-based recursive partitioning (MOB) already uses
parameter score instability, including slope changes,
to select split variables. A tree-based varying-
coefficient model (2025) also exists. TSP decision
trees were published by Czajkowski/Kretowski (2011).
Accordingly, **do not present conditional coefficient
trees or pair rules as automatically novel.**
Potential contribution is a distinct structured
pipeline: original RLR leaf expert + genuine
context-conditional within-patient ordinal gene pair
interaction + hierarchical compactness and rigorous
selection null/fallback, proven on new data.

References:
- Zeileis, Hothorn, Hornik, "Model-Based Recursive
  Partitioning" (2008), https://www.zeileis.org/papers/Zeileis%2BHothorn%2BHornik-2008.pdf
- Czajkowski, Krętowski, "Top Scoring Pair Decision
  Tree for Gene Expression Data Analysis" (2011),
  DOI:10.1007/978-1-4419-7046-6_3
- Zakrisson & Lindholm (2025),
  "A tree-based varying coefficient model",
  https://doi.org/10.1007/s00180-025-01603-8
