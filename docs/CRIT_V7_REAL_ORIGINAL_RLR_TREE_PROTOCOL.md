# CRIT-v7 — ORIGINAL_RLR anchored real conditional model tree

**Status 2026-10-23: REAL R implementation and 2×3 guarded experiment
prepared on GitHub; NO R execution or biomedical improvement yet confirmed.**
Repository: mczajkowskipb/relative-logistic-regression
Branch: feature/crit-v7-regime-switch-probe
Scientific stage: EXPOSED_RETROSPECTIVE_DEV_ONLY

## Why this particular tree

The full independently audited 19 cohorts ×5 folds negative evidence:
- v4 ACTIVITY_ROUTE +0.001233117 logloss vs original paper RLR;
  SOFT_GAP_ROUTE +0.050673367;
- v5 honest router 0/6 accepted, exactly original RLR;
- v6 ONEPAIR_MARGIN +0.005470103 (8/11 cohort wins/losses),
  ALLPAIR_MARGIN +0.010726693.
All development data exposed. Global additive rank-margin repair
was not a robust solution. The user has experience with linear model
trees whose leaf models generally use only 2–5 features; start with ONE
conditional rank comparison per leaf, NOT 10 separate leaf RLR fits.

Novel falsifiable mechanism: a relation may REVERSE its predictive
effect depending on a second within-patient relation, despite having
little/no marginal score. Need an actual GENE-ORDER split, not an
arbitrary threshold on a root probability.

## New REAL R algorithm

Source: R/experimental/crit_v7_original_rlr_two_leaf.R

Root is the LITERAL ORIGINAL published RLR:
switchBox top1000, K=10, m=10, d=2,
cv.glmnet alpha=.5, lambda.min, seed42. Reference predicted
probabilities are independently compared against the 285-model
archived original-RLR patient-level split with tol 1e-9.
RLR_ALIGNED_V2 IS NOT USED AS THE BASELINE.

A confirmed tree:
  z(x)=I(geneA>geneB)
  logit P_CRIT(x)=logit(P_ORIGINAL_RLR(x))
          + (1-z(x))*(beta_false+gamma_false*r_false(x))
          + z(x)*(beta_true+gamma_true*r_true(x))
with r_leaf(x)=+1 if geneC>geneD, otherwise -1.
Each leaf may select its OWN distinct gene pair.
Fixed ridge intercept=.10, leaf coefficient=.20.
Predicted full logit offset clamp ±2.5.
If not confirmed: **EXACT original RLR**.

Supervised steps occur only inside OUTER TRAIN:
1. split T into stratified discovery D=2/3 and heldout
   confirmation V=1/3 using seed 773, reproducible;
2. build original RLR parent on D, and 3-fold original
   parent OOF probabilities STRICTLY from D (not V);
3. construct only D-based pair pool:
   up to 66 label-free rank-variability-selected pairs
   from 12 genes, plus original RLR selected pairs,
   plus 256 deterministic randomly proposed pairs.
   This avoids excluding pure XOR interactions lacking
   a high marginal TSP score. NO OUTER TEST patients
   or labels influence any candidate;
4. residual interaction screen on D using
   |sum_D (y-p_D_OOF) * gate * candidate_relation|,
   requiring disjoint gate and feature genes, >=10 samples
   and >=2 per class in each leaf. Up to 6 candidate gates.
5. within each gate, independently select 1 leaf model
   relation in each branch by centered local residual
   correlation, requiring both signs >=3 in the branch.
   Fit ridge logistic leaf-specific intercepts and
   slopes with strict original-RLR D-OOF logit offset;
6. obtain parent_D predictions on disjoint V, compare
   fixed tree predictions to parent_D by per-patient
   paired heldout logloss. Conservative gain is mean gain
   minus 1.0 SE(gain) minus complexity cost 0.002.
   Accept a split only if this exceeds 0.002.
   The best D/V candidate may still suffer winner's
   curse due to multiple V-evaluated gates; report it.
7. refit selected tree coefficient vector on all OUTER
   TRAIN with strict 3-fold original RLR OOF probabilities,
   keeping selected gate and leaf pairs fixed.
   The exact original parent itself is independently
   trained on all T exactly as in the original publication.

No tiny-leaf independent cv.glmnet, no mutation of historical
RLR artifacts. The code may refuse trees if genuine strict
OOF or minimum class support cannot be satisfied.

## Predeclared four-way OUTER HELDOUT comparison

- ORIGINAL_RLR: exact published root;
- GLOBAL_CAL: root + training-only intercept calibration;
- GLOBAL_INTERACTION: explicit non-tree single global
  gate×relation logistic interaction using the same train
  selected gate / D-scored candidate; trained with exact
  original-RLR full train OOF offset, even when the tree
  itself is rejected (when candidate exists);
- CRIT_TREE: real two-leaf, 1 local rank relation per leaf,
  potentially DIFFERENT across the two leaves; only accepted
  with strict D/V evidence, otherwise EXACT ORIGINAL_RLR.

A one-split model with identical relation in both leaves
can be written as a global logistic interaction. Therefore
**a gain over ORIGINAL RLR alone is not a proof that the TREE
representation gives superior functional expressivity**.
Novelty requires comparison with GLOBAL_INTERACTION and
prior model-based recursive partitioning / TSPDT.

## Fail-closed tests and limited pilot

GATE 0: two real synthetic scenarios (planted strong XOR and
label-only global baseline), actual switchBox+cv.glmnet original
RLR parent, explicit independent original paper wrapper parity,
original RLR OOF D/V disjoint, true two-leaf fit and RDS replay,
rank-preserving affine invariance, planted XOR real split
confirmation and perfect-original baseline forced refusal.

GATE 1: previously EXPOSED chen-2002 and gordon-2002,
three original seed42 outer folds each: 6 genuinely fitted
training bundles, 4 prediction methods each=24 saved
heldout method-fold vectors. Frozen original dataset file
source SHA256, patient-ID equality for train and test,
label equality, and archived original root predictions
must match within 1e-9. New experimental output dir only:
results/crit_v7_original_rlr_2x3_v1.

GATE 2: independent Python patient-level replay:
6×4 scores and Brier values, real leaf assignments,
no-split exact-original fallback, named gate+two leaf pairs,
0/1/2 leaf ID validation, archived original RLR parity,
strict no-overwrite and status failure on any discrepancy.

Test code can reject legitimate runs through numerical
issues or class support; no partial pass or silent
substitution is permitted.

## Validation caveats

The entire 19-cohort population and all 95 CV fold locations
have been inspected in previous RLMT development. This pilot
is therefore explicitly **EXPOSED DEVELOPMENT ONLY**.
No implication of statistical significance or an externally
confirmed gain can be made from this 2×3 experiment.
Even a promising six-fold pilot must be followed by locked
19×5 descriptive robustness and crucially an **untouched
external independent patient/cohort/batch confirmation**.
Only if the 2-leaf design exhibits stable benefit vs ORIGINAL
RLR **and GLOBAL_INTERACTION** should depth 2 and 2–5
active relations per leaf be studied.

Also: under a strong XOR planted simulation a simple rank
interaction is expected to beat a global additive model;
that synthetic representational fact is not evidence for
real patient subgroup interactions.
