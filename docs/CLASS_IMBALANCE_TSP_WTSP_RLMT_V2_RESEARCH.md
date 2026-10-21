# Class imbalance in TSP / Weight-kTSP / RLMT-DC: source-grounded research decision

**2026-10-21 — SCIENTIFIC DESIGN ONLY. No v2 implementation, refitting, retuning,
or access to new cohort outcomes in this change.** The historical publication
`main` and the frozen 285 RLR/RLMT + 95 TSPDT-inspired results are untouched.

## 1. Distinct source methods that MUST NOT be conflated

### Classical TSP — Geman et al. (2004)

For a pair (i,j), `Z_ij=I(x_i<x_j)`, with binary classes c=0,1,
class-specific empirical frequencies `q_c=P_hat(Z_ij=1 | Y=c)`:

`S_TSP(i,j)=|q_1-q_0|`.

The denominators `n_0` and `n_1` are *already different*, so this
score does **not** naively favor the numerical majority merely by using
more majority-class samples. It can still be unstable because the minority
`q_c` has higher sampling error; it does not solve all imbalance problems.
https://pmc.ncbi.nlm.nih.gov/articles/PMC1989150/

### Weight-kTSP of Czajkowski & Krętowski (IEA/AIE 2008)

https://www.aragorn.wi.pb.edu.pl/~mczajk/ieaaie2008.pdf
DOI: 10.1007/978-3-540-69052-8_48

**This is NOT the same as cost-sensitive WTSP**. The 2008 paper's primary
modification is a *relative ratio-to-estimated reference* relation:
`I[(x_i/x_j)<S_ij]`, rather than strict `I[x_i<x_j]`;
the secondary ranking and k-pair voting are also modified (paper §§2.3).
`S_ij` is based on training-only average expression ratios in its
specified procedure. This adds a train-dependent threshold and may
uncover non-reversing but class-discriminative relations. It is not an
automatic inverse-class-frequency-weighted TSP.

**Compatibility hazard:** ratios/relative percentage changes are ill-defined
or distorted on logged, zero, negative or baseline-shifted expression scales.
Do not blindly apply the historical 2008 ratio formula to pre-log-transformed
ARFF without verifying measurement scale and positivity.
If a new log-contrast threshold is proposed, it is a **new method**, not exact
historical Weight-kTSP.

### Cost-sensitive WTSP of Luo et al. (PRIB 2008)

Huaien Luo, Yuliansa Sudibyo, Lance D. Miller, R. Krishna Murthy Karuturi,
*Weighted Top Score Pair Method for Gene Selection and Classification*,
LNCS 5265, 323–333; DOI: 10.1007/978-3-540-88436-1_28.
Public accessible full text:
https://www.researchgate.net/publication/221336115_Weighted_Top_Score_Pair_Method_for_Gene_Selection_and_Classification

Their **distinct** optimization is class-prevalence/cost-aware:
`C = pi_1 * lambda_1 * P(error|1) + pi_0 * lambda_0 * P(error|0)`.
With an oriented pair classifier, one of the published weighted rank
forms is equivalent to preferring
`pi_1*lambda_1*q_1 - pi_0*lambda_0*q_0` (for the corresponding
orientation). Consider both pair orientations and use the actual
training-only error counts for the cost definition; label and
tie conventions must be explicit.

At effective costs `pi_1*lambda_1=pi_0*lambda_0`, the contrast is
a scaled classical TSP criterion. Equal misclassification costs
`lambda_1=lambda_0` instead optimize the empirical
prevalence-weighted **classification error**, not balanced error.

### TSPDT 2011

https://wi.pb.edu.pl/p/wp-content/uploads/sites/3/2024/07/AEMB2011.pdf
DOI: 10.1007/978-1-4419-7046-6_3

The article specifically distinguishes TSPDT, k-TSPDT, Weight-TSPDT
and Weight-k-TSPDT. Root test matches its corresponding TSP-family
algorithm; rules are recomputed on node-local samples and C4.5-like
pessimistic pruning is used. Experiments include supervised Relief-F
f1000 (ten neighbors), not the isolated v1's MAD100. The paper refers
to the authors' Weight-kTSP [8], NOT necessarily cost-sensitive Luo
WTSP. Do not rewrite the historical 2011 identity.

## 2. Three separate imbalance questions

**A. Pair/test selection at a node.**
If targeting balanced classwise error, classical TSP conditional
score is already class-normalized. For tiny minority nodes, add
training-only uncertainty/stability screening (e.g. conservative
estimates for class-conditional pair probabilities; independent
minimum per-class support) or a prespecified cost-sensitive
pair criterion with explicit clinical costs. Do **not** insert an
unmotivated factor 1/n_minority in the original difference;
it double-counts class conditioning.

Cost-aware pair selection may help accuracy under known deployment
class prevalence and misclassification costs. It does NOT
directly optimize probabilistic logloss, nor does it guarantee
better minority recall.

**B. Node feasibility and statistical support.**
Record (n0,n1), pair-route counts, minority support in BOTH
children, relative minority prevalence, effective sample
support, selection stability and internal train-only prediction
uncertainty. A minimum total leaf size alone is insufficient.
But a hard minority threshold may exclude truly useful splits
and can defeat original root-TSP equivalence; historical
TSPDT reconstruction must keep this distinction explicit.
For rare classes, the main issue may be lack of enough
observations to justify any distinct leaf model.

**C. Probabilities from a logistic RLR leaf.**
This is the primary documented RLMT-DC-v0 failure mode:
9/95 accepted splits; 6/9 worsened independent outer logloss,
mean across 19 datasets +0.017217934; worst breast2 fold01
+1.698460, where two leaves (26 and 35 fitted samples) both
contributed to severe highly confident wrong predictions.

Blind inverse-frequency weighting / rebalancing of glmnet's
binomial loss changes the *estimated posterior* and can worsen
calibration and logloss even if balanced accuracy improves.
See van den Goorbergh et al., *The harm of class imbalance
corrections for risk prediction models* (2022):
https://pubmed.ncbi.nlm.nih.gov/PMC9382395
Publication-case-control class prevalence may not equal
real clinical prevalence. Treat these as different estimands.
Population prior correction requires a justified, externally
known deployment prevalence; do not infer it from an arbitrary
case-control dataset:
https://gking.harvard.edu/publication/inference-in-case-control-studies/

## 3. Recommended ONE first successor method, not a combinatorial grid

**Provisional identity: RLMT-DC-PS-v1 (parent-shrunk logistic leaves).**
A separate *new* method proposal, **not an evaluated improvement**.

First retain the same high-level data-separation principles,
but avoid learning two fully independent local predictive models
without an explicit shared baseline. Propose:

`p_PS(y=1|x,r) = sigmoid[logit(p_parent(x)) + delta_r(x)]`,

where r is the test routing result and `delta_r(x)` is
an explicitly regularized, *small* train-only residual logistic
correction. Strongly penalize both local intercept changes
and local pair coefficients when node support is limited.
Alternatively compare a **proper convex mixture**
`p_mix = (1-g)*p_parent + g*p_leaf` with a **single**
shared gating/shrinkage parameter `g in [0,1]`, learned
strictly by nested internal CV. The mixture is a simpler
implementation reference, but its mixing parameter must not
be chosen from the 19 exposed test fold outcomes.

Choose **one** of these mathematical forms as the v2 primary,
then specify intercept/regularization treatment and complete
training procedure. A priori scientific preference is the
shared-parent logit residual form because it is parsimonious,
but complexity may justify testing the simpler mixture only
on synthetic data first. Neither is guaranteed to improve
logloss or represent original TSPDT.

Crucial invariant: `delta_r(x)=0` or `g=0` must
return **exact** original-parent predictions; no accepted split
must return exact published RLR output (including sample IDs,
factor ordering, 1e-12 probability tolerance).

**Confirmation-model alignment.** The original RLMT-DC
confirmation estimates gain for child models trained on
discovery D with D-specific pair pools, but final accepted
models are refit on outer training T with a different
full-T supervised pair pool and child re-ranking. This is
valid train/test separation but the confirmation estimator
is not the identical final fitted model. A separate method
must internally validate the **complete training procedure**,
not only a surrogate: within each inner train partition
repeat all supervised pair selection, parent fitting,
local fitting, and internal lambda selection; assess on
disjoint validation samples. Nominate candidates without
using confirmation labels; use confirmation to assess only
the nominated fixed candidate. A positive 1-SE margin is
not a generalization guarantee.

**Class imbalance integration:** by default retain
**unweighted proper logloss** if the objective is calibrated
patient-level probability on the sampling target. Do not
introduce cost-sensitive pair selection and class-weighted
logistic leaf training in the same unidentifiable change.
Class weights, cost-weighted pair selection and a changed
clinical decision threshold are scientifically distinct
comparators that require separately justified estimands.
For case/control sampling, report balanced accuracy and
macro-F1, but do not claim population-risk calibration
without known sampling-to-population adjustment.

## Existing tooling to reuse, not reinvent

- **Relief-F:** R `CORElearn::attrEval` explicitly supports
  ReliefF and cost-sensitive ReliefF; investigate the exact
  10-neighbor estimator configuration, tie and missing-value
  conventions on synthetic fixtures rather than silently
  treating it as identical to historical Weka Relief-F:
  https://search.r-project.org/CRAN/refmans/CORElearn/html/attrEval.html
- **Penalized logistic offsets:** `glmnet` supports fitting
  with `offset` and predicting with `newoffset`;
  this offers a path to a parent-logit offset plus sparse/ridge
  local correction without inventing a new logistic optimizer.
  Check `cv.glmnet` fold/offset semantics and train-only
  selection before implementation:
  https://github.com/cran/glmnet/blob/master/R/predict.glmnet.R
- **Class-weighting:** `glmnet` observation weights are
  available but should be a separately justified experiment,
  not an automatic remedy for biased sample prevalence.
- **Tree pruning:** a historical 2011 reproduction should
  follow the relevant documented C4.5-like method
  rather than silently substituting contemporary CART
  cost-complexity pruning or fixed max-depth.

## 4. Pre-implementation and prospective evaluation gates

1. Write a frozen mathematical contract: objective,
   prevalence meaning, candidate pair criterion,
   tie convention, fold hierarchy, feature provenance,
   shrinkage estimator, class-order handling and
   hard-class decision threshold.
2. Synthetic-only controlled counterexamples for:
   strong class imbalance, tiny minority support in
   a child, constant/degenerate pairs, complete separation,
   extreme wrong probabilities, and confirmed candidate
   shifting its refit behavior. Ensure stable finite
   probability output and deterministic replays.
3. Unit test classical TSP ranking against reference;
   test Luo cost-sensitive pair ranking separately;
   reconstruct the historical 2008 Weight-kTSP **only
   on legitimate positive raw-expression data**.
4. Run an isolated, train-only simulation suite with
   independent synthetic test data and varying class
   prevalence. Track true expected logloss/Brier,
   calibration-in-the-large, per-class recall,
   balanced accuracy and accepted-tree complexity.
5. Keep previously audited 19 publication ARFF results
   frozen and labeled **exposed, exploratory**. No
   choosing v2 priors/weights/lambda or stopping
   rules by optimizing on breast2 or leukemia.
6. Identify an actually independent public clinical
   or omics cohort with genuine new patients and
   defensible binary end point; freeze source license,
   class labels, patient grouping, preprocessing,
   target sampling prevalence and outer splits
   **before** any v2 outcome inspection.
7. Compare new v2 against unmodified RLR and the
   already frozen RLMT-DC-v0. Primary **unweighted
   logloss** unless the scientific use case
   explicitly predefines different weighted risk.
   Prespecify Brier, per-class recall, balanced
   accuracy, macro-F1, calibration and interpretable
   model complexity. Report every cohort and failure.

## Priority and explicit NON-decisions

**P0:** No change to v0 and no "fix breast2" by new
post-hoc weights or a stricter numerical 1-SE threshold.

**P1:** Verify implementation of 2011 TSPDT separately:
historical 2008 Weight-kTSP is NOT Luo WTSP; do
not silently relabel historical TSPDT as cost-sensitive.
Historical fidelity requires train-only Relief-F1000,
root pair equivalence, tie rank, proper pruning and
scale/ratio checks for Weight-kTSP.

**P2:** Implement one parsimonious parent-shrunk RLMT
successor with inner validation aligned to the
full refitting mechanism, then synthetic-only evidence.

**P3:** Test genuinely cost-sensitive pair selection
only when the deployment class prevalence and clinical
misclassification trade-off are clearly specified.
If those are unknown, use classical class-conditional
TSP score as the clean ranking reference and report
classification/probability objectives separately.

No algorithm was modified or external patient outcome observed
by creating this research note.
