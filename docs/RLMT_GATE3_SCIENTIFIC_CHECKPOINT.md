# RLMT GATE-3 scientific checkpoint — 2026-10-21

**Status:** technical gates passed; predictive superiority **not** demonstrated.

This records the *observed* GATE-3 pilots and the rules for subsequent
research. It is not a methods approval or a change to the frozen RLR paper.

## Immutable code and experiment contexts

- RLR publication `main`: `4079a5ac400f1d2cfded523c7000f655470925ce`.
- Original RLMT GATE-3 (default random internal `cv.glmnet` fold allocation):
  `results/rlmt_gate3_pilot_2x3`.
- RLMT GATE-3 with deterministic class-stratified *leaf-internal* glmnet CV:
  `results/rlmt_gate3_stratcv_2x3`.
- Code version of rerun (recorded in its `CODE_SHA.txt`):
  `7b9f716c932a65fdd359c1f5e83facae13b3ff73`.
- Both pilots: datasets Chen-2002 and Singh-2002, first 3 outer folds from
  one 5-fold stratified split, seed 42; fixed K10/m10/d2/pool1000,
  alpha=0.5, lambda.min, top B=20 and conservative 1-SE pruning.

All result paths are local to the compute server and are **not** tracked as
GitHub experimental outputs. Reported findings below are from terminal
audits of those local outputs.

## Technical audit results

| Measure | Original pilot | Stratified leaf-CV rerun |
| --- | ---: | ---: |
| Completed outer evaluations | 6/6 | 6/6 |
| Screened candidate assessments | 120 | 120 |
| Valid assessments | 96 | 116 |
| `glmnet` class-absence errors | 17 | 0 |
| `cv.glmnet` dimension errors | 3 | 0 |
| Inner-child support failures | 4 | 4 |
| Candidates with positive mean gain | 4 | 10 |
| Candidates with positive lower 1-SE gain | 0 | 1 |
| Genuinely accepted two-leaf splits | 0 | 1 |

The three dimension errors were independently reproduced and resolved by
explicit stratum-balanced `foldid`, with unchanged input matrices.
RLR-to-RLMT-L1 parity, frozen RLR outputs across pilots and outer-test
logloss/Brier reconstruction passed on both datasets in all folds.

One accepted split occurred in Singh-2002, `rep00_fold00`, on
`f_216_at > f_39758_f_at`.

- Inner OOF mean logloss gain = `+0.13849446`; SE = `0.12668074`;
  lower-bound gain = `+0.01181372`.
- Outer-test RLR logloss = `0.29002105`; RLMT-L2 = `0.51731647`;
  `RLMT - RLR = +0.22729542` (worse for RLMT).
- Outer-test RLR Brier = `0.07592205`; RLMT-L2 = `0.08449060`.
- Original class labels: `N` (negative), `PR` (positive).
- Child node 2 training: 31 N / 8 PR; test: 10 N / 2 PR.
- Child node 3 training: 9 N / 33 PR; test: 0 N / 9 PR.
- Child 2 per-sample-average change in test logloss = `+0.43717757`;
  child 3 = `-0.05254745`.
- On the entire outer test, RLMT improved individual logloss for 19/21
  samples and worsened it for 2/21; both models had 2 classification errors.
- `sample_92`, true PR and routed to child 2: RLR P(PR)=`0.032270`;
  RLMT P(PR)=`0.000193`; loss change `+5.118134`.
- `sample_80`, true PR and routed to child 2: RLR P(PR)=`0.574324`;
  RLMT P(PR)=`0.205025`; loss change `+1.030061`.
- Excluding `sample_92` produces a descriptive leave-one-out average
  change `-0.017247` over the other 20 observations. **This exclusion is
  not an admissible change to the reported evaluation.**

These observations suggest instability in minority-class probabilities in
a small leaf; no single fold establishes that the method is generally
inferior or superior.

## Methodological limitations

1. The global top-1000 pair pool is label-supervised and computed once on the
   entire outer-training split. The inner OOF folds therefore reuse a
   supervised representation influenced by the inner validation labels.
2. Candidate screening also uses residuals derived from the same inner OOF
   observations later employed to compare the shortlisted candidates.
3. Choosing the best of up to 20 candidates based on its estimated
   `mean_gain - SE` can introduce selection optimism. The **1-SE rule is a
   pre-pruning heuristic, not a family-wise hypothesis test** or a
   statistical guarantee of positive improvement on unseen data.
4. The outer-test is independent of the fit for its own fold, but
   `singh-2002/rep00_fold00` has now been examined in detail.
   **No threshold, feature, calibration, shrinkage or pruning tuning using
   those outer-test outcomes may be presented as independent confirmation.**
5. Stratified glmnet CV fixed the documented class-fold failures but did not
   independently establish well-calibrated class probabilities.

## Decision for the next research gate

**GATE 3 technical: PASS. GATE 3 scientific superiority: NOT ESTABLISHED.**

Freeze these pilots and the corresponding code SHA as *hypothesis-generating
evidence*. Do not modify them, alter the historical 1-SE outcomes, discard
the two difficult samples, or relax pruning.

Before any further results are inspected:

1. Document a fixed evaluation protocol, complete code commit, eligible
   external datasets and audit criteria. Do not pick datasets/splits
   based on favorable RLMT performance.
2. Clearly separate a **frozen original RLMT method** (for apples-to-apples
   reproduction) from any **new candidate algorithm** for honestly
   evaluating split-selection performance. The latter may require
   independently reconstructing pair rankings and screening *within each
   training partition*, rather than reusing them across validation.
3. If sample size makes such nested discovery/confirmation unstable,
   treat this as a limitation and report the evidence; do not lower
   class-support rules merely to force additional trees.
4. Measure outer-test logloss (primary), Brier, AUC, F1, realized leaf
   counts, selected-candidate gain versus independent gain, class support
   in leaves, warnings/failures and probability tail risk.
5. All future candidate tuning motivated by this observed Singh holdout
   must be assessed on **new, previously uninspected outer-test cases or
   genuinely independent datasets**. The remaining two folds of these
   same datasets may be reported as exploratory continuation only:
   their train/test partitions overlap data used during method development.

No new algorithm or prospective confirmatory benchmark is authorized by
this checkpoint alone; it freezes what is known and records the scientific
protocol constraints.
