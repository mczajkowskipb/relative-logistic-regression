# RLMT-CA-PS-v1: class-aware, parent-shrunk relational logistic model tree

**Research prototype; built 2026-10-21 after inspection of the 19 published
RLR development cohorts. NOT prospective/pre-registered evidence.**
Model identity: `RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2`.
New branch `feature/rlmt-class-aware-ps-v1`, isolated from published
RLR `main=4079a5ac400f1d2cfded523c7000f655470925ce`
and the 285+95 frozen model results.

## Motivation and IMPORTANT guardrail

A classifying split need NOT leave every child in the 50:50 class
ratio. Useful trees may deliberately enrich one child for a disease
subgroup. Forcing 50:50 can destroy informative rules and change
the target scientific question. Our intervention is to **require
both classes to be estimable** in each child, discourage highly
one-sided nodes when fitting a local logistic model, and give small
or imbalanced children **less authority to override the parent**.

Historical standard TSP scores already compare separate
class-conditional frequencies. The present CA-PS test is NOT
the 2008 historical Weight-kTSP or the distinct Luo
cost-sensitive WTSP. See
[CLASS_IMBALANCE_TSP_WTSP_RLMT_V2_RESEARCH.md](CLASS_IMBALANCE_TSP_WTSP_RLMT_V2_RESEARCH.md).

## 1. Data hierarchy and separation

On every outer-training set T:

1. Build or receive the **T-only** original RLR relation pool;
   no outer test labels.
2. Partition T with the original deterministic RLMT-DC
   discovery/confirmation split D/V.
3. Inside D, rebuild supervised RLR top1000 relation pools
   independently inside each discovery OOF training fold
   via unchanged `rlmt_dc_parent_oof_strict()`.
4. Reuse 4096 deterministic unordered discovery pair proposals.
   Rank only the feasible candidates as specified below.
5. Nominate **one** candidate using D labels and D-only
   out-of-fold parent probability residuals.
6. On V, use the existing D-fitted parent and split leaf models
   to obtain two probabilities; mix their predictions using
   **weights computed from D labels only**. Confirm the
   nominated candidate using paired validation logloss
   `mean(parent_loss - mixed_loss) - SE > 0`, min nV=20.
   No second candidate is selected based on V outcomes.
7. On full T, fit the **original unmodified RLR parent**
   with its original seed and hyperparameters. If the split
   was not confirmed, its test probabilities must match that
   parent bit-exactly.
8. If confirmed, refit both local RLR leaves on full T,
   using only T's original relation pool and unchanged
   leaf `glmnet::cv.glmnet` fitting code.
   Calculate NEW shrinkage weights using **full T class
   support** and convexly combine the T-fitted parent
   and child predictions. Outer-test labels enter ONLY
   post-fit evaluation, not weights or candidate nomination.

**Known design limitation:** The D-fitted confirmation parent
uses the original DC's strictly stratified inner CV and
discovery-specific supervised pair pool, while the full-T
RLR reference uses the original publication's CV procedure.
The D->T estimator/pool shift is still present. CA-PS-v1
tests a *conservative mixture* and improved class-aware
screening, **not a proven complete solution to H1**.
A fully aligned estimator-selection pipeline remains a
separate follow-up methodological hypothesis.

## 2. Training-only class-aware split ranking

Retain the frozen original min child support:

`n_child>=16; n_child,class0>=6; n_child,class1>=6`.

The original residual-contrast screen is

`S_residual(z) = |mean(y-p_oof | z) - mean(y-p_oof | not z)| * sqrt(nL*nR/n)`.

Define node class balance (a reliability proxy) as

`B(n0,n1)=4*n0*n1/(n0+n1)^2`.

Then rank feasible candidates by

`S_CA(z) = S_residual(z) * sqrt(min(B_left,B_right))`.

This factor lies in [0,1] and is equal to one for both
50:50 child compositions; it does NOT REQUIRE 50:50
and does not modify the final class decision threshold.
Candidate ties resolve by the deterministic original
`pool_index`.

**Non-claims:** This is a newly proposed reliability
heuristic and may suppress useful disease-enriched
splits. The factor is not a formally derived confidence
interval for a TSP rule and has not been validated
on independent clinical cohorts. The existing
min-per-class support condition is retained unchanged.

## 3. Parent-shrunk RLR probability prediction

For leaf l with full-training counts n0_l,n1_l and
n_l=n0_l+n1_l, define a predetermined, data-support-only

`g_l = B(n0_l,n1_l) * n_l/(n_l+2*K), where K=10`.

Predict

`p_CA(x) = (1-g_l)*p_RLRparent(x) + g_l*p_RLRleaf_l(x)`.

Properties:

- Every probability remains in [0,1] if both inputs do.
- A small or single-class leaf gets limited authority.
- For a balanced leaf with n=20 and K=10,
  `g=0.5`; with n=100, `g≈0.8333`.
- If no tree split is accepted, **p_CA=p_parent bit-exactly**.
- This shrinks probabilities TOWARD the parent; it
  **does not** guarantee calibrated probabilities,
  logloss improvement, or protection if parent itself
  predicts extremely incorrectly.
- K=10 is inherited unchanged from original RLR,
  not chosen by optimizing the examined breast2
  outer test. The 2*K form is a **heuristic**
  requiring independent study.
- The regularized leaf fit (alpha0.5, lambda.min)
  remains the original RLR-style `glmnet` leaf:
  **no class weights or additional probability clipping**.

## 4. Exact implementation boundaries

Files:

- `R/experimental/rlmt_ca_ps_v1.R`:
  complete research fit, D/V candidate selection,
  full-T parent and child training, prediction, metrics,
  and provenance files.
- `R/experimental/test_rlmt_ca_ps_v1.R`:
  unit tests plus mock-callback synthetic end-to-end
  acceptance, rejection and outer-test-label permutation.
  Test doubles are clearly distinguished from real fits.
- `R/experimental/test_rlmt_ca_ps_v1_real_synthetic.R`:
  true unchanged SWAP/switchBox pool, glmnet RLR
  and leaf fits on **random synthetic numeric genes**.
  Reports acceptance/fallback and direct probability
  mixture consistency. No ARFF or real patient data.
- `scripts/rlmt_ca_ps_v1_synthetic_gate.sh`:
  isolated, fail-closed unit+real synthetic gating and
  R source parse.

Each run saves its own original RLR `rlr_reference`
subdirectory, complete leaf model RDS (if accepted),
per-sample parent/leaf/g/mixed probability file, split
candidate audit, discovery/confirmation sample IDs,
confirmation gain metadata and final metrics. Never writes
the original 285/95 results or experimental branch.

## 5. Test status and interpretation

### First actual Ubuntu gate: observed FAILURE, mock infrastructure only

The first synthetic gate was executed at
`cc40166f85b2b6ed534aa1560c3f30d2d8b6439d`.
It passed `RLMT_CA_PS_R_PARSE=PASS` and **16** primitive
arithmetic/candidate/confirmation checks, but the mocked
end-to-end test stopped at:

```
Error in ensure_dir(out_dir) : could not find function "ensure_dir"
```

The source under test calls `ensure_dir()` from the original
`07_world_protocol_utils.R`, while the isolated mock test
**deliberately does not source** that file to avoid loading
production packages. Source review also identified that the
original `rlmt_dc_defaults()` requires its `%||%` operator
from the same file.

Both **pure helpers have now been explicitly supplied as
test doubles in the synthetic mock test**. The main
`R/experimental/rlmt_ca_ps_v1.R` algorithm source was not
changed by this repair.

The failed result directory
`results/rlmt_ca_ps_v1_synthetic_gate/` must remain
untouched. The corrected attempt uses a separate
`results/rlmt_ca_ps_v1_synthetic_gate_v2/` directory.
This is **attempt v2** for the same algorithm v1, not an
algorithm redesign or a change in train/test protocol.

**Status after this patch: the amended R test and real
SWAP/glmnet smoke have NOT YET been executed on Ubuntu.**
A commit alone is never proof of passing the development gate.



At GitHub commit creation:

- **Implemented:** R code, static-bracket structural check,
  mock integration tests, and real integration smoke scripts.
- **NOT YET VERIFIED:** actual R parser, real R tests
  and returned empirical outcomes. The assistant runtime
  does not include Rscript; these gates must be run
  with the user’s Ubuntu R environment and results
  preserved. A committed test is not a passed test.
- **No 19-dataset results** for CA-PS-v1; no original
  methods retrained and no new patient data accessed.
- **No new empirical superiority claim** is authorized.

Next gate after a successful synthetic test should
be a separately locked, new public omics cohort with
verified independent patients, fixed class outcome,
grouped splits, primary unweighted logloss and
calibration secondary metrics, and a direct comparison
with frozen RLR/RLMT-DC baselines. The already
exposed 19 published ARFFs are only DEVELOPMENT
data; any retrospective v1 comparison there must
be labeled explicitly as such.


## Verified Ubuntu v2 synthetic gate — passed

Operator provided the actual Ubuntu output at
`cc40166…` -> `2a1271d0e21377a084fe1056809b2b5beaf76ff9`.
The corrected v2 gate ran successfully, showing:

`RLMT_CA_PS_R_PARSE=PASS`  
`RLMT_CA_PS_SYNTHETIC_UNIT_AND_MOCK_INTEGRATION=PASS`  
`RLMT_CA_PS_REAL_SYNTHETIC=PASS`  
`RLMT_CA_PS_DEV_GATE=PASS`

The true unmocked SWAP/glmnet run on 42 synthetic features,
180 training and 64 holdout samples produced
`ACCEPTED_SPLIT=FALSE`,
`REALIZED_LEAVES=1`,
`HELDOUT_SYNTH_LOGLOSS=0.56640815`.
These outcomes **prove the real no-split pipeline runs**, but
**not** the real accepted two-leaf route.

## Next LARGE reproducibility gate — sources committed, outcome pending

The fixed v3 development gate now groups *all* of the following:

1. R parse and original synthetic unit/mock tests.
2. Unmodified real RLR/switchBox/glmnet fallback smoke.
3. **Actually fitted** 2-child RLR/glmnet, using a
   **TEST-ONLY forced confirmation Boolean** so the
   complete accepted branch is deterministically
   exercised without tuning threshold or cherry-picking
   a synthetic case. All test-side probabilities are
   **independently replayed** from saved original RLR
   parent and both new child model RDS objects.
   A permutation of *only synthetic test labels* must
   not affect nomination, confirmation nor predicted
   probabilities. Forced acceptance is visibly labeled
   and is **NOT** a true confirmation decision, nor
   proof of scientific improvement.
4. Three fixed synthetic data-generating scenarios:
   `balanced_heterogeneous`,
   `minority_heterogeneous`,
   `balanced_null_tree`; two predeclared
   seeds 1139 and 2267, 240 training and 96
   independently simulated test samples per case.
   Each case uses **real, unchanged** RLR source
   for its parent/comparator, and the NEW CA-PS
   method for fitted/predicted probabilities.
   All **six** cases must finish; no adaptive choice of
   favorable seeds, no performance-dependent early stop.
5. Independent Python audit of all six saved model
   prediction pairs: patient/sample identity alignment
   (synthetic IDs), probability mixing identity,
   realized split state, exact original-parent parity
   if no split, raw heldout logloss/Brier/class accuracy/
   balanced accuracy and scenario aggregate deltas.
6. New durable result root
   `results/rlmt_ca_ps_v1_synthetic_gate_v3`,
   unique per-run code SHA, original paper SHA and
   fail-closed status marker. V1 failure and V2
   PASS result roots remain untouched.

**Scientific caveat:** Simulation cannot prove a
calibrated clinical model. Both the shrinkage weight
and class-balance screening factors remain ad-hoc
fixed assumptions introduced after inspecting
the previously used 19 publication datasets.
Do not use even a favorable v3 synthetic result as
independent biomedical generalization evidence.

**Status:** v3 R/simulation tests have **not**
been executed on the operator's Ubuntu host.
No original 19 ARFF / clinical patients or
publication outputs are read by any v3 synthetic test.
Any unsuccessful gate must be debugged without
overwriting its source SHA and failure evidence.

## Benchmark next stage — only after v3 PASS

After v3 succeeds, develop **one complete experimental
runner** that matches all 95 original frozen outer
splits and read-only RLR/RLMT baseline predictions.
For source-matched prior cohorts this is still
EXPLORATORY, not new validation; no choice of
weights or thresholds based on breast2 outcomes.
For independent scientific evaluation, additionally
register new biologically defensible cohorts and
patient-grouped train/test splits **before** looking
at model outcomes.
