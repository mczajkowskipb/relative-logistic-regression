# RLMT GATE 4 — Gordon-2002 experimental freeze

**Status:** frozen *before* running GATE-4 RLMT. Technical/exploratory
benchmark only; not independent external cohort validation.

## Eligibility was decided without model metrics

Performance-blind screen:
`docs/RLMT_GATE4_METADATA_RULES.md` and
`R/new/37_rlmt_gate4_metadata_screen.R`.

The 2026-10-21 local audit reports all 9 CompCancer data checksum
matches and absence of nonfinite numeric values / exact numeric
row duplicates. Among the 6 datasets not used in RLMT development,
**only `gordon-2002` passes the prespecified thresholds**.

- Gordon-2002: 181 samples, 1626 numeric genes/features.
- Binary labels: AD (150), MPM (31). Minority class is MPM.
- Canonical processed ARFF SHA256 (publication manifest):
  `42bb0900e7c22851244c27243b626514917c0674d1bbc1`.
- Processed ARFF input: `/opt/2026/RLR/data/CompCancer_truth_9x10/gordon-2002.arff`.
- Patient independence cannot be established from the processed ARFF.

Other eligible datasets **must not** be added after model results
are inspected by relaxing the sample/class thresholds.

## The frozen run

Configuration: `configs/rlmt/gate4_gordon_frozen.json`.

- RLR K10/m10/d2/pool1000, RLMT-L1, RLMT-L2.
- Model options: K=10, m=10, d=2, pair pool 1000, alpha=0.5,
  lambda.min, max_leaves={1,2}, max_depth=2, screening top B=20,
  min_child_n=16, min_child_class=6, inner_folds=3.
- Deterministic class-stratified **leaf-internal** glmnet CV enabled.
- 1 repeat × 5 stratified outer folds; seed=42, same protocol and
  row assignment generator as existing publication RLR code.
- All 5 outer folds `rep00_fold00` through `rep00_fold04`
  are evaluated, without early stopping after seeing an attractive result.
- New result root: `results/rlmt_gate4_gordon_frozen_5fold`.
  Existing experiments are not overwritten.
- The launch must verify branch, publication-main SHA, code SHA, input
  checksum, eligible set, and absence of an existing GATE-4 output.

## Predefined outcomes

1. Primary descriptive outcome: **paired outer-test logloss difference**
   `RLMT-L2 - RLR` for each of all five folds (lower is better).
2. Secondary: Brier score, AUC, F1, realized leaf count, accepted split
   fraction, class-support and error/warning counts.
3. Validity: exact RLR-to-L1 probabilities, class predictions and
   coefficients; unchanged RLR; independent recomputation of metrics
   from the original class labels; model-to-leaf routing and
   serialized model prediction reconstruction for accepted splits.
4. Every split is determined from its own outer-training partition.
   **Never** use outer-test metrics to decide whether its split should
   be retained or what algorithm parameters should be used.

## Interpretive limits

This is *one historical dataset*, not five independent cohorts;
5-fold partitions share observations and do not establish a
cross-dataset statistical superiority claim. The original RLR
publication already used Gordon-2002, but the RLMT algorithm has not
been piloted on it in this conversation.

The current candidate-selection method computes the supervised
global pair pool on the entire outer-training fold and reuses it in
inner CV, and uses inner OOF residuals to screen candidates.
Thus inner-OOF selection gain is not an unbiased nested estimate
of generalization, even when the outer-test is untouched.

Do not change the model, threshold, calibration, eligible dataset
list, or fold selection after analyzing Gordon's outer-test.
A truly confirmatory study requires genuinely independent datasets
and a separately frozen (ideally leakage-robust) selection algorithm.
