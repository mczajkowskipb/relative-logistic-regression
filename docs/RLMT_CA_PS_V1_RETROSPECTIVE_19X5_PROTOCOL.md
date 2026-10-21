# RLMT-CA-PS-v1: complete original-paper 19×5 retrospective comparison

**Scientific status: EXPOSED / RETROSPECTIVE DEVELOPMENT. Not external validation.**
**Data and methods were previously examined in the original 2026-10-21
RLR/RLMT-DC/TSPDT diagnostic program. No prospective improvement claim.**

## Previous reproducibility gates observed (not speculative)

The Ubuntu operator executed full v3 large synthetic test at
`63bb63586434f9e557bb1a1f0bd8ca44a09d3394`. Recorded:

`RLMT_CA_PS_LARGE_R_PARSE=PASS`  
`RLMT_CA_PS_ORIGINAL_UNIT_GATE=PASS`  
`RLMT_CA_PS_REAL_SYNTHETIC=PASS` (naturally no split)  
`RLMT_CA_PS_REAL_FORCED_ACCEPTED_BRANCH=PASS`  
`REAL_PARENT_AND_TWO_LEAF_RDS_REPLAY=PASS`  
`HOLDOUT_LABEL_PERMUTATION_INVARIANCE=PASS`  
`RLMT_CA_PS_SIX_CASE_SYNTHETIC_EVALUATION=PASS`  
`RLMT_CA_PS_SYNTHETIC_INDEPENDENT_METRIC_REPLAY=PASS`  
`RLMT_CA_PS_LARGE_SYNTHETIC_GATE=PASS`  
`LARGE_BATCH_COMPLETE=PASS`.

Six fixed **real** synthetic `glmnet/switchBox` experiments,
three data-generating mechanisms by two seeds, accepted 1/6
splits. Mean paired CA-PS minus original-RLR logloss
**−0.00231333**, Brier from one positive case improved
**−0.007862**; all other five samples of scenarios were
exact parent fallbacks. In the balanced heterogeneity
scenario one of two fits accepted, mean Δlogloss
**−0.006939983** and mean Δbalanced accuracy
**+0.02738359**. Both unbalanced heterogeneity tests
and both null-tree tests did not accept any split.
This is engineering and synthetic plausibility only, NOT
a statistical demonstration of superiority.

## Next complete scientific-data run (code prepared, NOT YET executed)

Objective: **95 NEW** RLMT-CA-PS-v1 train/test fits over
19 PREVIOUSLY USED binary microarray cohorts and all five
original RLR publication folds. As independent comparators,
read **285 ORIGINAL** completed RLR, RLMT-v1 and RLMT-DC-v0
model records from frozen reference output.

- Reference repo:
  `/opt/2026/RLR/relative-logistic-regression`
  at commit `4873878f57f806bf5b52144756f8087b17fa540c`.
- Original publication `main`:
  `4079a5ac400f1d2cfded523c7000f655470925ce`.
- New worktree:
  `/opt/2026/RLR/relative-logistic-regression-rlmt-ca-dev`
  on `feature/rlmt-class-aware-ps-v1`.
- Reference complete benchmark and immutable evidence:
  `REF/results/rlmt_dc_published_data_5fold_v1`.
- NEW private result root (must not previously exist):
  `results/rlmt_ca_ps_v1_retrospective_19x5`.
- 19 source files: SHA256 verified against the exact
  independently audited publication input manifest.
- Same stratified original five folds per dataset
  (seed42, rep00_fold00..rep00_fold04).
- Each original training split produces a fresh
  train-only top1000 RLR pair pool and one NEW
  RLMT-CA-PS fit; **no original 285 method is retrained
  as a comparator**. Note that the new model requires
  its own fresh original-form RLR **parent** fit; this is
  an integral and unavoidable part of the new
  predictor, not an overwrite/retrain of the archived
  original comparator.
- Each new RLR parent must reproduce the saved old
  RLR probabilities to 1e-12 and the exact original
  hard labels for the same sample IDs. Abort on
  mismatch instead of quietly inventing new controls.
- For every **rejected split**, new method prediction
  must be exactly the same as its own RLR parent.
- For every **accepted split**, persist full parent
  RLR model, both refit local RDS objects, confirmation
  metadata, route and convex-shrink mixing components.
- Save 95-row new fitting manifest incrementally,
  preserving any failure evidence and refusing a
  second run in the same output root.
- Independent Python audit rereads all **380 model/fold
  records** (95 × four methods), all 95 paired
  original test sample IDs/labels, factor class-level
  order from archived \`class_support.csv\`, complete
  probability mixtures, original-parent identity and
  all five metrics: accuracy, balanced accuracy,
  macro-F1, logloss and Brier.
- Dataset-mean W/L/T summary uses **19 datasets as
  descriptive units** rather than pretending the
  correlated 95 folds are independent datasets.
- No change to weights, score heuristics, minimum
  leaf class counts, validation threshold, CV split,
  data preprocessing or post-hoc class threshold.

One-command launcher:
`bash scripts/rlmt_ca_ps_v1_retro_batch.sh /opt/2026/RLR/relative-logistic-regression`.
It includes source-hash preflight, real 95-fold execution
and independent audit, with robust stop-on-first-evidence-
inconsistency behavior. The R run may require substantial
CPU time; no background running, auto-install or result
replacement is performed by the protocol.

## Evaluation and criteria for further work

**Primary descriptive outcome:** per-dataset mean
difference of unweighted held-out logloss against RLR.
**Secondary:** Brier, class-balanced accuracy,
macro-F1, hard accuracy, accepted-tree rate, min
per-class leaf support and adverse fold tails.

If all 95 models equal the original RLR, do NOT
call the method superior: record the no-split
behavior. If there is harm or a breast2-like
extreme-logloss outlier, preserve it and do not
change \`g\` or the balance ranking to remove the
failure on these exposed holdouts.

For real publication claims, separately freeze
ONE OR MORE *new* independent patient groups/cohorts,
biological outcome definitions, preprocessing and
group-aware splitting **before** seeing their
held-out metrics. The above 19 datasets cannot
be laundered into a new blind test by rerunning CV.

Status now: **new 95-fit runner code and independent
audit have been added, but have not yet been run on
the user's Ubuntu server**.
