# RLMT-v4 — frozen 19×5 exposed development replication (2026-10-22)

Source branch descends unmodified from
\`feature/rlmt-v4-rank-repair\` at
\`e4276d337ef4674f52bbacc398e2dea14e7c84ba\`.
No changes to algorithm, feature nomination rules, penalties, or seed.

## Prior 2×3 pilot — real finished and audited

Original RLR patient-level parity **PASS**, maximal absolute
prediction deviation 4.44e-16, all six source fold IDs exact;
synthetic 5 model gate PASS; real 30/30 outcomes audited PASS.
On exposed \`chen-2002\` and \`gordon-2002\`,
three paper folds per cohort, equal weight each fold:

| Method | Paired ΔLL vs original RLR | Fold W/L/T vs original | Real gene-pair routes |
|---|---:|---:|---:|
| ORIGINAL_RLR | +0.000000000 | 0/0/6 | 0/6 |
| GLOBAL_CAL | −0.000489915 | 5/1/0 | 0/6 |
| CONF_ROUTE | −0.002551260 | 5/1/0 | 6/6 |
| ACTIVITY_ROUTE | −0.000816796 | 5/1/0 | 6/6 |
| SOFT_GAP_ROUTE | −0.009452910 | 5/1/0 | 6/6 |

SOFT_GAP vs global-calibrated control: −0.008962995
mean heldout LL (paired). This is an INTERESTING DEVELOPMENT
SIGNAL, not proof of method superiority, because method ideas
and datasets were repeatedly exposed. Source: operator final
auditor stdout on 2026-10-22. No final 19×5 results yet.

## New fixed evaluation, without reselecting winning folds

All 19 previously analyzed publication datasets, the exact
original seed-42 5 folds = 95 matched outer-training bundles
per method, for five modes:
ORIGINAL_RLR, GLOBAL_CAL, CONF_ROUTE,
ACTIVITY_ROUTE, SOFT_GAP_ROUTE.

- 475 method-fold prediction files. Original archived 285-model
  references are READ ONLY and verified against exact same
  source-patient row IDs and training/test partitions.
- Parent is \`fit_rlr_train_test\` literal original
  switchBox RLR (K10/m10/d2/pool1000/alpha.5/lambda.min)
  rather than the previously different aligned-v2 CV root.
- For every outer fold, strict original-RLR 3-fold OOF
  train-only predictions are used for supervised pair nomination
  and calibration; 4096 random candidates plus the 10 original
  RLR pairs. Exact parameters from v4 pilot:
  λβ=.02, λγ=.10, min leaf support=8, rank-gap scale=10.
- Hard-stop on any parent prediction difference >1e-9 from
  archived original, source SHA change, fold/sample/label
  mismatch, prior run output, or a failed independent audit.
- Independent Python auditor checks heldout per-patient
  probabilities, rescored LL and Brier, plus balanced accuracy
  and macro-F1 at fixed 0.5. Reports equal-dataset means,
  comparison both against ORIGINAL_RLR and GLOBAL_CAL,
  dataset-level W/L/T, 19-fold LODO diagnostics, gene-route
  counts and fitted gamma.

## Scientific caveats

The original observed six-fold pilot was **EXPOSED**,
as are all nineteen cohorts. Even a favorable broader
retrospective result is hypothesis-generating. Choosing
SOFT_GAP because of its 2x3 results is selection;
the broader 19x5 retains ALL four controls/alternatives
to avoid reporting only a lucky screen.
Data-derived candidate choice of 4096 pair
routes and coefficient estimation on the same
outer-training sample creates a significant
selection-overfitting risk *despite honest original-RLR OOF*.
Do not claim it is eliminated by cross-fitting.

A differentiable gene rank-gap correction is a
**one-rule soft relational tree / logistic calibration model**,
not a conventional hard-split decision tree.
Name and comparator claims must reflect this.

If the broader descriptive analysis shows meaningful
stable gains beyond parent/calibration, the next
publication-grade step is a pre-declared exact
outer-training nested procedure and a NEW
untouched patient cohort/lab batch/data source
for independent replication. Original finished RLR
paper must remain unchanged until separately
validated novel results exist.
