# RLMT-DC — 19-dataset × 5-fold exploratory benchmark v1

**2026-10-21 pre-run frozen protocol.** Prepared *after* fully audited
first-fold pilot. This is further analysis on the **already exposed**
published RLR / RLMT development corpus, not prospective independent
validation or inferential proof of method superiority.

## Locked scope

- Exactly 19 original RLR-publication processed ARFFs:
  11 CompCancer and eight legacy8, in the exact registered order
  of `configs/rlmt_dc/published_data_5fold_v1.json`.
- Match **every source-file SHA256** to the frozen publication
  `main` checksum manifest. Do not substitute transformed/new data.
- Repeat the *original* `make_repeated_stratified_folds` protocol
  with seed=42, n_folds=5, n_repeats=1, for each dataset. Use
  rep00_fold00, rep00_fold01, rep00_fold02, rep00_fold03,
  rep00_fold04 **for every one** of the 19 datasets.
  First folds are **recomputed and rerun** in a separate result
  tree for a uniform complete five-fold record; the audited first-
  fold pilot is retained unchanged and checked for parity separately.
- Exactly **95 paired outer splits × 3 methods = 285 model attempts**.
  Each method uses the **identical** train/test assignment and the
  shared outer-training top1000 supervised predictor pool.
- Fixed comparison methods:
  - RLR K10/m10/d2/pool1000, alpha=0.5, lambda.min, seed42,
    unchanged original implementation.
  - RLMT-v1 max_leaves=2, unchanged original implementation.
  - RLMT-DC-v0 max_leaves=2, discovery/confirmation and 4096
    label-agnostic proposals, unchanged DEV-v3-passed implementation.
- All class support (minimum 16 per child, six of each class)
  and confirmation 1-SE acceptance thresholds remain fixed.
  No posthoc revision after seeing Gordon, Singh, breast2, leukemia.
- All small datasets remain in scope. Report failed fits explicitly;
  do not silently remove, alter fold generation or improve
  imbalances by discarding difficult splits.

## Technical gate

Runner is an opt-in NEW `R/new/47_rlmt_dc_published_5fold.R`
with only protocol/output-root/count adaptations of original
`R/new/44_rlmt_dc_published_pilot.R`. The newly registered
launcher is `scripts/rlmt_dc_published_5fold_launch.py`.
An independent read-only R auditor
`R/new/48_rlmt_dc_published_5fold_audit.R` reconstructs
each test fold's labels/probabilities/metrics from its exact ARFF.
A separate accepted-tree replay `R/new/49_rlmt_dc_published_5fold_replay.R`
must verify all accepted serialized two-leaf models, not only
favorable datasets.

Runner and audits may write only under
`results/rlmt_dc_published_data_5fold_v1/`.
Fail closed before starting if the input files, paper SHA, Git HEAD,
configuration or source parser are inconsistent. All 285 attempts
and errors are logged in an incremental manifest. The run must never
rewrite historical publication/GATE/pilot data.

Freeze the executable code SHA *before* launching. A technical
PASS requires 285 triplet records, independent per-fold metric checks,
no unexpected sample/fold mismatches, exact one-leaf RLR probability
parity, and complete two-leaf replay for accepted models.

## Scientific analysis order (fixed BEFORE new fold results)

1. For every dataset compute **mean across the five paired
   outer-test fold logloss differences** `DC−RLR`
   (negative = better); do the same for RLMT-v1 and Brier.
   Fold-level differences are correlated and **not 95
   independent replications**.
2. Summarize all 19 dataset-level contrasts with mean, median,
   wins/losses/ties, accepted trees per dataset, number of
   model failures, and worst observed paired test loss / tails.
   Show full fold details without cherry-picking and include
   the known first-fold breast2 degradation.
3. Use original labels and sample IDs. Patient group independence
   cannot be assumed from processed ARFF without subject metadata.
   Legacy8 upstream raw-data provenance is incomplete.
4. Any cross-dataset comparison is exploratory; do not report
   a discovery-confirming p-value or superiority claim from
   these previously examined cohorts.
5. If improvement is inconsistent, stop/interpret it, rather than
   changing the support threshold, candidate count, seed, or
   confirmation rule retrospectively.

This **5×1** protocol is not the entire publication repeated-CV
schedule: original CompCancer had 10×5 for nine datasets and
2×5 for two; legacy8 10×5. PTB fixed/LOCO is separate and
not represented here. A complete published-protocol replication
is a later opt-in stage; do not mislabel this five-fold analysis.

## Existing pilot reference

First-fold pilot execution code SHA:
`814fee71fe6306db1e62b2ca47b1038df82d77a8`.

Two independent auditors passed on that prior result tree:
`PUBLISHED_DATA_AUDIT=PASS`,
`PUBLISHED_ACCEPTED_REPLAY=PASS`.
Three accepted DC splits: Gordon (gain), breast2 (harm),
leukemia (gain), with 16 unsplit ties. Original RLMT-v1 accepted
Singh and harmed heldout logloss. These observations motivate
a complete all-fold *description*, not new parameter choices.

GSE72094 remains on hold as an untrained source-intake candidate;
no GEO download or new mutation-target modeling is involved.


## Reproducibility enforcement

The complete auditing workflow is
`scripts/rlmt_dc_published_5fold_audit.sh`:

1. `R/new/48_rlmt_dc_published_5fold_audit.R`: read-only
   reconstruction of **all 285** paired test predictions/metrics,
   test sample identities and unsplit RLR exact parity.
2. `R/new/49_rlmt_dc_published_5fold_replay.R`: replay every
   actual accepted L2 split, routing, saved `cvfit.rds`
   and logloss from original ARFF feature matrices.
3. `scripts/rlmt_dc_published_5fold_pilot_parity.py`:
   compare **all 57 first-fold model outputs** against the audited
   prior smoke execution SHA
   `814fee71fe6306db1e62b2ca47b1038df82d77a8`,
   requiring matching test IDs, classes, split decisions and
   probabilities to 1e-12. If this fails, do not interpret scores.
4. `scripts/rlmt_dc_published_5fold_summary.py`:
   after all three checks pass, aggregate **paired
   logloss/Brier differences at dataset level** and report
   all nineteen groups without inferential p-values.

The orchestrator verifies the full run's code SHA from
`CODE_SHA.txt` and compares exact model source bytes with the
executed git commit; no post-run model changes are permitted.
Output directories are separate and refuse overwriting.
Each auditor is new and all are subject to live R/Python parse
verification by the actual Ubuntu server.

**Status at protocol authoring: scripts created, experiment NOT
yet executed.** No full five-fold performance estimates are claimed.
