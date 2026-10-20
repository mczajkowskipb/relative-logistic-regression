# RLMT GATE 4 — metadata screening rules (pre-outcome)

This protocol was authored *before reading any RLMT performance on the six
remaining CompCancer matrices*. It extends the scientific checkpoint in
`docs/RLMT_GATE3_SCIENTIFIC_CHECKPOINT.md`.

## Source and scope

The original RLR publication already used the CompCancer data in the
`data/CompCancer_truth_9x10/` shard; the canonical expected processed-file
SHA256 values are in `checksums/DATA_SHA256_EXPECTED.txt` under
`data/CompCancer/<dataset>.arff`. The nine entries are the filtered matrices
documented in `docs/DATA_ACQUISITION.md`.

Exclude `alizadeh-2000-v1`, `chen-2002` and `singh-2002` due to their
documented use in RLMT development. The remaining pre-specified inspection
set consists of `armstrong-2002-v1`, `bittner-2000`, `golub-1999-v1`,
`gordon-2002`, `pomeroy-2002-v1`, `shipp-2002-v1`.

This is a **historical-benchmark, not an independently recruited external
cohort**. Rows are typically assigned synthetic sample IDs in the existing
RLR loader. No patient-level identifiers or cross-source repeat-sample
relations can be inferred from numerical ARFF matrices alone.

## Fixed, performance-blind inclusion screen

Include **all** and **only** unseen datasets that satisfy every criterion
below; do not inspect predictive results when evaluating eligibility.

1. An exact corresponding file exists, with SHA256 matching the canonical
   entry `data/CompCancer/<dataset>.arff` in the repository manifest.
2. The loader recognizes exactly one binary class column; no missing class
   labels, no ambiguity in class levels.
3. At least **80 samples**, **25 instances of each class**, and **20 numeric
   predictors**. These conservative feasibility thresholds are metadata-only
   safeguards for small-leaf CV, *not tuned hyperparameters*.
4. No non-finite values in the numeric predictors used by RLR and no
   duplicated complete numeric feature rows. For any failure, flag it;
   do not silently deduplicate, impute, or choose a different file.
5. Record whether sample IDs are real or synthetic. Neither synthetic names
   nor row-level duplicate absence proves patient independence. An
   eligible dataset can remain only an exploratory benchmark if patient
   dependencies cannot be ruled out.

The audit is descriptive and read-only for source data. It writes only
`results/rlmt_gate4_metadata_screen_v1/`. Failure of manifest integrity
or parse checks blocks any RLMT training in GATE 4. If no dataset passes
eligibility, stop and find a genuinely new cohort instead of weakening
thresholds after observing class counts.

## Planned later modeling; NOT authorized by this inventory

Once metadata screening is complete and documented, freeze the exact code
version and a new output root before starting fits. If eligible datasets
exist, apply one stratified 5-fold run (`seed=42`), unchanged
K10/m10/d2/pool1000/alpha0.5/lambda.min, B20, 1-SE,
`min_child_n=16`, `min_child_class=6`, inner folds 3, and
class-stratified leaf-internal CV. Include all screen-eligible datasets,
with no performance-dependent filtering. Report logloss as primary,
Brier/AUC/F1 as secondary, realized leaves, failures/warnings, and
between-dataset summaries as such. Do not treat repeated folds within
the same source dataset as independent cohorts.

Because the current RLMT pair pool and screening share training labels
with inner validation decisions, their inner-OOF gains are *model selection
scores* and not unbiased estimators of generalization benefit.
Additionally, if this protocol is used after the held-out Singh fold was
examined and method choices were informed by it, it is exploratory
method-development evidence. Independent scientific confirmation needs
unseen cohorts and a truly frozen selection algorithm.
