# TSPDT-inspired v1 — first matched 19 × 5 publication-data run

**Registered before any real-cohort TSPDT-inspired model fits.**
Status at writing: the separate synthetic-only DEV gate **passed**
on compute-server checkout SHA
`aa4dae241ba583394501f2f11efd757c8d74ad46`.
BigTSP optional smoke was **SKIPPED_NOT_INSTALLED**.
No TSPDT real-dataset results have yet been computed.

## Baseline and fixed parameters

The model label is
`TSPDT2011_INSPIRED_RECONSTRUCTION_V1` (**not** the
original 2011 source implementation). Independent implementation:
`R/baselines/tspdt_reconstruction_v1.R`, unchanged since synthetic
gate. Code semantics are specified in
`docs/TSPDT_RECONSTRUCTION_V1_METHOD.md`.

- Training-only MAD f100 instead of historical Relief-F f1000.
- Node-local strict-inequality TSP pair search on each node's
  **training** subset; no global dense BigTSP design.
- min_leaf=5, max_depth=3, Wilson error approximate pruning
  (cf=0.25), fixed tie convention and Laplace alpha=1 class
  probability reconstruction.
- No parameter optimization or use of outer-test labels
  in feature selection, splitting, pruning or leaf fitting.
- `BigTSP::tsp.tree` is a **different** pair-feature CART
  model; it is **not** part of this first run. No automatic
  package installation or replacement of the original method.

## Matched historical experiment

- Original 19 ARFFs: CompCancer11 + legacy8, exact SHA256 of
  publication `main=4079a5ac400f1d2cfded523c7000f655470925ce`.
- Exactly 5 canonical outer folds, seed42, on each original
  dataset: **95 TSPDT fits** total.
- Original RLR, RLMT-v1 and RLMT-DC models are **never
  recomputed** here. They were or are being calculated in
  a separate 285-fit run:
  `results/rlmt_dc_published_data_5fold_v1`.
- This baseline **must not launch** before the 285-fit
  runner exits successfully AND the four independent audits
  are finished:
  `published_5fold_audit_v1`,
  `published_5fold_accepted_replay_v1`,
  `published_5fold_firstfold_parity_v1`,
  `published_5fold_summary_v1`.
- Each TSPDT fit must match exact train/test sample IDs and
  class levels from the corresponding completed RLR fold
  **before any TSPDT fitting takes place**.
- Original RLR predictions and audited logloss are used
  **read-only**, from their existing output directories.
- Record each fit, serialized tree, node structure,
  held-out predictions, metrics and failures. Never discard
  inconvenient datasets/folds or restart only failed records.
- Output folder on independent TSPDT worktree:
  `results/tspdt_recon_19x5_v1/`. It must not exist
  before the strict launcher starts.

Frozen setup:

- `configs/tspdt/tspdt_reconstruction_19x5_v1.json`
- `scripts/tspdt_recon_19x5_launch.py`
- `R/baselines/50_tspdt_published_19x5_run.R`

Start with Python `--preflight-only` to validate **all 19**
source SHA256 values, original manifest, 285-fit completion,
source DEV gate, immutable algorithm file and R syntax.
The actual launcher uses a detached process, recording
`CODE_SHA.txt`, `REFERENCE_MODEL_SHA.txt` and
incremental `run_manifest.csv`.

## Independent confirmation of the 95 baselines

`scripts/tspdt_recon_19x5_audit.sh` runs
`R/baselines/51_tspdt_published_19x5_audit.R`
**without retraining any model**.

For all 95 folds, the auditor:

1. Rebuilds the exact train/test split from original ARFF.
2. Rechecks test sample IDs, labels and class ordering
   against original frozen RLR predictions.
3. Re-loads each `model.rds` and replays all held-out
   predictions and leaf routing.
4. Independently reconstructs logloss/Brier from Laplace
   scores, balanced accuracy, macro-F1 and accuracy from
   class decisions and published test labels.
5. Computes all corresponding reference RLR accuracy,
   balanced accuracy, macro-F1, plus logloss from its
   independently audited 285-fit prediction artifacts.
6. Aggregates **95 within-fold contrasts into 19 dataset
   means**, refusing partial output and preserving errors.
7. Requires `TSPDT_95FOLD_AUDIT=PASS` before treating
   the comparison as technically complete.

**Interpretation:** the selected 19 datasets are historically
known to the PI and previously used for RLR/RLMT development;
five folds share patients and are not independent draws.
The 2011 TSPDT paper used different procedures and feature
selection. All estimates are exploratory, not evidence of
prospective external validation or a reproduced original TSPDT.
Laplace leaf-score logloss is a **new reconstructed metric**
and should be labeled separately from standard probabilistic
RLR/RLMT outputs.

## Noninterference boundary

The dedicated branch
`feature/tspdt-baseline-evidence` is checked out as a
**detached separate worktree**. Do not run its launcher in
the original RLR worktree. Do not merge experimental
code into publication `main`, change the 285-fit
RLMT benchmark, install unverified CRAN archives, or copy
source from BigTSP into the original TSPDT reconstruction.
