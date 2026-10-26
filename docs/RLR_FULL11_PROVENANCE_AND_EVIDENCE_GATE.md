# RLR publication-only same-pair component evidence: 11-cohort gate

**Date:** 2026-10-26. The original RLR manuscript and `R/new/07_world_protocol_utils.R` on `main` remain unchanged.

## Experimental question

Do separately learned individual RLR coefficients yield improved *held-out logloss* relative to one equally regularized coefficient on the arithmetic mean of **exactly the same signed ternary m=10 gene-pair features**?

Fixed pair nomination uses the exact original `compute_top_pairs_switchbox`, top-1000 shortlist, strict gene disjointness d=2, K=10/20/50. The original RLR is fitted using **the original literal function** `fit_rlr_train_test()`. The control uses `glmnet` with the original alpha=.5, 5-fold inner CV, `lambda.min`. A constant all-zero second column is excluded and verified zero to allow the univariate penalized vote in `glmnet`. No new data or decision trees are introduced.

**Primary target:** 11 original binary CompCancer datasets, 470 original outer folds per method (nine 5-fold × 10 repeats; West and Yeoh 5-fold × 2 repeats), seed 42, 3 K values, 2 methods per K; **2820 fits** if completed. Dataset, not fold, is the unit for exploratory paired Wilcoxon inference; Holm correction across the three K configurations is used only for this new exploratory comparison.

## Provenance STOP / outstanding limitation

The 11 original preprocessed ARFF files are **not distributed in this public GitHub repository**. The original locked SHA256 values are in `checksums/DATA_SHA256_EXPECTED.txt`. The publicly sourced Schliep filtered text matrices were downloaded and independently checked by [source gate #40325493713](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40325493713).

All 11 original dataset names had two-class labels and verified sample counts. The catalog and downloaded probe count agree for 10/11. **Golub-1999-v1 source has 1868 filtered probes whereas the source catalog lists 1877**. The source first header row for Golub and West is a repeating subtype annotation, so the pipeline uses stable ordinal sample IDs while preserving exact column order. This DOES NOT show original processed ARFF equivalence.

- `scripts/rlr_original_arff_sha_gate_v1.py`: fail-closed checker for 11 local **original** ARFF files. Only SHA-matched originals are eligible for a claim of exact full paper reproducibility.
- `scripts/rlr_full11_source_gate_v1.py`, `scripts/rlr_full11_fetch_one_v1.py`: full source inventory plus single-DS public reconstruction. Exact ARFF parity remains explicitly `NOT_VERIFIED`.
- `R/experimental/rlr_full11_same_pairs_v1.R`: literal parent RLR + equally regularized vote in same CV folds, with 1-fold smoke mode and complete 470-fold mode.
- `scripts/rlr_full11_same_pairs_audit_v1.py`: patient-level independent replay, complete fivefold coverage per repetition.
- `scripts/rlr_full11_aggregate_v1.py`: aggregate by dataset and enumerate all 2^11 signed-rank signs; descriptive Holm adjustment.
- `.github/workflows/rlr-full11-same-pairs-smoke-v1.yml`: [11-cohort CI smoke #40326452733](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40326452733) on source gate branch, launched.
- `.github/workflows/rlr-full11-470-cv-v1.yml`: gated 470-fold matrix + automatic independent summary on this full branch. **The full workflow MUST NOT be claimed completed until a source-controlled launch trigger and audit PASS exist.**

The reconstructed-source 470-fold run, if successfully executed, can only provide an **additional exploratory diagnostic**; the original manuscript results cannot be retroactively labelled reproduced unless the original ARFF hashes match. If the reconstructed results disagree materially, investigate source/probe/preprocessing and positive-class orientation **before modifying the science**.

## Stop rule

Avoid further search over alternate thresholds, tree structures, rank encodings, or hyperparameters on the already examined data. If the fully audited equal-regularization control proves competitive, refine the novelty claim for the *integrated* RLR pipeline instead of claiming that unequal weights alone explain the advantage. Preserve negative findings. No automatic manuscript update.
