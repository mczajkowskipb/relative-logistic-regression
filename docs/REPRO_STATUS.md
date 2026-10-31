# Reproducibility status

## Canonical state

The publication code in this release was collected from the RLR working tree on 2026-09-10 at git commit `b3c4d1316b0ab4a5e481dd50e3dd9e0c3a07db07`. The stage-6 collector reported `missing_files=0`.

Canonical reference counts at collection time:
- CompCancer: 4,230 `metrics_test.json` and 4,230 `predictions_test.csv`;
- PTB fixed: 9 + 9;
- PTB LOCO: 36 + 36;
- legacy8: 3,600 + 3,600.

The primary source/config files in this repository are stage-6 captures with only portability edits: machine-specific paths were replaced by repository-relative paths or environment/CLI parameters, and scripts that source `07_world_protocol_utils.R` now resolve it relative to the script location. No model logic, metric definition, statistical family, or benchmark schedule was intentionally changed. Original and publication-source checksums are provided under `checksums/`.

## PTB fixed correction

The canonical PTB fixed predictive outputs are unchanged. They correspond to the 103/44 discovery-data train/test split. GSE83456 is the 106-sample external-validation cohort in the source benchmark, but it is not the test set used for the reported canonical PTB fixed metrics. The corrected ORA uses only the 103-sample training data. See `PTB_FIXED_PROTOCOL_CORRECTION.md`.

## Environment captured

R 4.6.0 on Ubuntu 22.04.4 LTS. Package versions are listed in `environment/package_versions.tsv`; `sessionInfo.txt` and system details are included. A full dependency-closure `renv.lock` was not reconstructed from the historical run and is therefore not fabricated here.

## Reproducibility boundary

The historical legacy8 upstream provenance cannot be reconstructed more precisely from the archived materials. The exact processed inputs used for the locked state are therefore identified by filename and SHA256 rather than by a newly asserted upstream source. This limitation is explicit and accepted for the historical/descriptive benchmark.

The public repository does not mirror every split-level frozen output. Those files remain in the authors' private frozen archive. Publication-facing provenance, canonical configurations, checksums and the complete corrected PTB ORA are included publicly.

## Final external verification

Fresh-clone quick verification of the public repository passed on 2026-09-10 (`python3 scripts/check_release.py --quick`; `scripts/smoke_test.sh` is the equivalent wrapper). No additional benchmark reruns are required for the numerical lock.

## Original ARFF verification (owner-performed, 2026-10-30)

The authors ran SHA256 against the frozen `checksums/DATA_SHA256_EXPECTED.txt` manifest on the original local processed data directory: **45/45 matching files**. The full eleven original CompCancer ARFFs and their 9x10/2x2 execution-shard counterparts also passed a source/shard preflight (**11/11**). The exact byte identity preserves original sample, probe and class-factor ordering; it does **not** license redistribution of the ARFF matrices.

The repository now includes the public, read-only `scripts/verify_original_arff_preflight.py` to independently redo the checksum verification on any authorized local copy. This is a provenance confirmation for the original **processed input files**. Do not retroactively label historical public-Schliep-reconstructed exploratory GitHub Actions fits as having run against those exact original local ARFFs; the original-paper CompCancer algorithm and 470-fold numerical lock have not been changed.

## Supplement S7.1 separate matched encoding ablation

The 11 original ARFF datasets (one five-fold outer repeat = 55 held-out folds) were processed on the author's server with original RLR20 and two controlled binary-pair penalized logistic refits (ENET and MCP), totaling 165 method-fold fits. The output was independently replayed from per-patient predictions, and a second SHA256 check across original and sharded ARFF was **PASS**. This experiment is separate from the primary full-470 protocol and is explicitly **exploratory**. The binary models are methodological controls, **not** complete published PurIST or DeCAF clinical models. The author-run audit verified dataset-equal RLR logloss 0.237943894 versus binary ENET 0.290442282 and binary MCP 0.291233444. For the reproducible code and detailed conditions, see `docs/COMPCANCER_TERNARY_BINARY_ABLATION_REPRO.md`. Warnings of MCP model saturation were seen in some fits; no claim of superiority over original DeCAF/PurIST is warranted.

## Manuscript handling

The old `manuscript_lock/RLR3_submission_final.*` PDFs/TeX were an **obsolete manuscript snapshot** and have been removed from the current public branch together with their stale PDF hash list. Current submission PDFs and source are exchanged via the journal's editorial system, not published in this source-code release. Historical repository commits still preserve provenance of the obsolete snapshot.
