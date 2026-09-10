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

After the repository is pushed to `https://github.com/mczajkowskipb/relative-logistic-regression`, run one smoke check from a fresh clone. No additional benchmark reruns are required for the numerical lock.
