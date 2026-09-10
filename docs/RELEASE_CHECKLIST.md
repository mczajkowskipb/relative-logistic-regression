# Release checklist

- [x] Canonical source/config capture: stage-6 missing files = 0.
- [x] Primary configs rebuilt from canonical `truth_rebuild` configs.
- [x] Machine-specific paths removed from publication-default code/configs.
- [x] PTB fixed protocol corrected to discovery 103/44 split; predictive outputs unchanged.
- [x] PTB ORA rerun training-only and outputs/GMT/checksums bundled.
- [x] R/system/package versions captured.
- [x] Paper-number manifest updated for PTB protocol and ORA.
- [x] MIT License selected and included.
- [x] legacy8 checksum/provenance limitation accepted and documented; no further source search required.
- [x] Final clean repository URL locked as `mczajkowskipb/relative-logistic-regression`.
- [x] Final main manuscript and Supplement compiled with complete figure assets.
- [x] Final PDFs visually inspected.
- [x] Quick release validation passes.
- [x] Clean public GitHub repository created and publication release pushed.
- [x] Fresh-clone quick verification passed (`python3 scripts/check_release.py --quick`; `scripts/smoke_test.sh` is the equivalent wrapper).
- [ ] Make `relative-logistic-regression-working` private after the clean public repository is verified.

No further scientific/benchmark computation is part of the release checklist.
