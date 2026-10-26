# Original CUDA TSPDT parity bridge V3 — pinned reconstructed RLR source, 470 splits, automatic class mapping

**26 October 2026. Original unpublished RLR paper and core `main` code remain unchanged. Original researcher's `cudatsp4(1).zip` CUDA source remains PRIVATE — not committed to GitHub.**

## Goal and scope

Before claiming that any tree extension beats TSPDT, execute the ACTUAL original CUDA code, not the previous simplified greedy TSPDT proxy, on matched heldout samples from the same **11 original binary CompCancer cohorts** as the unpublished RLR paper. The archived original RLR source-reconstruction protocol has **470 stratified train/E folds** and **8564 literal original RLR K50 patient predictions**; this is *not* an exact original processed ARFF reproduction.

The v3 adapter is delivered as a **code-only private conversation attachment** `RLR_original_CUDA_TSPDT_parity_adapter_v3.zip`, size 41138 bytes, SHA256 `1c95ae20943109c811e90f537b7333573d012c3bbaee51b0266310d4769b9246`. Contains 20 code/docs/test files, NO original CUDA source, no patient matrices or private CUDA executables. The earlier 470 frozen score/split bundle remains `RLR_11_original_paper_cohort_CV_splits_470_reconstructed_source.zip`.

## What was actually implemented and tested this cycle

1. **Double-SHA source acquisition** `fetch_pinned_compcancer_source_v3.py`: downloads only the exact public Schliep text of one or all eleven manuscript cohorts; verifies the text SHA256 then reconstructs the positional patient IDs, within-source column order, gene/probe names and conversion format used by the previously completed literal original RLR CI run. The final converted CSV SHA256 must match the archived `input_source_manifest.json` byte-for-byte. Rejects different public source versions, probe counts or preprocessing. No source downloaded in this CPU environment: the network acquisition mechanism was tested on synthetic pinned data only and must pass its real source hashes on GPU host.
2. **Patient and positive-label preflight** `preflight_reconstructed_rlr470.py`: checks the actual source CSV SHA against RLR's recorded converted matrix, every K50 RLR patient prediction and its y truth vs source class, and exact disjoint fivefold test coverage for each repeat. Infers original RLR positive/negative biological class strings independently of CUDA output. Rejects mismatches before GPU training. Emits `reconstructed_rlr470_preflight.json` and prepares only provenance-verified train/test/config split files.
3. **Stricter actual CUDA evidence**: the real-GPU runner records input source SHA256; the independent `aggregate_paired_cuda_rlr.py` compares against the exact previously used RLR input-matrix SHA256 and the verified class mapping. It also checks actual CUDA per-fold predictions file hashes against GPU evidence index before scoring. Full-schedule claims now REQUIRE the positive-class/source binding certificate.
4. **Two real ORIGINAL private source C++ CPU integration checks**: actual original `structures.h` C++ config parser compiled and read v3-produced compact JSON paths and omics-threshold correctly — `ORIGINAL_CPP_JSON_TWO_FILES_AND_OMIC_THRESHOLD=PASS`. Separately actual original `load_files.h::readCSV` compiled and read v3-produced original-style train/test ARFF files with correct patient sample counts, two class coding/order, finite features and exact dimensions — `ORIGINAL_CPP_CUDA_LOADER_TRAIN_TEST_CLASS_ORDER_AND_DIMENSIONS=PASS`.
5. **29/29 Python tests PASS**, including negative cases: byte-changed source data, wrong patient, flipped RLR binary truth, wrong manually supplied positive label, RLR/CUDA source SHA mismatch, CUDA E score tampering, falsely claiming single fold is complete paper schedule, overwriting science evidence, and numerical artifact integrity.
6. **One-command PRIVATE GPU smoke/full runner** `run_private_one_original_paper_cohort_v3.sh`: on a Linux/WSL2 NVIDIA CUDA machine, validates nvcc/nvidia-smi, re-downloads pinned public matrix, checks 470-split evidence, prepares either one original RLR heldout fold or all 50/10 folds for one dataset, runs the real original CUDA on a SHA-locked private copy instrumented ONLY to print E-patient predictions, and compares matched hard BA/macro-F1 to literal original RLR K50. Fail-closed on CPU-only environment: `STOP: actual NVIDIA nvcc and nvidia-smi are required. No GPU run was performed.`

The original researcher CUDA code, gene expression data, full CUDA executable and runtime logs **must stay in private local directories**, never in this public repository. Some original source CUDA settings (`algorithmType=4`, `splitSize=3`, `rankTest=1`, `topSelection=12`, `minNodeSize=7`) remain PROVISIONAL, not independently matched to original published TSPDT experiments.

## Scientific status — no invented actual baseline

| Gate | Status |
|---|---|
| Eleven original cohort identities, 470 RLR heldout folds, 8564 archived RLR K50 heldout probabilities | PASS (public reconstructed source) |
| Reconstructed source pinned text+converted CSV SHA check implementation | PASS in synthetic reproducibility tests; real download/pinned SHA checks still require target host |
| Patient identity and class orientation verification algorithm | PASS in synthetic tests; real-source verification still requires input CSV |
| Python 29 isolated tests | PASS |
| Original private C++ settings + data loader on adapter-produced synthetic inputs | PASS |
| Original CUDA compiler/runtime on NVIDIA GPU | **NOT RUN** — no nvcc/NVIDIA GPU in current assistant runtime |
| Original processed paper ARFF checksum parity | **NOT VERIFIED** |
| Original CUDA TSPDT E-patient predictions, real BA/macro-F1 or superiority claim | **NOT AVAILABLE** |

**Recommended operational next step:** run v3 `smoke` on the researcher's RTX-equipped private Linux or WSL2 workstation, inspect true CUDA compile logs / signed patient evidence, then launch original 50/10 folds cohort by cohort. For original-paper ARFF numerical parity, separately obtain exact ARFF bytes, the original patient index map and refit RLR on those bytes; do not compare old reconstructed-source RLR probabilities against different original-ARFF CUDA inputs. No manuscript edits before a scientifically successful tree extension is supported.
