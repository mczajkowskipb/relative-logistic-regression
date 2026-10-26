# Real original CUDA TSPDT baseline — readiness and parity audit

**Date:** 2026-10-26. **User supplied original source**: `cudatsp4(1).zip`, SHA256 `d788812836097aa8bae342381c2aa7727e0c53279b35e1b1447b12c897ac0232`, 9 files. Original kernel SHA256 `2b47a5834f5fe1d0dc3bab40b4c4578f2b60567ba1c40b97c21db119c802a770`. Original source is private researcher material. **Do not upload, release, commit, paste or otherwise distribute that CUDA code in this public RLR repository.** No CUDA source was committed in this branch.

## Why this gate takes priority over another new tree hypothesis

All previously evaluated TSPDT branches were **new simplified greedy proxy trees**, not the author's actual CUDA TSPDT implementation, so no valid real original-TSPDT superiority claim was made. The original RLR paper is **not yet published**. The user wants a single genuinely new tree-relevant extension of RLR for the habilitation, on **only its 11 original binary CompCancer datasets**, with matched original controls and preserved prior manuscipt. Before proceeding, actual CUDA TSPDT baseline reproducibility must be established on the same patient-train/E split as RLR.

## What was actually verified

The original source archive was inspected without publishing it. The actual original `kernel.cu` has a `main(argc,argv)` accepting a JSON config file and optional GPU index, and reads external TRAIN and TEST files through original `load_files.h::readCSV`. Its dense ARFF-like parser infers integer class labels from **order of training class appearances**; this must be saved and remapped. The original `foldsCV=-1` means external train/test, `normalization=0` prevents independent train/test normalization; for one omics block, `multiomicsthresholds` must be EXPLICITLY set to gene count because original CUDA allocation uses the vector's size before model selection fills the default.

A standalone **private** adapter and source-patch toolkit was generated in the active research session, consisting of `prepare_original_tspdt_benchmark.py`, `patch_private_original_cuda_tspdt.py`, `collect_original_tspdt_predictions.py`, and isolated Python unit tests. They are delivered as a private ChatGPT conversation ZIP named `RLR_original_CUDA_TSPDT_parity_adapter_v1.zip` (ZIP SHA256 `127d81143b7a4e2ff831761b38a4f4ed66b21b781325132e0456ba1cdb274b89`). The toolkit **does not contain original or patched CUDA source**, and is not claimed to be persisted within the GitHub repo.

Verified in local controlled tests:
- **9/9 independent Python tests PASS** on fictional binary tabular inputs and simulated CUDA stdout, including strict original ARFF checksum enforcement, patient ID alignment, precise split preservation, class mapping, class truth mismatch refusal, reconstruction labeling, private-source patch anchor, and refusal of missing/duplicated GPU predictions.
- The archive SHA and original `kernel.cu` hash were verified; `patch_private_original_cuda_tspdt.py` generated a separate, local, SHA-locked instrumented copy of the supplied kernel and inserted one per-patient `TSPDT_E_PRED,<row_index>,<true_class_idx>,<predicted_class_idx>` line into the **existing test confusion matrix loop**, leaving its prediction and tree-induction logic unchanged.
- **Actual user-supplied original `load_files.h` compiled with `g++ -std=c++17` and successfully read generated TRAIN+TEST ARFF**, reporting `REAL_ORIGINAL_CPP_READCSV_LOADER=PASS`, 4 synthetic training patients, 2 heldout patients, two recognized class IDs and correct numeric gene counts. This verifies the real original C++ data loader, not the original CUDA classifier on GPU.
- `collect_original_tspdt_predictions.py` strictly checks E-row order, patient identity, actual source truth class and maps the CUDA class integer back to the original binary label. Reports BA/Macro-F1 on **hard predictions only**, not fabricated probabilities or calibrated logloss.

## What is still blocked / NOT a result

1. No verified original **processed ARFF bytes** present in the public repository: they are ignored by `.gitignore`. The exact 11 original source SHA256 checks are recorded in `checksums/DATA_SHA256_EXPECTED.txt` but must be provided locally. The public Schliep reconstruction is NOT exact original ARFF parity; exploratory mode is separately and unmistakably labeled.
2. The original RLR locked patient-row ordering mapping and exact original heldout fold membership are not yet obtained for all original processed 11 ARFF datasets. The adapter **refuses** original-ARFF comparability without such mapping.
3. Original CUDA source has NOT been fully built, instrumented and run on any GPU during this work session: the model's container lacks `nvcc` and NVIDIA GPU, and the GitHub Actions jobs used for RLR are CPU Ubuntu runners. Do not invent original CUDA performance, passing compilation or patient-level predictions.
4. Original CUDA baseline parameters/configuration from the author's TSPDT publication (including `algorithmType`, `splitSize`, `rankTest`, `minNodeSize`, `topSelection`, all normalization/threshold settings) still require scientific parity inspection. The bridge's default parameters are provisional, not a claimed exact replication of that paper.
5. Original CUDA program did NOT previously emit per-patient E-predictions; original aggregate `CV FINAL` accuracy/macro-F1 is not directly comparable to foldwise balanced accuracy. The local private patch resolves extraction *if* actual GPU execution passes.

## Gate for a defensible habilitation-related RLR+TSPDT conclusion

Run the **original** CUDA TSPDT privately on at least the EXACT SAME heldout patient-splits and original processed 11 paper datasets as frozen RLR. Compare matched E-patient BA/Macro-F1, tree size and pairs on route. Use dataset-level units for inference. Preserve dataset source SHA, original config SHA, original CUDA source SHA, split patient IDs, true class mapping and real predictions. Never substitute the greedy proxy, a reformulated method or a synthetic test for actual original CUDA execution.

**Scientific decision:** Do not introduce more post-hoc tree structures to the unpublished RLR manuscript until the original actual CUDA TSPDT baseline and original data provenance gate can be evaluated. There is no evidence yet of a new tree beating original TSPDT. Existing certified RLR early-exit results are meaningful computationally but not yet established as a novel tree learner.
