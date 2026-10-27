# Yeoh hard-cohort RLR tree mechanism — COMPLETE D/V audit, 2026-10-27

**Status: technical CI SUCCESS, promising conditional LOGLOSS only, NO F1 improvement or heldout external confirmation.**

- Actual [GitHub Actions run #40496639422](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40496639422) **SUCCESS**, fitted source SHA `fb220fe59e88359bf95d104c57636f948561618a`.
- Dataset `yeoh-2002-v1` from public Schliep-source reconstruction, original paper-named cohort with same original outer first 5-fold seed42. The original processed **ARFF checksum parity is NOT verified**. Only outer TRAIN T (198 patients) is used; D=149, V=49, original external E labels **not used**.
- Frozen prior mechanism and hyperparameters, no retuning: patientwise rank m10, original RLR K50/d2/alpha .5/lambda.min, at most six D-selected ordinal anchor root gates, up to eight branch ternary pairs, class-TSP or residual-selected pair candidates, conditional ridge-like elastic-net correction alpha=.5 lambda=.10. Pure rank-based within-patient comparisons.
- **Two independent R experiments** in one CI job, with archived `metrics_V.csv` and `predictions_V.csv` per mechanism:
  - class-selected vs OOF residual-selected local pair candidates, fitted against parent in-sample offset;
  - identical residual-selected pairs NOMINATED ONCE per root, local vs equal-feature global, each with in-sample vs OOF parent fitting offset.
- Independent standard-library Python replay downloaded actual GitHub archived artifact zip (ID `12730766758`), checked **54 root×method metric records** (24+30), **2646 V patient-method probability rows** (1176+1470). Recomputed logloss, Brier, positive-F1, Macro-F1, BA to absolute tolerance 1e-9, checked patient IDs, truth, n=49, all six roots, original RLR baseline parity, and real frozen inner V partition: PASS.

## Six proposed-root mean inner V scores

| Method | Mean V logloss ↓ | Mean V Macro-F1 |
|---|---:|---:|
| Original RLR50 | .01161792 | 1.000 |
| Local class-TSP-selected pairs, IN-SAMPLE parent offset | .01161792 | 1.000 |
| Local residual-selected pairs, IN-SAMPLE parent offset | .01161792 | 1.000 |
| Same nominated pairs globally, IN-SAMPLE offset | .01161792 | 1.000 |
| Same residual pairs globally, OOF offset | .00967637 | 1.000 |
| **Local conditional residual pairs, OOF parent offset** | **.00350969** | **1.000** |

- Tree OOF offset improves mean V LL by **.00810823** versus original RLR50 and **.00616668** versus the equal-feature global OOF control.
- The local OOF-offset tree has a lower V logloss than its corresponding OOF-offset global SAME-pair control for **6/6 D-proposed roots**, but all six use the *same 49 V patients*; they are NOT independent corroborations.
- Aggregate active child coefficients across six candidates: 0 with IN-SAMPLE offset, 47 with OOF offset. This counts multiple fitted models, **not** 47 distinct biological genes, and not one patient-specific rule card.
- All six methods already score Macro-F1=1.0 on this inner V: **NO classification improvement** is demonstrated. In the already previously inspected first **outer E fold**, original RLR50 positive-F1 was .75 versus kTSP50 .90 (different model-training sample); the current inner V result cannot fix this external evidence. No new E outcomes were accessed to select the method.

## Scientific interpretation

The new matched feature-controlled D/V diagnostic shows that (i) in-sample RLR offsets can hide residual conditional information, (ii) using crossfitted offsets can activate it, and (iii) **in Yeoh's inner V, conditional routing adds real logloss benefit beyond global application of exactly the same pairs**. This is the clearest training-only evidence for actual model-tree structure in the project, but it has no F1 ceiling left to improve and the cohort was previously explored.

A ***scientific next gate*** is to freeze the exact new offset/anchor/candidate specification and check **all 11 original paper named cohorts on D/V** with predefined inner acceptance against original RLR50 AND samepairs global. The old original RLR manuscript should remain frozen. Only then consider the expensive full 11×470 outer E reproducibility run, ideally original processed ARFF SHA checked, with real switchBox kTSP and value LR. No result from these repeated development folds supports a publication claim yet.
