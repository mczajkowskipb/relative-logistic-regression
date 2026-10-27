# Frozen cross-fitted offset mechanism comparison — RLR residual-local model trees V1

**Study stage:** 2026-10-27 — D/V-only development experiment, **NOT independent outer E benchmark**. Original RLR K50 source and locked unpublished manuscript unchanged. No change to train-only gene-universe ranks, candidate features or study labeling.

## Scientific issue

Earlier local conditional RLR candidates used training-only original RLR50 **out-of-fold residuals** for selection but fit local logistic coefficients using **full-in-sample fitted RLR50 probabilities** as offsets. This inconsistency can make local coefficients zero if the original RLR learns the D labels almost perfectly, even where OOF errors contain relevant conditional signal.

The ablation changes only the parent offset used to fit the local coefficients:

- **INSAMPLE_OFFSET:** `eta0_D=logit(p_RLR50_fullD(D))`
- **CROSSFIT_OFFSET:** `eta0_D=logit(p_RLR50_3foldOOF(D))`

At prediction time, **both use the identical fitted original RLR50 probability on heldout inner V** as the base `eta0_V`. The same D-selected **root gate, anchor gene IDs, threshold, and branch-specific ternary m10 rank-pair IDs** are frozen for both methods.

## Equal-feature controls and no-tuning requirement

- Original RLR50 probability for V.
- Local residual-gene expert at lambda .10 with in-sample original RLR offset in D.
- Same local expert with *crossfitted* RLR offset in D.
- Same globally fitted non-tree correction on **all identical nominated pair IDs** at lambda .10 with in-sample original RLR offset.
- Same globally fitted non-tree correction on **all identical nominated pair IDs** at lambda .10 with crossfitted original RLR offset.

Do not filter global model pair union by which features happen to be active under one fitted offset: that's an unfair comparator and can depend on regularization. Both controls are allocated the same nominated gene-pair union, with only offset training method changed. When local penalties shrink all coefficients to zero, that is a meaningful negative result rather than a changed feature set.

The outer original E fold is never used; candidate roots are nominated exclusively on D, and all quoted model comparisons use the same heldout inner V patients. Root shortlist max 6, candidate roots max 240, local nominated pairs max 8 disjoint pairs per branch, fixed alpha .5 lambda .10 and original global parent RLR50. The three original paper-named reconstructed cohort cases are predeclared: **Chen, Gordon, Singh**. Use only the original first split seed42. Nested crossfit uses original RLR50 with supervised features selected afresh per OOF subfold.

## Output and interpretation

Files per case: `metrics_V.csv`, `predictions_V.csv`, `manifest.csv`, exact `SOURCE_SHA.txt` and source checksum/provenance. Metrics: logloss, Brier, positive F1, Macro-F1, BA, count of active child regression coefficients for each offset method. The original RLR parent scores can be studied relative to the new ones.

The initial GitHub run [#40495243640](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40495243640) failed for Singh at a **comparison invariant that mistakenly equated active-after-shrinkage pair sets with frozen nominated pair sets**. Distinct fitted offsets are expected to change whether a fixed coefficient is zero; thus the old assertion invalidated the experiment. That failure is not biological evidence. Fixed version explicitly nominates **once**, freezes local pair IDs for all offsets and non-tree controls, and checks exact structural identity before fitting; [updated GitHub Actions workflow](https://github.com/mczajkowskipb/relative-logistic-regression/actions/workflows/rlr-crossfit-offset-DV-v1.yml) is the result to trust when it completes successfully.

**Science go/no-go:** if OOF-offset fit activates more conditional relations but does not improve D-independent V logloss/F1 beyond both original RLR50 and equal-feature global correction, it is NOT a winning model and should not be promoted to full 11/470 evaluation. If it does improve, freeze exact code and evaluate on the full original 11 paper cohort splits using identical real switchBox kTSP, RLR10/20/50 and LR baselines, ideally SHA-checked original processed ARFF from the author's server. Do not reinterpret repeatedly exposed original outer E as fresh confirmation.

## Publication novelty caution

Crossfitted offset construction, sparse logistic residual corrections and local model trees have precedent. New methodological value would require **conditioned pair selection plus an interpretable single-sample rank gate**, useful on independent real datasets relative to matched non-tree, kTSP and RLR baselines. This ablation alone is a model-training consistency check, not a new paper claim.
