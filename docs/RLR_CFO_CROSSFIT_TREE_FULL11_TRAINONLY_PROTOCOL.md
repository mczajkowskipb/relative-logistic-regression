# CFO-RLR model tree across 11 original CompCancer cohort names — inner D/V protocol V1

**Study type: discovery/training-only mechanism feasibility, NOT external confirmation.** Original RLR manuscript remains unpublished and untouched. Code and fixed hyperparameters were frozen before this 11-cohort run was launched. Several cohorts, especially Yeoh, had been studied previously, so even this 11-cohort D/V screen is still exploratory and cannot be called a new untouched sample test.

## Why

Original RLR K50 has favorable logloss on published-source reconstructions but in the unpublished manuscript's complete CompCancer analysis its AUC and positive-class F1 means slightly trail real kTSP50. Earlier branch-RLR modifications produced inconsistent/nonexistent local coefficients when fitter uses parent **in-sample** probabilities as offset. Switching only to a 3-fold crossfitted parent offset in D dramatically changed local coefficients and produced a true conditional logloss gain beyond an equal-feature global model on Yeoh's inner V. Chen had degradation and Singh's global same-feature control was better; therefore need balanced across-11 proof before any expensive 470-split study.

## Model and original samplewise invariance

- Original global RLR50 uses supervised top1000 switchBox candidate selection within each training fold, strict original gene-disjoint `d=2`, `m=10` three-level within-patient rank comparison encoding, binomial elastic net `alpha=.5`, original `cv.glmnet` `lambda.min`. Frozen original helper from `R/new/07_world_protocol_utils.R`.
- Root gate `1[sum_{a in A}1(r_g(x)>r_a(x))>=t]`, A 3 or 5 D-only selected reference anchors, up to one anchor swap per position and up to 240 unique proposals. D-only crossfitted original RLR50 residuals score conditional root; <=6 roots evaluated.
- Child local feature IDs selected **ONCE** per root and branch from D-only rank-ternary candidate features, using local OOF residual gradient minus global main-effect proxy, at most eight disjoint pairs per branch.
- Fit branch-specific penalized logistic **correction** with original RLR logits as offsets; **offset for D fitting uses true 3-fold D-OOF parent RLR50**; fixed correction alpha=.5/lambda=.10; branch may be zero and then patient receives exact original parent score. For prediction on unseen V, use original parent RLR50 fit on D.
- Fit a **same-nominated-pair non-tree correction**, using identical D OOF parent offset and the union of **all** nominated local pairs, regardless of which coefficient was shrunken to zero in the branch models. Equal-feature control is mandatory, not optional.
- No raw expression values, cross-patient abundance ranges, expression ratios, learned gene-dependent calibrations or cohort-scale normalization at prediction. Same fixed assay probes/universe required. Any strictly increasing transform of all values in a single patient preserves ranks, gates and conditional pair features.

## Partition and model choice

For each of the exact original 11 paper **cohort identities**, the 5-fold stratified outer split `rep00_fold00` from original protocol (seed42) is taken. **External E patients and outcomes are withheld entirely**, no E metric is calculated. On outer T only, pre-existing `rg_inner` defines D (~75%) and V (~25%). Small cohorts whose D or 3-fold crossfit subtraining sets cannot support at least four internal glmnet folds produce explicit original-RLR safe fallback, NOT a spurious local model or crash.

On D, build up to six fixed candidate root/child-pair sets. Evaluate original RLR50, branch conditional RLR50 with OOF-offset-fitted local corrections and global samepair correction **on identical V subjects**, preserving probabilities. Choose minimum V logloss root (this choice makes its V loss optimistic; **no outer E validation yet**).

One trained root is accepted for further research only if all four validation conditions are met:

1. V original RLR50 LL − tree LL >= **0.005**.
2. V global same-pairs correction LL − tree LL >= **0.002**.
3. Tree V BA − original RLR50 V BA >= **−0.02**.
4. Tree V Macro-F1 − original V Macro-F1 >= **−0.01**.

The first three gates reuse previous SA-RLR training-only protocol. Macro-F1 criterion was a **post prior cohort analysis safety response to observed loss of classification in Singh**, introduced BEFORE the new 11-cohort run but NOT independent of the development evidence. All of these are engineering/discovery choices; **not pre-registered external tests**. Even accepted cases should not be presented as out-of-sample superiority unless frozen method survives original full independent outer E protocol.

## Real execution and audit

[GitHub Actions 11-cohort run #40498147762](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40498147762) — 11 separate full original-cohort first-split D/V jobs plus independent 11-result patient-level auditor. Output per cohort: `SOURCE_SHA.txt`, original RLR helper SHA, source manifest, `selection_DV.csv`, and (when roots fitted) `patient_V.csv` containing ALL six root designs×five methods×inner V patient probabilities. Auditor recomputes LL, Macro-F1, BA, exact root chosen and all acceptance rules from patient-level V predictions. It aborts if even one result is a genuine computational error.

**Stop rules:** If only the previously seen Yeoh or other already exposed cohorts improve, or gains are confined to logloss with no discrimination benefit / absent equal-feature global superiority, DO NOT claim a publishable enhanced tree. The full paper goal is to beat original real kTSP in positive F1/AUC, preserve RLR probability quality, and have a stable interpretable patient-specific tree rule on original 11 exact processed ARFF datasets and 470 originally used outer folds. These GitHub public source reconstructions do NOT establish ARFF parity.

## Prior evidence not to misquote

- Completed [3-cohort D/V double ablation](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-crossfit-offset-local-experts-v1/docs/RLR_RESIDUAL_PAIR_AND_CROSSFIT_OFFSET_DV_COMPLETE_AUDIT.md) (Chen, Gordon, Singh): independent 162 patient-derived metrics; OOF offset activates features but does not produce consistent tree/ F1 advantage over global controls.
- Completed [Yeoh D/V](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-crossfit-offset-yeoh-dv-v1/docs/RLR_YEOH_CROSSFIT_OFFSET_CONDITIONAL_TREE_DV_COMPLETE_AUDIT.md): 49 V, six D-nominated root candidates, identical local pairs; tree OOF logloss mean .00350969 vs parent .01161792 and samepair global OOF .00967637, **tree better in all 6 roots**. However Macro-F1=1.0 for every arm, so no classification improvement, and six candidates are correlated.
- Previous 11 cohort FIRST heldout E fold threshold-decision pilot showed mean positive F1 real kTSP50 .865621, original RLR50 .867924; conditional threshold tree **0/11 accepted**, global threshold worsened average F1. It is not a true full 470-fold reproduction.
