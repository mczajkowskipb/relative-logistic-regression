# RLR vs binary logistic pair variants — frozen prospective pilot, 2026-10-30

## Hypothesis
Is original RLR20's margin-aware ternary rank encoding beneficial over two existing-method-inspired binary pair models when pair choice and training set are held fixed?

## Source
Eleven **public Schliep reconstructed** CompCancer cohorts; names and sample counts match the paper. **Exact original publication ARFF SHA256 parity NOT VERIFIED.** The pilot must not be substituted for the main paper benchmark.

## Analysis world
Use the original repeated stratified 5-fold seed42 schedule, **only first original outer training T**. Divide T into D/V using existing `rg_inner` seed420. Supervised selection, coefficients, all hyperparameter decisions learned using D only. Outer E remains unread. One split per dataset, 11 independent dataset units. Do not retune after viewing results.

## Matched shared-pair experiment
D-only original `compute_top_pairs_switchbox` top1000 candidate ranking, `apply_disjointness` strict d2, K20. The same ordered gene pairs are used for all models.

- **RLR20_ORIGINAL**: original `rg_parent` from frozen original source; ternary within-patient rank-gap m10, glmnet elastic net alpha .5, training-only CV and lambda.min.
- **BINARY_ENET20_PURIST_INSPIRED**: binary gene pair indicator, glmnet elastic net alpha .5, training-only CV. This is an *adaptation* inspired by PurIST's weighted binary pairs, **NOT** a reproduction of original PurIST fitted pancreatic model or its full learning algorithm.
- **BINARY_MCP20_DECAF_INSPIRED**: binary gene pair indicator, ncvreg penalty MCP alpha .05, training-only CV. Based on published DeCAF training logic; **NOT** original DeCAF fitted CAF classifier, its gene selection or original disease labels. Constant D predictors removed from BOTH binary arms and counted.

## Endpoints
Primary exploratory descriptive: mean of 11 dataset-level validation logloss values (smaller better). Secondary: Brier, positive-class F1, AUC, ECE, number of active coefficients. Independent Python auditor replays all patient-level logloss and Brier and computes exact two-sided 2^11 dataset sign-flip tests for paired RLR vs each binary arm.

## Reporting rules
Always report all arms, all 11 datasets, all metrics including unfavorable comparisons. No claim of preregistration before prior RLR results; this comparator experiment is specified prospectively **before its own predictions**. Do not present any positive pilot outcome as confirmatory novelty or publication-ready without original ARFF parity and independent holdout evidence. Do not overwrite original manuscript or publication results. No changes to historical Git commits.

## Original method references
- https://github.com/naimurashid/PurIST
- https://github.com/jjyeh-unc/decaf
- https://github.com/jjyeh-unc/decaf_manuscript (includes MCP training code)
