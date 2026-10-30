# Controlled RLR vs binary pair comparators — five outer folds / 11 cohorts, v2

Scope: exploratory comparative evidence on **public Schliep reconstructed** datasets, not verified identical to paper-processed original ARFF.

## Motivation and controls
The prior one-split 11-cohort D/V pilot favored original RLR20 numerically but did NOT establish an advantage versus binary elastic net (p=.1543). This v2 protocol is fixed before obtaining any v2 outer-fold results, but the design is not blind to v1 findings. Original manuscript code and frozen paper data unchanged.

## Design
One full fixed-seed original 5-fold outer CV repeat on each of 11 Schliep source cohorts, using original stratification seed42. All outer fold test E samples are untouched until predictions. Within each outer training T, run D-only switchBox top1000 supervised pair ranking, strict gene-disjoint shortlist d2 K20, then the identical selected relations for three arms:
- ORIGINAL_RLR20: original frozen source rg_parent, m=10 signed rank-gap ternary, glmnet alpha=0.5, 5-fold training-side CV and lambda.min.
- BINARY_ENET20_PURIST_INSPIRED: binary raw gene pair ordering with glmnet alpha=0.5, same CV/lambda.min, constant train pair columns dropped for both binary arms.
- BINARY_MCP20_DECAF_INSPIRED: binary ordering and ncvreg MCP penalty alpha=0.05, same CV fold budget, as in DeCAF public training helper.

These latter two controls are **method-inspired ablations, not independently reproduced complete PurIST/DeCAF clinical models**. Particularly DeCAF originally filters genes by cross-study class-stable rankings and fits MCP; this benchmark deliberately holds pair selection equal to isolate binary/ternary and penalization.

## Outputs and inference
For each of 11 cohorts, five folds and three arms, record 165 method-fold summaries and every out-of-fold sample-level probability. Report dataset-equal means for logloss (prioritized), Brier, AUC, F1 and ECE. Exact paired two-sided signflip 2^11 permutations on dataset-level fold means, not 55 independent folds, and exploratory BH across eight paired comparison/metric tests. A separate Python auditor must independently verify all patient scores/metrics and test-split disjointness.

## Caveats and success criteria
Do not assert dataset input parity before validating SHA256 against published paper manifest. All results, including unfavorable and null, must be reported. No manuscript selection/updates based solely on pilot results, and no post hoc hyperparameter change. Findings remain exploratory because the previous pilot and extensive publication benchmark have already been analyzed.

## Future canonical replay
The **same code can be adapted** to use exact original ARFF only after checksums pass, with training/test row/column/class mapping checked first. No direct use of original ARFF in GitHub Actions is claimed.
