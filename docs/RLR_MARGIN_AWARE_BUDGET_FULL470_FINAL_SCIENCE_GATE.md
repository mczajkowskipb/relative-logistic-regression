# Full original-design 470-fold audit: margin-aware RLR vs budget-matched real kTSP — SCIENCE GATE

**Date:** 2026-10-27. Original source/method and manuscript are **UNPUBLISHED and locked, untouched**.

**Primary real run:** [GitHub Actions #40589807042](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40589807042): **95/95 SUCCESS** (94 actual original cohort/repeat matrix jobs, five outer folds each = 470, plus independent patient-level auditor). **5640 method-fold summaries replayed from predictions**, original 5fold test coverage per repetition checked, native switchBox K and `K`-budget matching from outer T only checked. Full source: `R/experimental/rlr_margin_budget_full470_onefold_v1.R`, independent `scripts/rlr_margin_aware_budget_full470_audit_v1.py`, original `rg_parent` unchanged. [Reproducible inferential Python code](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-margin-budget-full470-v1/scripts/rlr_margin_budget_full470_inference_v1.py). Original Schliep public reconstructed cohort matrices are **NOT checksum-verified identical to original processed CompCancer paper ARFF**; cohorts already used for method development, so all inference exploratory/non-confirmatory.

**Important averaging distinction:** CI fit log `MB_RLR_MEAN ...` printed 470 equal-*fold* weighted means; these are **NOT** 11 equal-*dataset*-weighted means because West and Yeoh have only 10 folds and other nine have 50. All manuscript-ready science summary below recalculated from the auditor-exported `all11_dataset_mean_matched_budget.csv` assigning equal weight to the 11 datasets.

## Correct equal-weight 11-cohort summary

| Method | Positive F1 | Macro-F1 | AUC | Logloss | Fitted global nonzero rel. | Active local rel./patient |
|---|---:|---:|---:|---:|---:|---:|
| Original RLR20 | **.866941** | .895826 | **.955219** | **.267208** | 16.63 | 15.10 |
| Margin-aware RLR20 | .869632 | **.900086** | .954014 | .274552 | 16.66 | 15.26 |
| Original RLR50 | .868613 | .897200 | .957076 | **.226026** | 28.09 | 25.55 |
| Margin-aware RLR50 | .870235 | .899364 | .959469 | .231992 | 28.21 | 25.76 |
| Real kTSP10 | .860293 | .887541 | .954313 | .378239 | 10.00 | 10.00 |
| Real kTSP15 | .866800 | .890910 | .954437 | .338869 | 15.00 | 15.00 |
| Real kTSP20 | **.871680** | .894456 | .958268 | .346522 | 20.00 | 20.00 |
| Real kTSP25 | .874176 | .895991 | .959438 | .355291 | 25.00 | 25.00 |
| kTSP T-matched global K from original RLR20 | .868213 | .892371 | .957819 | .337768 | 16.63 | 16.63 |
| kTSP T-matched average active K from original RLR20 | .868631 | .893155 | .957551 | .337150 | 15.34 | 15.34 |
| kTSP matched RLR50 globally nonzero K | .879546 | .901182 | .961977 | .354809 | 28.09 | 28.09 |
| kTSP matched RLR50 T-active K | .875223 | .897233 | .960284 | .351685 | 25.83 | 25.83 |

Original RLR20 **must evaluate** about 16.63 globally retained pair relations to discover that ~15.10 of them are nonzero for a specific patient; do not call 15.10 actual tests. These categories are different from gene/probe assay footprint. Real kTSP15 performs 15 pairs, F1 .866800 vs original RLR20 .866941, no significant accuracy difference; mean logloss .338869 vs RLR20 .267208, descriptive advantage in probabilities but 11-dataset paired exact p=.165. For kTSP10 LL .378239 vs original RLR20 .267208, exact dataset-level p=.00879, unadjusted, so still qualify multiple comparisons and kTSP10 cost gap.

## Exact inference on 11 dataset-mean pairs

Two-sided exact 2^11=2048 sign flip test; nonparametric bootstrapped 95% CIs from 30,000 seeded dataset resamples. Primary six F1 comparisons Holm-corrected; all exploratory:

| Paired contrast | Δ mean | Bootstrap 95% CI | unadjusted exact p |
|---|---:|---:|---:|
| MarginRLR20 − original RLR20 Positive F1 | **+.002691** | [−.02030,+.02252] | **.8887** |
| MarginRLR20 − original RLR20 AUC | **−.001205** | [−.00583,+.00355] | **.6152** |
| MarginRLR20 − original RLR20 logloss (negative better) | **+.007344** | [−.01258,+.03205] | **.6367** |
| MarginRLR20 − real kTSP10 F1 | +.009339 | [−.01601,+.03270] | .4639 |
| MarginRLR20 − real kTSP15 F1 | +.002832 | [−.02523,+.02861] | .8555 |
| MarginRLR20 − real kTSP20 F1 | −.002048 | [−.02988,+.02353] | .9092 |
| MarginRLR20 − real kTSP matched original RLR20 train global K F1 | +.001419 | [−.02827,+.02785] | .9199 |
| MarginRLR20 − real kTSP matched original RLR20 train local K F1 | +.001001 | [−.02898,+.02926] | .9443 |
| MarginRLR50 − original RLR50 F1 | +.001622 | [−.01787,+.02257] | .7324 |
| Original RLR20 − real kTSP15 F1 | +.000141 | [−.03131,+.02670] | .9863 |
| Original RLR20 − real kTSP20 F1 | −.004739 | [−.03551,+.02131] | .7412 |
| Original RLR20 − K matched local real kTSP F1 | −.001690 | — | .9268 |
| Original RLR20 − matched local kTSP logloss | −.069942 | — | .1738 |

**All Holm-adjusted PRIMARY F1 comparison p values = 1.0.** No meaningful stable margin-aware F1 or AUC superiority across full11. Firstfold AUC +.022047 for MarginRLR20 vs original was a misleading small-fold signal; 470 repeated folds delta AUC **−.001205**. Robustly abandon this *specific pair-resorting* as a claim of improved classification. Do not retroactively change lambda, support penalty, m, or K on seen E.

## Scientific publication consequence

- **Do not postpone or rewrite original unpublished RLR as if it lost to kTSP50:** K50 should be historically retained in supplement. PRIMARY comparators: real kTSP10/15/20, kTSP training-budget-matched to original RLR global nonzero + patient active K. In all plots explicitly distinguish nominal shortlist, global nonzero `β`, patient-active `r≠0`, actual tested comparisons, unique genes.
- Existing original RLR20 has similar positive F1 to real kTSP15, and better mean probability loss (though 11-level LL signflip p=.165); this is much cleaner research framing than cherry-picking 50-pair extremes. Original RLR50 has substantive calibrated probabilistic advantage vs 25–30 comparable kTSP panels, at cost of more tested relations. Do not claim statistically proven F1 equivalence from p>0.05 or automatic clinical benefit.
- **No NEW decision tree success** resulted from the margin-aware global selector, earlier SAVS1–4 small conditional trees or sparse 4-rule voting. The older certified exact-prediction early-exit interpretation remains technically credible but its learned contextual tree did not outperform fixed coefficient ordering; don't rebrand it as novel tree learning without new evidence.
- Next research should target training-only improvements that are not an after-the-fact choice among exposed E baselines; e.g. **cross-fitted sparse ordinal margin panel selection** with explicit non-tree and simple kTSP controls, potentially model rule sparsity/stability instead of retuning m. For a habilitation-targeted publication, reserve tree claims for a truly conditional, independently better structure, or use original completed RLR as companion relational non-tree model to existing published tree works.
