# FULL 11-COHORT, ORIGINAL-DESIGN 470-FOLD CFO-RLR TREE — AUDITED FINAL RESULT
2026-10-27. **The original RLR manuscript is UNPUBLISHED and was not modified.**

## Provenance / strict scope

- **REAL, completed GitHub Actions #40501213574**: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40501213574, run SHA `842a3769ead23067ec8b997c124d045a7c472b8a`, status **SUCCESS**, 95/95 jobs passed (94 cohort×repetition jobs with 5 original-design outer folds each, then independent Python audit).
- Python independently replayed **470 E fold IDs ×11 actual comparable model arms = 5,170 scores** from E patient probabilities, true labels and native switchBox predicted classes, with full fivefold-repeat test-patient partition coverage for every cohort. All passed.
- **Original paper-named public reconstructed Schliep expression matrices** with source SHA manifests, **NOT original processed ARFF checksums**: external comparison is research/development evidence, *not proven exact reproduction or an untouched confirmation*. Many cohorts had already influenced exploratory tree design.
- Original study exact **9 datasets×10 repetitions×5 folds + West and Yeoh 2 repetitions×5 folds**. Inferential unit is the **11 datasets**, not 470 correlated test folds.
- Real `SWAP.KTSP` K10/20/50, TSP1, original RLR K10/20/50 d2/m10/alpha.5/lambda.min, conventional value LR, new CFO sparse within-patient rank-only tree (D-only root/feature selection + V acceptance + full T OOF-offset refit), and SAME-PREDICTOR UNION global non-tree correction. Original article RLR and code untouched.

## Dataset-equal-weighted primary outcomes

| Method | Positive-class F1 ↑ | Macro-F1 ↑ | AUC ↑ | Logloss ↓ |
|---|---:|---:|---:|---:|
| kTSP50 | **0.883625** | **0.900266** | **0.962416** | 0.377593 |
| RLR50 ORIGINAL | 0.868613 | 0.897200 | 0.957076 | 0.226026 |
| **CFO conditional RLR tree** | **0.868825** | 0.897004 | 0.956765 | **0.223235** |
| SAME-pair global non-tree | 0.867404 | 0.895848 | 0.956921 | 0.223804 |
| Conventional value LR | 0.829073 | 0.869138 | 0.953675 | 0.269454 |

Other real controls: kTSP10 F1 0.860293 / AUC 0.954313 / LL .378239; kTSP20 F1 .871680 / AUC .958268 / LL .346522; RLR10 F1 .852929 / AUC .949333 / LL .305127; RLR20 F1 .866941 / AUC .955219 / LL .267208.

**Headlines:** CFO tree minus original RLR50 positive-F1 **+0.000212**, Macro-F1 **−0.000196**, AUC **−0.000311**, logloss **−0.002791**. CFO minus actual kTSP50 positive F1 **−0.014800**, AUC **−0.005650**, logloss **−0.154358**. CFO minus equal-feature global positive F1 **+0.001421**, logloss **−0.000569**. Thus **NO meaningful classification improvement** from the tree and **NO demonstrated statistical tree-specific LL superiority**. Original RLR probabilistic advantage over kTSP is far stronger and already belongs to the unpublished original paper.

## 11-dataset paired uncertainty, NOT 470 independent folds

Independent Python exact enumerated two-sided **2^11 sign-flip tests** on equal-weight paired dataset differences (the mean as statistic), plus reproducible 20k dataset-level nonparametric bootstrap 95% intervals (seed 470+metric offset). Exploratory; no multiplicity adjustment. All intervals refer to resampling 11 datasets, not individual patients or correlated E folds.

| Outcome | CFO − comparator | Mean Δ | Bootstrap 95% CI | Exact sign-flip p |
|---|---|---:|---:|---:|
| Positive F1 | RLR50 | +.000212 | [−.002764, +.003313] | .8828 |
| Positive F1 | kTSP50 | −.014800 | [−.069646, +.025076] | .7520 |
| Positive F1 | samepair global | +.001421 | [−.001225, +.004321] | .3906 |
| Macro-F1 | RLR50 | −.000196 | [−.002566, +.001966] | .8594 |
| Macro-F1 | kTSP50 | −.003262 | [−.034569, +.024394] | .8535 |
| AUC | RLR50 | −.000311 | [−.002231, +.001278] | .7969 |
| AUC | kTSP50 | −.005650 | [−.017366, +.006672] | .4004 |
| Logloss | RLR50 | −.002791 | [−.009674, +.003629] | .4551 |
| Logloss | kTSP50 | −.154358 | [−.225984, −.087273] | .0020 |
| Logloss | samepair global | −.000569 | [−.008308, +.007385] | .8691 |

**Interpretation:** there is an important mean LL advantage vs kTSP, but **neither positive F1 nor tree-specific LL differs convincingly**; do not claim statistically significant tree improvement. Mean F1 below kTSP but high heterogeneity and CI across 0 imply we also should not claim significant deterioration.

## Per-dataset positive F1 and logloss differences

Δ is CFO tree minus ORIGINAL RLR50 (negative LL helps).

| Dataset | Real kTSP50 F1 | Original RLR50 F1 | CFO tree F1 | Δ F1 | Δ LL | Tree internally accepted |
|---|---:|---:|---:|---:|---:|---:|
| Alizadeh | .9619 | .9198 | .9178 | −.0021 | −.0064 | 36/50 |
| Armstrong | .9774 | 1.0000 | 1.0000 | .0000 | .0000 | 6/50 |
| Bittner | .9019 | .8653 | .8626 | −.0027 | −.0077 | 12/50 |
| Chen | .8726 | .9096 | .9140 | +.0044 | −.0125 | 44/50 |
| Golub | .9250 | .9331 | .9335 | +.0004 | +.0014 | 36/50 |
| Gordon | .9757 | .9824 | .9824 | .0000 | .0000 | 2/50 |
| Pomeroy | .7121 | .4653 | .4653 | .0000 | .0000 | 0/50 |
| Shipp | .8031 | .8550 | .8566 | +.0016 | +.0062 | 19/50 |
| Singh | .8495 | .9184 | .9182 | −.0003 | −.0025 | 34/50 |
| West | .8505 | .8390 | .8284 | −.0106 | +.0193 | 7/10 |
| Yeoh | .8903 | .8668 | .8784 | +.0116 | −.0283 | 9/10 |

- Internal V accepted **205/470 (43.6%)**, rejected **188**, exact RLR fallback **77**. **3/205 accepted** trees later had **0 active local coefficients** in final T refit (all Armstrong).
- **Crucial failure mode: West:** V accepted 7/10 but test F1 AND LL degraded in aggregate. Internal criterion cannot guarantee E improvement.
- **Pomeroy:** all 50 fallbacks, very small minority class in nested OOF; kTSP50 F1 .7121 versus RLR/CFO .4653, despite RLR AUC .838 (vs kTSP .861). A global fixed .5 RLR probability decision threshold may be part of the minority-class F1 problem, but that is a **hypothesis** requiring train-only testing, not a result.
- Class-F1 gap CFO−kTSP is dominated by Pomeroy (−.2468 in that dataset); the mean over the OTHER 10 cohorts is **+.008396**. This is a **post-hoc sensitivity**, not a new subset-based superiority claim and cannot justify excluding Pomeroy.
- Yeoh local tree improved F1 +.0116 over RLR and LL −.0283; small cohort repetition n=10 and earlier experimental exposure mean caution.
- The historically highest value LR for Yeoh F1 .9517 is a **value-based comparator**, not a rank-only alternative.

## Decision / next research integrity steps

1. **Do NOT proclaim this conditional model tree better than kTSP**. The main objective remains unmet. Publishing tree-specific superiority on this evidence would not be defensible. The original unpublished RLR article remains complete and independently useful for calibrated probabilistic omics predictions.
2. Address two actual bottlenecks as *new scientific hypotheses*: **(a)** small-class support and decision rules (Pomeroy) and **(b)** inner V acceptance that can harm outer performance (West). Do not hyperparameter-fit to observed outer E.
3. A genuinely compact interpretable tree perhaps should use a rank-only anchor-vote gate and **small-sample non-logistic leaf decision** as alternative to unsupported per-branch RLR, with balanced-class support, regularization and internal validation. Compare equal-feature global no-tree and real kTSP. But the formula must be frozen before its next E outcomes and previous cohorts are no longer untouched confirmation.
4. **Full original paper includes PTB fixed, PTB four-cohort LOCO and Legacy8**, which were NOT run here. Original processed ARFF not yet available in GitHub. See `docs/RLR_FULL_PAPER_STUDY_WORLDS_AND_470_CV_CONTROL_2026_10_27.md`. Preserve original train-only MAD2000 for PTB, original cohort LOCO boundaries, and historical eight datasets. Do not silently claim all manuscript worlds completed.
5. Repeat critical confirmatory experiments using SHA-verified **original processed ARFF on the author's server** or external untouched omics datasets. The current reconstructed-source 470-fold analysis is evidence but source parity with original material is not proven.

## Technical artifact

GitHub run: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40501213574

Download `cfo-full470-all11-dataset-means-40501213574` (artifact id 12738443946); files: `audit.log`, `all11_dataset_level_means.csv`, `all_cohort_method_fold_metrics.csv`, `all_cohort_selection_manifest.csv`. Further dataset-level sign-flip/bootstrap report produced from those precise exported artifact files with frozen seed. Original RLR code/manuscript are unchanged.
