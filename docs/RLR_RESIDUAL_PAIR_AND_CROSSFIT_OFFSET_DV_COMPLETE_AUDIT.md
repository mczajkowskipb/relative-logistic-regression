# RLR local conditional pairs — COMPLETE 3-cohort D/V mechanism audit (2026-10-27)

**Overall status: IMPLEMENTED, RAN, and INDEPENDENTLY REPLAYED; NO SCIENTIFIC MODEL-TREE WIN.**

The original RLR manuscript is **UNPUBLISHED** and unchanged. These experiments intentionally use **D/V partitions within original outer training T** only, no original outer E labels or outcomes. They examine whether local RLR corrections are worth pursuing, not clinical or paper-level independent test accuracy.

## Sources, immutable algorithm, real CI

- [Residual-selected versus class-TSP-selected branch pairs, GitHub Actions #40492939641](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40492939641) **3/3 jobs SUCCESS**; branch [`analysis/rlr-residual-conditioned-child-pairs-v1`](https://github.com/mczajkowskipb/relative-logistic-regression/tree/analysis/rlr-residual-conditioned-child-pairs-v1).
- [Residual-pair branch experts trained with OOF versus in-sample RLR50 offsets, GitHub Actions #40495707887](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40495707887) **3/3 jobs SUCCESS**; branch [`analysis/rlr-crossfit-offset-local-experts-v1`](https://github.com/mczajkowskipb/relative-logistic-regression/tree/analysis/rlr-crossfit-offset-local-experts-v1), fitted commit `97c60d1e104456eca49b5ee75a7272fb7dd416bd` (subsequent documentation does not modify model).
- Independent **Python CSV replay directly from the six GitHub archived artifacts**: verified 3 cohorts × 6 root candidates × 4 comparison arms = **72** residual-pair fold-arm scores; plus 3 × 6 × 5 = **90** matched-offset fold-arm scores, **162 records total**. Recomputed *from each patient prediction* every logloss, Brier, positive F1, Macro-F1, balanced accuracy, exact original parent parity across six roots, class identity, V sample count and per-root method alignment. All numeric metrics match archived R outputs to absolute error <1e-9.
- Original gene schema and original RLR K50 (`d=2,m=10,alpha=.5,lambda.min`) untouched; conditional correction glmnet `alpha=.5,lambda=.10`, max 8 gene-disjoint ternary rank relations per child. All feature tests are within-patient rank comparisons; no feature value cutoffs, expression ratios or normalization ranges.
- Input matrix is reconstructed public Schliep source for Chen-2002, Gordon-2002, Singh-2002; **original processed paper ARFF SHA parity still NOT verified**. A single frozen 5-fold outer train split per cohort, D/V inner split. Six roots chosen from D-only anchor one-swap candidate screen per cohort; identical local ranking and fit comparisons on same V patients. **Do not treat six roots as six independent patient cohorts or statistical replicates.**

## Experiment 1 — choose local rank relations by ORIGINAL RLR OOF residuals vs local class TSP scores

Second-order branch-specific score `G_branch²/[2(H_branch+λ)]` minus proportionate whole-D global feature gain, with original RLR50 K50 3-fold-D OOF residuals. No E use. Train both locally using in-sample parent offset to isolate only the nomination effect.

| Dataset | Original RLR50 V logloss | Mean LL of class-TSP-selected child corrections | Mean LL of residual-selected child corrections | Residual better than class, roots |
|---|---:|---:|---:|---:|
| Chen | .164126 | .164126 | .164126 | 0/6 |
| Gordon | .007420 | .007420 | .007420 | 0/6 |
| Singh | .151361 | **.144263** | .146623 | 1/6 |

Both methods fitted no active branch correction in Chen/Gordon at fixed λ=.10. In Singh both improved mean V logloss against original RLR but class-TSP-selected pairs performed better than the residual-selected local tree. No macro-F1 gain: Singh parent and both local variants Macro-F1=1.0 on the selected inner V fold.

**Decision:** selection of local pairs based on this crossfitted residual proxy is technically valid but has NOT demonstrated a predictive advantage over existing class TSP nomination on these three development cohorts.

## Experiment 2 — SAME frozen residual-selected root and pair IDs, different OFFSET at child training

Scientific defect detected in initial implementation: features were nominated using original RLR50 **OOF residuals** but conditional glmnet coefficients fitted with a different, more optimistic **in-sample parent logit offset**. New matched ablation freezes each root and local reference gene IDs on D once, then only changes child fitting offset. Global no-tree control uses **exact same full nominated pair union**, not just post-penalty surviving pairs.

| Cohort | Original RLR50 LL | Local, full-D in-sample offset LL | Local, D OOF offset LL | Global SAME pairs, D OOF offset LL |
|---|---:|---:|---:|---:|
| Chen | .164126 | .164126 | **.193778** | .168764 |
| Gordon | .007420 | .007420 | .007410 | .007420 |
| Singh | .151361 | .146623 | .126708 | **.125139** |

Mean over six root designs in each cohort, evaluated on identical heldout *inner V* patients. Smaller is better.

- **Chen:** Child coefficients activate with OOF offset where old full-D offset yielded zero, but V logloss deteriorates vs parent and vs global same-pairs control; no Macro-F1 improvement (.97105 for all arms). Thus extra nonzero coefficients do **not** establish useful local corrections.
- **Gordon:** Numerically tiny average LL gain (.0000099), 1/6 roots locally improved, baseline Macro-F1 already 1.00. Not meaningful clinical or paper-level evidence.
- **Singh:** OOF offset substantially improves mean V LL .151361→.126708 relative to parent; but the **non-tree same-pairs global correction yields .125139** and much better Macro-F1 .99165 vs local .95798; parent .151361 has Macro-F1 **1.0000**. Although in 3/6 individual root designs the local OOF model beats its global comparator, the average is inferior to the global comparator. This is NOT reliable evidence that routing by an anchor stump adds benefit.
- Active coefficients summed over six candidate trees at λ=.10: Chen full-D 0 / OOF 84, Gordon 0 / OOF 8, Singh full-D 9 / OOF 83. These are **counts across many candidate fits**, not numbers of distinct genes, not stability statistics and not fit coefficients in one final accepted classifier.

### Technical integrity correction

Initial [run #40495243640](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40495243640) incorrectly compared **post-shrinkage active pair lists** between the two offset models, causing a false assertion failure in Singh: a different offset should reasonably zero coefficients. This was NOT clinical evidence. Fixed [run #40495707887](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40495707887) strictly **nominates pair IDs ONCE before both fits**, freezes those same pair IDs, and constructs the same-feature global controls from the **unfiltered nominated list**. All three final CI jobs PASS.

## Publication decision

**STOP this specific combination at development gate: neither class-residual nomination nor OOF parent-offset adjustment has shown a robust, F1-improving TREE benefit over same-feature no-tree controls.** The original complete RLR paper should not be delayed indefinitely for such negative experiments. Especially, Singh with V Macro-F1 1.0 for original RLR cannot evidence improved classification by more local coefficients.

The next only defensible hypothesis would address manuscript's still weak datasets (not merely the nearly perfectly separable Gordon/Singh), and use a clinically understandable *rank-invariant selective conditional comparison mechanism* on **D/V only** with a prespecified limited search. Candidate dataset **Yeoh** is already present in the original manuscript (RLR50 vs kTSP50 F1 gap) but has been inspected in previous experiments, so any result from it is strictly exploratory. Use original exact processed ARFF when available before claiming reproducible original paper numbers.

A full 11-cohort × original 470 repeated folds and real kTSP/LR RLR comparisons remains the eventual **confirmation** gate, after one tree variant shows honest D/V predictive benefit and we verify original processed ARFF checksums. We are not entitled to claim beating kTSP/ LR on the basis of any of these D/V mechanics.
