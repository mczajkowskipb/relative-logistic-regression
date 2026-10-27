# SA-RLR train-only full-11 mechanism audit — interim progress and decision gate

Date: 2026-10-27, Europe/Warsaw. **This is NOT an efficacy or heldout E validation report.**

## Code and actual reproducible GitHub runs

1. New 11-cohort, canonical first 5-fold outer split, **T-only/D-V only** original manuscript research procedure:
   - Workflow: [#40476059174](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40476059174)
   - Model: `R/experimental/rlr_selective_anchor_gated_v1.R` (SA-RLR v1)
   - Runner: `R/experimental/rlr_sa_trainonly_full11_feasibility_v1.R`
   - Independent 11-cohort artifact auditor: `scripts/rlr_sa_trainonly_full11_aggregate_v1.py`
   - True original paper cohort identifiers and original split assignment method; public Schliep reconstructed sources rather than checksummed original processed ARFF.
   - E outcomes not provided to any model. Full 470 outer-fold benchmark has NOT been executed.
2. Independent local-expert penalty-path mechanism diagnosis on the **D part only** for Chen, Gordon, Singh: [SUCCESS #40476671294](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40476671294).
3. Technical sample-support fix and regression test: [SUCCESS #40477269866](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40477269866).

## Preliminary D-only findings, no outer E assessed

- **Chen (D n=108):** 240 eligible *unique* ordinal root partitions generated; 225 differed from baseline reference-gene panels through actual one-anchor exchange. Inner D/V search found no candidate with acceptable active conditional RLR. Local D-only penalty-path check over six screened roots/twelve branches: 0 nonzero fitted coefficients with `lambda=.10`, 6 with `.03`, 46 with `.01` (sums over local fits, not a final selected tree).
- **Gordon (D n=108):** 240 roots, 205 actual swaps, but all 240 partitions only had ONE branch with sufficient mixed-class D support; no inner V-accepted model. In the D-only six-root/twelve-branch diagnosis, 0 active coefficients for `.10/.03/.01`; lowering the penalty alone does not create useful local coefficients for these candidate branch/pair choices.
- **Singh (D n=61):** 240 roots, 212 actual swaps, 123 with only one class-mixed viable branch and 117 with two viable branches. Inner V **accepted one** one-branch model with **3 active local rank relations**, V logloss improvement over D-trained RLR50 **+.0093202** (defined baseline LL minus SA LL) and same-pairs-global control **+.0093202**. No outer E data used. In the separate D-only six-root/twelve-branch diagnostic, total active coefficients 19 at `.10`, 47 at `.03`, and 60 at `.01`. This demonstrates some local residual capacity but **does not prove generalization**.
- **Armstrong:** 74 unique root proposals, 52 swaps; 0 V-selected conditional local models.
- **Alizadeh:** no eligible root partitions with the predeclared branch-support criteria, exact original RLR fallback.
- **Pomeroy:** first pre-fix pilot returned a controlled **diagnostic ERROR**, not a negative biological result: the inner D group was too small for `cv.glmnet`, yielding `nfolds must be bigger than 3`. GitHub artifact includes the exact error in `train_only_selection.csv`. The algorithm now explicitly returns exact original RLR fallback if either D or an OOF subtraining fold would have <4 folds via `safe_nfolds`; the focused deterministic mechanical R test for this guard PASSED. This fix does not alter original RLR or tune based on any E outcomes. The previously running cohort matrix was launched before the fix; its archived Pomeroy entry must be labelled **pre-fix**.

Do not claim all 11 are complete until the relevant GitHub workflow run and independent aggregation have completed successfully. Do not reinterpret the independent path diagnostic as an outer-E search, or present smaller lambdas as prevalidated effective hyperparameters.

## Implications for next scientific choice

1. SA-RLR's anchor swaps are generating diverse **patient partitions** on suitable training cohorts. The method is not failing simply because all root panels are duplicates.
2. The key bottleneck is **conditional relational expertise beyond RLR50**, strongly dataset-dependent. Some candidate branches have too few mixed-class samples (Gordon), others are penalized to zero at the frozen .10 penalty (Chen) and some have nonzero informative local pair candidates (Singh).
3. A single fixed conditional logistic penalty may be too rigid. Before changing penalty or tree depth on the same exposed data, require D-only inner validation across all 11 cohorts to reveal whether any setting reproducibly improves BA/F1 rather than only LL.
4. Another *distinct, complementary* publication-oriented angle is train-only crossfitted **decision-threshold selection** for original RLR50 and a rank-anchor-conditional decision threshold. Such postprocessing keeps the original continuous RLR probabilities, logloss, and AUC unchanged but **can change positive F1/BA**. It must be compared against a global tuned threshold at equal patient folds, so any tree-specific merit is demonstrable rather than falsely attributing ordinary threshold calibration to tree learning. Do not claim that this is novel merely because it uses a tree: cost-sensitive decision rules/threshold calibration are prior art.
5. Original paper compiled on the canonical 11 cohorts: RLR50 AUC .957 / positive-class F1 .869 / LL .226, kTSP50 .963/.884/.378. The scientific target remains discrimination (F1/AUC) while not losing RLR's large probability-calibration advantage. Performance of a new method must be measured on the full 11/470 original protocol, ideally checksummed original processed ARFF on Ubuntu, and at dataset level (11 independent datasets), not by selecting the best repeated fold outcome.
