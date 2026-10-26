# RLR original-manuscript component evidence: same pairs, equal votes, learned coefficients

**Audit date 2026-10-26 — conservatively negative/mixed component-specific result.** Original locked RLR code and manuscript in `main` have not been modified.

## Source inventory correction

Only the following 11 datasets are in the original publication's primary CompCancer inventory: `alizadeh-2000-v1`, `armstrong-2002-v1`, `bittner-2000`, `chen-2002`, `golub-1999-v1`, `gordon-2002`, `pomeroy-2002-v1`, `shipp-2002-v1`, `singh-2002`, `west-2001`, `yeoh-2002-v1` (`docs/DATA_ACQUISITION.md`). Historic exploratory Chowdary-2006 experiments **are not paper-inventory experiments**; do not report them as a valid subset of the original 11.

The corrected experiments ran strictly on original-inventory `armstrong-2002-v1` (72 samples, 1081 filtered probes) and `singh-2002` (102 samples, 339 probes), using source matrices downloaded from the same Schliep dataset collection as the original benchmark. **Exact equality of the reconstructed source-text matrices to the locked processed ARFF checksum files has not been independently proved.** Thus this remains a **paper-inventory, exposed development pilot** rather than exact reproduced full publication benchmark.

## Reproducible algorithm and evidence

All methods use the SAME supervised `switchBox` top-1000 pair pool nominated from outer training only, the SAME greedy strict gene-disjoint `d=2` selected pairs at K=10,20,50, and the SAME original signed ternary rank margin m=10. All training and regularization selection are inside the outer training set; six heldout outer E folds: one 3-fold stratified repetition per dataset. Source `R/new/07_world_protocol_utils.R` is unchanged.

Comparators:

1. **ORIGINAL_RLR:** learned independent per-pair elastic-net coefficients via `glmnet::cv.glmnet`, alpha .5, 5 inner CV folds, lambda.min.
2. **CALIBRATED_EQUAL_VOTE (V2):** compute `v_i=mean_j r_{ij}` with *identical* original ternary features; fit conventional unregularized one-dimensional logistic Platt scaling of v on the outer training set, then evaluate E. This is a legitimate simple calibration baseline, but it differs from original RLR in regularization.
3. **REGULARIZED_EQUAL_VOTE (V3):** SAME `v_i`, trained by the SAME `glmnet::cv.glmnet` alpha .5, 5 inner folds and lambda.min. The implementation supplies one varying vote column and a second **exact constant zero column** to satisfy `glmnet`'s two-column API requirement; confirms the second coefficient remains zero. The predictor therefore has exactly one learned vote slope and intercept. This comparison much better isolates whether learning independent relational coefficients provides a substantial advantage beyond regularization of a single relational vote score.

The V2 and V3 experiments are separately versioned, complete GitHub-hosted runs with complete source SHA and independent patient-score recomputations:

- [V2 Platt PASS #40320814865](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40320814865), evidence artifact `rlr-same-pairs-paper-v2-40320814865`, ID `12648972399`.
- [V3 equal regularization PASS #40321243077](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40321243077), evidence artifact `rlr-same-pairs-fair-v3-40321243077`, ID `12649440935`; `SAME_PAIRS_36_METHOD_STUDIES_REPLAYED=PASS`, `SAME_PAIRS_1044_PATIENT_E_SCORES_REPLAYED=PASS`, `SAME_PAIRS_ALL_E_SPLITS_DISJOINT_AND_COMPLETE=PASS`. Each run has 36 model×fold studies and 1044 E patient×model probability records.

## All heldout Δlogloss = equal-vote minus original RLR (positive favors RLR)

| K | Unregularized Platt vote ΔLL | Folds RLR better | Identically regularized vote ΔLL | Folds RLR better |
|---:|---:|---:|---:|---:|
| 10 | +0.23771409 | 6/6 | +0.00626566 | 4/6 |
| 20 | +0.21260464 | 4/6 | −0.01920431 | 4/6 |
| 50 | +0.24962449 | 5/6 | +0.03145070 | 3/6 |

Unweighted mean of six outer fold values (3 from Armstrong, 3 from Singh), **not independent cohort-level replicates**. Large losses for Platt on Armstrong folds 1 and 2 illustrate its unregularized model can become overconfident on limited data; they must not be interpreted as a causal advantage of independent weights. Under the equal-regularization comparison, apparent RLR advantage is much smaller and not consistent across K and folds. For K20, regularized equal voting strongly improves Singh fold1 (ΔLL −0.348630), dominating the negative mean.

## Manuscript decision

1. **Do not use unregularized Platt differences alone** to assert learned per-pair coefficients cause the original observed kTSP/RLR logloss improvements.
2. **Do not claim a validated contribution from independent weights in isolation.** A properly regularized equal-vote model is competitive within this exposed two-dataset pilot.
3. Original RLR remains a legitimate *integrated* method combining learned pair nomination/orientation, margin-aware within-patient rank relations, strict disjointness, elastic-net fitting and inspectable probabilities. The existing main 11 CompCancer, PTB fixed/LOCO and supplemental calibrated-kTSP analyses still require critical review as an integrated pipeline claim.
4. A true component-level article claim would require running the same-pairs fair regularization ablation under the **frozen original 11 CompCancer resampling protocol** (9 datasets × 10×5 plus 2 datasets × 2×5), with original processed ARFF SHA verification. The current GitHub CI experiments have not done that.
5. STOP further feature/margin/tree hyperparameter search on reused E partitions. Better next publication work: original manuscript novelty/fair-baseline audit, independent full-world reproducibility check, and journal readiness. Do not silently change source or paper based on exploratory pilots.

This note documents the evidence in the project workspace. **Not submission-ready inference.**
