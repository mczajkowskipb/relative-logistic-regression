# RLR pair-specific rank evidence — verified developmental result

**Decision date:** 2026-10-26. **Evidence: PASS; no convincing overall improvement.**

- Source: original locked RLR publication code at `R/new/07_world_protocol_utils.R`, unchanged.
- Experiment branch: `analysis/rlr-pair-specific-rank-confidence-pilot-v1`.
- Original implemented modes: `TERNARY_M10` (identical original rank coding), `ADAPTIVE_Q25` (per-pair train-only threshold), `GRADED_M10_Q75` (rank-gap evidence), K={10,20,50}.
- Run: [successful GitHub Actions #40311987962](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40311987962).
- Artifact: `rlr-pair-specific-v1-40311987962`, ID `12645689714`.
- Cohorts **only from the original RLR paper**: Chowdary-2006 (104 samples, 182 filtered expression features), Singh-2002 (102 samples, 339).
- Six development 3-fold outer splits: one 3-fold repetition per cohort (seed 42). Pair nomination and adaptive parameters used outer training only; no test labels were used during model fitting.
- Nine method settings per split (K×encoding), 54 fits, 1854 saved E patient/method probability rows.
- Independent audit confirmed all 54 study metrics and 1854 heldout patient probability entries; every cohort E test partition is nonoverlapping and complete. `RLR_MARGIN_54_SPLIT_METHOD_METRICS_AUDITED=PASS`, `RLR_MARGIN_1854_E_PATIENT_PREDICTIONS_INDEPENDENTLY_REPLAYED=PASS`, `RLR_MARGIN_E_PATIENTS_PARTITIONED_WITHOUT_OVERLAP=PASS`.
- Status: `RLR_MARGIN_PILOT_EXIT=0`.

## Observed heldout E logloss differences (method minus fixed m10, negative better)

Each entry is the **unweighted mean of six heldout fold-level differences** across the two datasets. “Wins” counts folds with negative logloss difference.

| K | Adaptive-Q25 mean ΔLL | Adaptive wins | Graded-Q75 mean ΔLL | Graded wins |
|---:|---:|---:|---:|---:|
| 10 | +0.0474099 | 0/6 | −0.0017704 | 2/6 |
| 20 | +0.0553414 | 2/6 | −0.0680709 | 3/6 |
| 50 | +0.0735034 | 4/6 | +0.0917927 | 4/6 |

The `GRADED_M10_Q75` K20 apparent mean win is dominated by **one Singh-2002 fold1 ΔLL −0.436190**. Removing that one exposed result for a *sensitivity explanation*, not redefining the experiment, yields ≈**+0.00555** mean ΔLL on the other five folds. On Singh-2002 fold1 at K50, graded rank weighting instead produced **+0.594940** deterioration. The adaptive pair threshold reached mean margins of roughly 25–79 rank positions, depending on the dataset and K, instead of baseline m10; there is no consistently better E probability quality.

## Scientific interpretation / stop decision

1. This experiment addresses the user's original question of **strengthening the existing unpublished RLR article**, *without* constructing a separate XOR/LSH or tree paper.
2. **Do not include pair-specific Q25 or graded-Q75 as claimed superior extensions** on the basis of these data. There is no consistently replicable positive signal across K and cohorts; some cases have large losses.
3. **Do not tune more adaptive percentiles, sigmoid shapes or thresholds on these inspected six outer folds and then call the result independent.** The two datasets are original-paper development cohorts, not fresh confirmation.
4. The current results do not prove that no rank-strength weighting could ever be useful. They show **the two tested formulations** do not meet a practical evidence gate for expanding the RLR article.
5. Preserve original manuscript, baseline RLR, and all frozen results. The existing article's actual scientific emphasis is **probabilistic and calibrated interpretation of pairwise gene-order relations**. Conduct a rigor/novelty/manuscript readiness audit before deciding journal submission.
6. Reuse only the paper's already specified cohorts and protocol. Additional experiments require a **single predeclared mechanism with clear statistical reason**, not another unrestricted hyperparameter sweep.

## Caveats for citing

- This is a small **exploratory** pilot with exposed original-paper cohorts; cross-validation fold results are not independent biological cohorts.
- Exact original `m10` ternary encoding was parity-checked against `encode_pair_matrix()`, but this experiment is a focused RLR coding ablation, not independent replication of the entire published 11-cohort benchmark.
- It changes only rank feature encoding; elastic net fitting with CV may still select different coefficients and λ within the training folds. That is intended and does not constitute an unchanged fitted model.
- The tested `m_j` and graded gaps are simple statistical constructions, not automatically novel methods without a broader prior-art analysis.

**Evidence URL:** https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40311987962
