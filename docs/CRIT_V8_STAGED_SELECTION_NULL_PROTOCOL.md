# CRIT-v8 Stage D — paired validation selection stress test

Protocol locked 2026-10-25. All source cohorts and folds are **EXPOSED_RETROSPECTIVE_DEV_ONLY**.

## Motivation / evidence from completed stages

- v7: top six D-ranked candidate trees separately tested on same confirmation set V, best selected by conservative validation statistic. Accepted only chen-2002 rep00_fold01 (D rank 3) and fold02 (D rank 6), heldout Δlogloss −0.003734686 and −0.007803831 relative to archived original RLR.
- v8: nominated one candidate from D (D rank 1); all 30 model-arm validation tests across six folds were rejected, therefore 36 model-fold test vectors reproduced the exact original RLR. The complete patient audit passed.
- Stage C: original RLR, D/V sample IDs, D-first candidate gate and leaf relations and candidate 1 V gain, SE and conservative lower statistic were identical in v7/v8. Stage C independent auditor passed all 30 validation decisions. v7's V-selected D ranks were (3,3,6,1,4,5) for chen folds 00–02 and gordon folds 00–02 respectively. The observed difference between v7 and v8 **is caused by their different candidate-selection policies**; it does *not* prove the magnitude of v7 selection bias or that v7's two positive outer-test gains are artifacts.

## New Stage D protocol (model-specific but retrospective)

1. **Reconstruction.** Reuse the complete original v7 RDS files, exact archived original data and fixed original 3-fold inner D/V splits. Only re-estimate the necessary **D-only original-RLR parent and D-only OOF**; no outer-test prediction fitting and no use of E test labels. For each of the six saved v7 candidate designs in each of six folds, use fixed v7 ridge with original D OOF offsets to compute all V patient-specific logloss gains. Assert reconstructed mean gain, SE and conservative acceptance score match the saved candidate confirmation numbers **within 1e−7**. This validates that candidate scores can be decomposed into genuine per-patient contributions rather than inferring covariance from mere means and SE.
2. **Within-fold artificial zero-mean stress test.** Given observed patient × candidate gain matrix `G` (nV×6), subtract the **observed column mean** from each of the six candidates. For each simulation draw choose a shared independent Rademacher sign (±1) for **each patient row**, multiply that row's six centered gains by the same sign, and calculate all six one-SE conservative gain statistics. Centering imposes a zero empirical mean advantage for each candidate; applying the same patient sign preserves dependence across the six alternatives. Fixed seed 20261025+1009×fold index, **5000 draws per fold, 30000 draws total**.
3. Compare three **selection rules under the constructed null**, each using fixed minimum gain >.002 and subtraction of the same fixed cost .002:
   - **D-first single:** test only D-ranked candidate #1, original `mean−1SE−cost`.
   - **V-max six:** test the largest of six `mean−1SE−cost`, matching v7's multiple-candidate mechanism.
   - **V-max with approximate multiplicity penalty:** largest of `mean−z*SE−cost` where `z=Φ^{-1}(1−0.05/6)`. This is a *normal-approximation Bonferroni-style* sensitivity rule; the original one-SE rule is **not** a valid α=.05 test, and this approximate variant is **not** a statistically validated correction for the whole CRIT pipeline.
4. Store **all 30,000 draw-level scores and decisions**, six-fold summary and reconstructed patient gains. Independent Python auditor uses frozen original class orientation to recompute each patient logloss gain, then rechecks all 36 archived candidate score means/SE/lower values, original v7 accept flags, 30,000 draw-level threshold decisions and six summaries.
5. No model selection, inference, cross-study p-values or new scientific claim from this stress test. If the simulated null acceptance of V-max is appreciably larger than D-first, we will have **mechanistic evidence of extra selection opportunity under that artificial null**, not an estimate of true CRIT false-discovery probability. If results look innocuous, that still does **not** establish the safety of V-max under genuine biomedical null datasets.

## Immutable boundaries

- Original published RLR reference SHA `4873878f57f806bf5b52144756f8087b17fa540c`.
- v7 fix3 SHA `e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03`.
- v8 SHA `435b87ae527db82575d4b5cf800725250b53935d`.
- Stage C source SHA `c5708bb55feb35ebfcc8baec89ac2d32d90afacb`.
- All algorithms and dependencies carried forward unchanged; the new code only reconstructs existing validation patient scores and simulates reference distributions. New separate worktree/result directory; no overwrite.

## Output directory and completion criteria

`results/crit_v8_stageD_null_simulation_v1/`:

- `candidate_V_gains.csv` — saved fold/patient × candidate D-rank gains and parent/candidate V probabilities.
- `candidate_score_manifest.csv` — 36 candidate score reconstructions and v7 RDS parity error.
- `null_draws.csv` — 30,000 simulated null draws (first/max/multiplicity scores and acceptance decisions).
- `selection_null_summary.csv` — six sets of 5,000-draw acceptance rates and selection premium.
- `reconstruction.log`, `null_simulation.log`, `independent_audit.log`, `status.txt`; pinned source SHAs and `PROTOCOL.txt`.

Required completion signals:

`STAGED_ARCHIVED_V7_SCORE_RECONSTRUCTION=PASS`,
`STAGED_NULL_30000_PAIRED_SIMULATION_DRAWS=PASS`,
`STAGED_ARCHIVED_V_GAIN_METRIC_PARITY_AUDIT=PASS`,
`STAGED_30000_NULL_DRAW_SUMMARY_AUDIT=PASS`,
`STAGED_BATCH_EXIT=0`.

## Scientific next gate

After interpreting Stage D, **stop empirical tuning on the six exposed folds**. Predeclare a better D-only candidate priority and per-arm accept rules, evaluate exact algorithmic false acceptance and sensitivity on fresh synthetic null/alternative data with genuinely independent test sets, and reserve an untouched biological cohort or batch for eventual confirmatory evaluation. The existing 19 published cohorts and repeated 6-fold exercises remain development data regardless of whether their internal CV seeds change.
