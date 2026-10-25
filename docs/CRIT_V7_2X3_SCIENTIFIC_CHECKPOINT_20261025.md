# CRIT-v7 2×3 scientific checkpoint — 2026-10-25

## Provenance and status

- Frozen experiment commit: `e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03` (`feature/crit-v7-systemd-env-fix3`).
- Operator Ubuntu output directory: `/opt/2026/RLR/relative-logistic-regression-crit-v7-fix3/results/crit_v7_original_rlr_2x3_v1`.
- This note was transcribed from the operator's Oct 25 terminal outputs, with code-methodology cross-check against the pinned commit. It is **not** a second independent replay of the Ubuntu files.
- `CRIT7_REAL_SYNTHETIC_ORIGINAL_ROOT_AND_TREE=PASS`, `CRIT7_PILOT_ALL_24_METHOD_FOLDS=PASS`, `CRIT7_PILOT_EXACT_ARCHIVED_ORIGINAL_RLR=PASS`, `CRIT7_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS`, `CRIT7_24_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS`, `CRIT7_ALL_GATES=PASS`, `CRIT7_PILOT_BATCH_EXIT=0`.
- `systemd` unit exited cleanly, without an active CRIT process.
- Scientific designation: **EXPOSED_RETROSPECTIVE_DEV_ONLY**. These cohorts and fold choices have already been used during RLMT development.

## Locked original-RLR anchored evaluation

Exactly 2 datasets × 3 outer splits × 4 methods = **24 audited patient-level prediction vectors**.

| Method | Accepted trees / 6 | Fold wins/losses/ties vs original | Mean fold Δlogloss vs original |
| --- | ---: | ---: | ---: |
| ORIGINAL_RLR | 0 | 0 / 0 / 6 | +0.000000000 |
| GLOBAL_CAL | 0 | 5 / 1 / 0 | −0.000251335 |
| GLOBAL_INTERACTION | 0 | 3 / 3 / 0 | −0.002716681 |
| CRIT_TREE | **2** | **2 / 0 / 4** | **−0.001923086** |

Negative delta favors candidate. All means are *unweighted arithmetic means of six fold-level logloss deltas*, not pooled patient-level logloss and not independent-cohort evidence.

### CRIT_TREE observed split decisions

| Dataset | Outer split | Decision | CRIT Δlogloss vs RLR |
| --- | --- | --- | ---: |
| chen-2002 | rep00_fold00 | reject (exact RLR fallback) | 0.000000000 |
| chen-2002 | rep00_fold01 | accept | −0.003734686 |
| chen-2002 | rep00_fold02 | accept | −0.007803831 |
| gordon-2002 | rep00_fold00 | reject (exact RLR fallback) | 0.000000000 |
| gordon-2002 | rep00_fold01 | reject (exact RLR fallback) | 0.000000000 |
| gordon-2002 | rep00_fold02 | reject (exact RLR fallback) | 0.000000000 |

No generalized effect across datasets has been established. CRIT's two positive fold results occur only within `chen-2002`. The independent archive parity max absolute difference reported for each fold is of order 10^−16.

On the two accepted chen folds, the reported global interaction deltas were −0.012521440 and +0.004055763, respectively. Thus which model is preferable is fold-dependent. Across all six folds, the global-interaction *mean* delta is more favorable than CRIT's, but its 3/3 win/loss count is less consistent.

## Scientific interpretation and caveats

1. **Feasibility, not superiority.** The synthetic planted XOR test passed. Biomedical trees were accepted in 2/6 folds and improved their respective outer-fold logloss vs original RLR. This is evidence that the implementation sometimes produces useful context-specific fits on retrospectively exposed data; it is **not** external confirmation.
2. **Non-equivalent comparators.** The current `GLOBAL_INTERACTION` fits one screened relation and its product with the gate, even if the tree fails heldout confirmation. CRIT uses separately chosen leaf relations, but falls back exactly to RLR unless its heldout test approves. Therefore the 6-fold average compares different model classes *and different acceptance policies*.
3. **Function-class equivalence of a one-split tree.** Let `z = 1(gate positive)`. The accepted CRIT logit correction can be expressed as `β0 + (β1−β0)z + γ0 r0 − γ0 z r0 + γ1 z r1`. A GLM with the same four effective predictors `(1, z, (1−z)r0, zr1)` exactly represents the same fitted tree. Novelty cannot be claimed solely from the 1-split functional form; candidate discovery, conditional selection, regularization, interpretability and statistical efficiency must be evaluated separately.
4. **Selection reuse.** The code evaluates as many as 6 candidates on the *same* heldout validation subset `V`, picks the best estimated lower gain bound and then uses it to accept/reject. This is susceptible to multiple-candidate selection optimism, despite D/V separation.
5. **Encoding mismatch is deliberate but should be studied.** Published RLR uses signed ternary rank gaps with `m=10`; CRIT gates/leaf predictors are binary direct expression-order comparisons, with ties falling in the `−1` branch. These are not identical biological features.
6. **Patient-level audit scope.** The Python auditor verifies identifiers, labels, predictions, metrics, fallback parity and leaf ID constraints. It does not alone establish absence of all selection-induced bias, infer statistical significance, or prove generalization.

## Recommended next gate (read-only, no model refitting)

- Extract the complete `run_manifest.csv` and split-level structures including accepted/rejected status, `confirmation_lower`, gate and leaf gene identifiers, actual class support, n_train/n_test, and per-split logloss/Brier for all four methods.
- Separately summarize per-dataset delta and compare methods **on accepted folds**, without presenting these retrospective conditional contrasts as inference.
- Review split-selection bias (6 comparisons reusing V), ties, stability and original paper parity against the uploaded main paper and Supplement.
- Before more model iteration, freeze a prospective full method comparison including **matched-feature GLM**, same acceptance policy controls, and an untouched external dataset/cohort/batch. The already-exposed 19 × 5 material may be descriptive development evidence only.
- Preserve fix3 SHA and files; create new code/experiment branches and output directories for any future scientific changes.

**Scientific gate: EXPERIMENTALLY FEASIBLE; SCIENTIFIC ADVANTAGE UNPROVEN.**

## Full six-fold read-only manifest inspection (operator output, 2026-10-25)

The following values were extracted by an operator-run Python script from the existing 24-row `run_manifest.csv`, verified to reside in the pinned fix3 worktree. No models were refitted. `confirmation_lower` is the selection-statistic lower value AFTER one-SE and complexity-cost adjustments, **not** an external confidence interval.

| Dataset | Fold | Tree accepted | Confirmation lower | Global interaction ΔLL | CRIT ΔLL | Original Brier | Global interaction Brier | CRIT Brier |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chen-2002 | 00 | 0 | −0.001232921 | −0.008467017 | 0 | 0.053951 | 0.049628 | 0.053951 |
| chen-2002 | 01 | 1 | +0.013741061 | −0.012521440 | −0.003734686 | 0.074066 | 0.069716 | 0.072860 |
| chen-2002 | 02 | 1 | +0.020378850 | +0.004055763 | −0.007803831 | 0.017754 | 0.019760 | 0.015370 |
| gordon-2002 | 00 | 0 | −0.001933791 | −0.000222638 | 0 | 0.000027 | 0.000019 | 0.000027 |
| gordon-2002 | 01 | 0 | −0.002417689 | +0.000679455 | 0 | 0.009584 | 0.009918 | 0.009584 |
| gordon-2002 | 02 | 0 | −0.001965667 | +0.000175792 | 0 | 0.028356 | 0.028277 | 0.028356 |

Train/test sizes: chen-2002 = 143/36 for each fold; gordon-2002 = 144/37, 145/36, 145/36. All 6 model-fit `warning_count` values reported as zero.

### Dataset-level (means of the three fold deltas in each dataset)

| Dataset | GLOBAL_CAL | GLOBAL_INTERACTION | CRIT_TREE |
| --- | ---: | ---: | ---: |
| chen-2002 | −0.000377215 | **−0.005644232** | −0.003846172 |
| gordon-2002 | −0.000125454 | +0.000210870 | **0.000000000** |

Across the six folds the arithmetic mean CRIT-minus-GLOBAL_INTERACTION delta is **+0.000793595** (smaller logloss is better). On the two accepted folds alone, CRIT has mean −0.005769259 vs GLOBAL_INTERACTION −0.004232839 relative to original RLR. This is an observational, *post-selection* contrast; it is not a controlled efficacy comparison.

### Biological rule identifiers in the accepted folds

- chen-2002 fold 01: gate `f_16745 > f_6257`; FALSE leaf `f_23794 > f_5763`; TRUE leaf `f_20000 > f_4055`.
- chen-2002 fold 02: gate `f_3506 > f_5763`; FALSE leaf `f_23794 > f_7943`; TRUE leaf `f_20682 > f_6297`.

Only the FALSE-leaf relation contains the same feature identifier `f_23794` on both accepted folds; the paired partner changes. No stable complete gate or leaf-pair recurs across those two fits. These feature IDs are NOT automatically gene symbols; any biological annotation requires a documented mapping from the locked source data.

### Decision for the next research stage

1. **Do not run the entire exposed 19×5 benchmark yet.** Freeze the six-fold evidence and diagnose effect sizes, routing, and selection optimism from saved predictions only.
2. Make any literal GLM using the exact `(1, z, (1−z)r_false, z r_true)` model design a **representation/parity check**; it must produce the same fitted predictions under identical penalty coding and optimization, not be advertised as an independent predictive baseline.
3. A genuinely informative comparator is a **non-tree interaction learner** with an aligned feature-selection budget, a comparable regularization budget, and the **same confirmation/reversion-to-original policy**. Preserve a separate unconfirmed interaction control as a diagnostic.
4. Fix the repeated-use-of-V design: nominate a *single* gate/candidate from discovery D, and use V only for its predetermined acceptance test, or pre-register a statistically defensible multiplicity correction. With n≈143–145 in outer training, subdividing V again may be unstable.
5. Define and lock independent external dataset/cohort/batch validation before future method selection. The original 19 used datasets cannot serve that confirmatory role.

**Next scientific gate: paired patient-level diagnostic audit and a locked fairness protocol; avoid premature CRIT-v8 tuning.**
