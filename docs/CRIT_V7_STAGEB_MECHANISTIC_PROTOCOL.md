# CRIT-v7 Stage B — mechanistic ablation and exact global representation, 2026-10-25

## Why this is the next scientific gate

The independently audited fix3 six-fold original-RLR-aligned pilot passed all technical gates and accepted a two-leaf tree in 2/6 folds, both within **chen-2002**. The other four folds fell back exactly to RLR. Aggregated mean fold Δlogloss vs original RLR: GLOBAL_INTERACTION −0.002716681; CRIT_TREE −0.001923086. Both datasets and original folds are repeatedly exposed to development. **No external generalization or originality claim is permitted.**

A subsequent read-only operator diagnosis found 27/36 patient-level logloss improvements for CRIT in each accepted chen fold, with gains mainly in leaf 2:

| Fold | Leaf 1: n, positive, mean ΔLL | Leaf 2: n, positive, mean ΔLL |
|---|---|---|
| chen fold 01 | 20, 12, +0.001848557 | 16, 3, −0.010713739 |
| chen fold 02 | 21, 1, +0.000471172 | 15, 14, −0.019388835 |

For fold 01, deleting the patient with largest improvement eliminates the net heldout logloss gain; fold 02 is more robust to that one-case deletion. This is a **mechanistic diagnostic**, not a formal inference.

## Frozen boundaries

- Exact original RLR archival reference SHA: `4873878f57f806bf5b52144756f8087b17fa540c`.
- Prior completed CRIT-v7 fix3: `e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03`.
- Both prior worktrees and their results remain untouched.
- New branch: `analysis/crit-v7-stageB-mechanistic-v1`; one new git worktree and one fresh output root.
- The eight protected RLR/CRIT source, runner, and audit dependencies are verified by Git blob hash against fix3.
- No search over more gates and no modification to original CRIT thresholds. All stage B models reuse the **already selected fix3 gate and leaf relations**; those were selected with the earlier multi-candidate use of V and hence inherit its development-selection limitation.

## Predeclared comparison: same outer test patients, same gate/pair features

6 original outer folds × 6 methods = 36 saved patient-level prediction vectors:

1. `ORIGINAL_RLR`: exact published baseline, copied from the saved fix3 model; parity checked against original archive.
2. `GLOBAL_INTERACTION`: frozen fix3 single screening interaction (no fallback); copied, not retrained.
3. `CRIT_TREE`: frozen fix3 accepted/pruned tree; copied, not retrained.
4. `GATE_ONLY`: **two leaf intercepts only**, using the same selected split (no leaf-pair slopes).
5. `ADDITIVE_PAIRS`: two leaf intercepts plus two un-gated global pair slopes; no gate×pair product.
6. `UNCONDITIONAL_TREE`: same 4-parameter leaf model fitted on full outer train **without** the prior V-based acceptance requirement; diagnostic only, not a recommended pipeline.

All **new** coefficients use exact original-RLR 3-fold outer-training OOF probabilities as offsets, with the same `crit7_fit_offset` ridge penalties as the frozen tree. The root fitted probabilities are obtained from the archived complete CRIT model. The outer-test outcomes are used *only* to compute metrics and audit; never for fitting or selection.

**Important distinctions:**
- GATE_ONLY changes only the correction's structure, isolating the effect of the subgroup split.
- ADDITIVE_PAIRS keeps two relational predictors but removes their conditional activation.
- UNCONDITIONAL_TREE tests whether the conservative accept/reject policy sacrifices gains or protects against spurious corrections; it is not a prospectively selected method.
- The exact tree equals the **global** GLM feature expansion `(1, z, r_false, z*r_false, z*r_true)` with transformed coefficients. Predictive parity must be within 1e-12; this equivalence test does NOT establish any original model class.
- Because feature selection is inherited from fix3, this is descriptive *post-hoc mechanism ablation*, not fair independent evidence of one learning algorithm's superiority over another.

## Outputs and independent replay

- `results/crit_v7_stageB_mechanistic_v1/CODE_SHA.txt`
- `results/crit_v7_stageB_mechanistic_v1/FIX3_SOURCE_SHA.txt`
- `results/crit_v7_stageB_mechanistic_v1/PROTOCOL.txt`
- `results/crit_v7_stageB_mechanistic_v1/comparison_manifest.csv`
- `results/crit_v7_stageB_mechanistic_v1/patient_predictions.csv`
- `results/crit_v7_stageB_mechanistic_v1/stageB_R.log`
- `results/crit_v7_stageB_mechanistic_v1/stageB_independent_audit.log`
- `results/crit_v7_stageB_mechanistic_v1/status.txt`

The Python auditor independently recomputes logloss and Brier from the stored patient predictions, checks all original and fix3 identities and fixed counts, tests exact rejected-tree fallback and summarizes leave-one-patient-out robustness in the two accepted folds.

## Decision gate after execution

Assess **GATE_ONLY vs CRIT**, **ADDITIVE_PAIRS vs CRIT**, and **UNCONDITIONAL_TREE vs confirmed CRIT**, separately on chen/gordon. Also check the distribution of per-patient gains within each leaf. If GATE_ONLY explains the change, focus on robust relational subgroup detection rather than leaf-local models. If interaction terms add reproducible value, formalize a **new** method with a single D-only candidate nomination and one predeclared V acceptance test, plus equally gated comparator controls. Do not tune on these two exposed datasets, and do not describe further exposed 19×5 tests as fresh confirmation.

**Technical completion criteria:** `STAGEB_36_VECTORS_INDEPENDENTLY_AUDITED=PASS`, `STAGEB_EXISTING_METHOD_PATIENT_PARITY=PASS`, `STAGEB_GLOBAL_GLM_EXACT_REPRESENTATION=PASS`, `STAGEB_BATCH_EXIT=0`. These have not been observed on Ubuntu as of commit creation.
