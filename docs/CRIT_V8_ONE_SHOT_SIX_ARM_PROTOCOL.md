# CRIT-v8 — Single D-Nomination, Fair One-Shot Controls (development protocol)

Date: 2026-10-25. Stage: **EXPOSED_RETROSPECTIVE_DEV_ONLY**.

## Why v8

The frozen CRIT-v7 six-fold pilot passed the complete patient-level audit, with accepted conditional trees on two chen-2002 folds and none on three gordon-2002 folds. Post-hoc Stage B showed:

- GLOBAL_INTERACTION: mean fold Δlogloss −0.002716681, 3 wins / 3 losses.
- CRIT_TREE: −0.001923086, 2 wins / 4 exact original-RLR fallbacks.
- GATE_ONLY: −0.001019299, 4 wins / 2 losses.
- ADDITIVE_PAIRS: +0.000439279, 4 wins / 2 losses.
- UNCONDITIONAL_TREE: −0.002270963, 5 wins / 1 loss.

All data/folds were previously inspected during RLMT development. Stage B indicated a meaningful conditional effect in chen fold 02 but a largely intercept/split effect in fold 01; its mean numbers should **not** be described as independent estimates of external accuracy. A single split plus two separate pair rules is algebraically representable by explicit GLM interaction features, so representational novelty alone is not proven.

CRIT-v7 examined up to 6 candidates on the same inner validation set V and chose the best; v8 removes that specific multiple-nomination step. It does not solve every selection bias or guarantee model superiority.

## Exact frozen dependencies

- Original published RLR checkout `4873878f57f806bf5b52144756f8087b17fa540c`, untouched.
- Protected v7 reference algorithm `e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03` (Git blob checks for eight files).
- One separate feature branch/worktree and fresh result output directory.
- No mutation of previous CRIT7 results, paper methods, or archived patient files.

## Predeclared v8 algorithm

Outer train T, outer test E heldout. For each T:

1. Fit EXACT original publication RLR (switchBox top1000 oriented relations, K=10, m=10, d=2, glmnet α=.5, lambda.min), and strict full-T 3-fold OOF original-RLR probabilities. The original model remains unchanged.
2. Split T once into stratified discovery D and confirmation V (seed 773, approximate 2:1). Refit original RLR parent + strict 3fold OOF **using D alone**. Confirm support and disjoint patient IDs.
3. Construct unsupervised variability/relation candidate pool on D, plus original **D-only** RLR genes and deterministic D-only proposals. Nominate **at most one** gate and separate leaf-pair relations with the original v7 residual interaction procedure, setting `shortlist=1`. Screening uses D outcomes and D-OOF probabilities exclusively. No information from V determines which candidate or leaf pair is nominated.
4. Build five correction designs from the *same* nominated gate/relations (and intercept-only calibration independent of a gate). Fit each correction with frozen ridge coefficients on D against D-OOF logit offsets.
5. On V, apply each D-fitted correction to `parent_D` predictions and calculate a single per-method heldout logloss-gain statistic: `mean_gain - 1*SE - 0.002`. Each method independently accepts if the statistic exceeds 0.002. **Do not choose a winning arm using V.** This is one V test per arm, although there are five simultaneously evaluated arms; multiple-arm interpretation still needs multiplicity caution.
6. Refit coefficients on ALL T using strictly T-OOF original RLR probabilities, preserving the selected gene pairs. Predict E with the original RLR model refitted on T, and apply only corrections accepted in step 5. Non-accepted or infeasible methods return **exact original RLR** (not merely an approximate baseline).
7. Independently audit every E patient against frozen original paper RLR predictions, metric recomputation, exact fallback, fold/sample identity and V-acceptance flags.

## Six arms; 6 outer folds × 6 = 36 method-fold prediction vectors

| Arm | Model correction | Acceptance |
|---|---|---|
| ORIGINAL_RLR | none; exact published baseline | no |
| GLOBAL_CAL | 1 offset intercept | its own one-shot V test |
| GATE_ONLY | two gate-leaf intercepts | its own one-shot V test |
| ADDITIVE_PAIRS | two leaf intercepts + 2 global ungated pair slopes | its own one-shot V test |
| GLOBAL_INTERACTION | gate, one nominated relation and gate×relation | its own one-shot V test |
| CONDITIONAL_TREE | gate-specific intercept and one separate pair slope per leaf | its own one-shot V test |

All use fixed v7 ridge intercept=.10 and slope=.20 for corrections, same cap ±2.5 and one-SE/cost 0.002 decision threshold. Comparators share the **D-selected** candidate rather than tuning a separate search for their own benefit, which is appropriate for *mechanistic controls* but does not constitute exhaustive fair optimization of each model family.

The earlier v7 model can be compared descriptively to v8, but this is **not** a clean validation because v8 was designed after seeing these exact data and v7 results.

## Synthetic and patient audit gates

- Original RLR model identity and 32-gene synthetic (m=10 is feasible), six responses with exact fallback;
- D/V disjoint, max one nominee, five predeclared per-arm acceptance decisions, no validation-based candidate choice;
- fitted-coefficient and accepted flag consistency, affine rank invariance, bit-exact RDS replay, no-feasible-candidate fallback;
- identical original outer CV sample IDs, labels, SHA256 input files and archived baseline probabilities;
- 2 datasets (`chen-2002`, `gordon-2002`) × 3 original folds (`rep00_fold00`–`02`) × 6 methods; patient-level logloss and Brier independently replayed;
- no overwrite, separate SHA-pinned branch/worktree/results, fail-closed nonzero status.

Completion: `CRIT8_SYNTHETIC_ALL_GATES=PASS`, `CRIT8_36_METHOD_FOLDS=PASS`, `CRIT8_EXACT_ORIGINAL_RLR_PARITY=PASS`, `CRIT8_ONE_SHOT_36_VECTORS_AUDITED=PASS`, `CRIT8_PILOT_BATCH_EXIT=0`.

## Scientific interpretation and planned next gate

- Judge evidence by dataset and outer-fold consistency, accepted-arm rates, conditional versus gate-only separation, and per-patient robustness, not simply the lowest 6-fold mean.
- Every feature was nominated from D, but **the five arm-specific V decisions have not been multiplicity-corrected**. Do not select and advertise the best arm from V as if it were confirmed. If a single final algorithm is chosen, selection must be locked before testing an independent external cohort.
- One split still lies in the model class of a global GLM with interaction terms. Claim innovations in candidate discovery/selection, stability, interpretability or generalization only if independently demonstrated.
- The same published 19 cohorts or the same original split patients are **exposed development material**. An untouched independent dataset/cohort/batch is required for a novel algorithm performance claim.

**Stage v8 is an engineering/methodology hypothesis test, not a new publication-ready superiority claim.**
