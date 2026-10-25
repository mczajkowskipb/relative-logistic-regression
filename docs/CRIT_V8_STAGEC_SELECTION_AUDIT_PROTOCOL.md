# CRIT-v8 Stage C — archived v7/v8 candidate-selection audit

Date: 2026-10-25. **EXPOSED_RETROSPECTIVE_DEV_ONLY.** Code and original results are frozen. Stage C is a diagnostic protocol, not a new prospective outcome evaluation.

## Motivation

The user-operated Ubuntu pilot finished with `CRIT8_PILOT_BATCH_EXIT=0`, `CRIT8_SYNTHETIC_ALL_GATES=PASS`, `CRIT8_ONE_SHOT_36_VECTORS_AUDITED=PASS`. All five corrections were independently rejected in all six original outer folds. `nomination_count=1` in all six folds, and all 30 arm-level conservative V statistics were negative. This is technically correct fail-closed behavior but no observed gain in predictive performance relative to exact original RLR.

The v7 pilot previously evaluated up to 6 D-discovered candidates on the same V and accepted the best on two chen-2002 folds (V lower bound +0.013741061 and +0.020378850). v8 nominated only D rank one (same original D/V support and nomination rules, except shortlist=1) and rejected it on these two folds (tree conservative V score −0.054107015 and −0.018375153).

**Working hypothesis, to be tested rather than assumed:** the v8 nominee is exactly v7's first D-ranked candidate, and differences arise from v7 searching for the maximum V-scored candidate. This is a causal claim about the *algorithmic selection operation*, not proof of statistical overfitting or false discoveries.

## Stage C scope and strict acceptance

Read saved `complete_crit7.rds` and `complete_crit8.rds` for 2 datasets × 3 frozen folds. **No fitting, no bootstrap, no new RLR training, no use of original outer-test outcomes for model choices.**

- Validate inner D/V patient IDs and lengths are identical in v7 and v8 and do not overlap; verify exact original model identities.
- Extract every saved v7 D-ranked candidate and its V lower-bound, raw mean gain, standard error, nominated gate and two leaf relations. On each fold independently recompute that v7's V-selected candidate equals the maximum recorded candidate score.
- Assert that v8's *single D candidate* is identical to the **first** v7 D-ranked gate, leaf relations and nominated interaction. Assert its V conditional-tree **mean, SE, and conservative lower** are numerically identical to v7 candidate 1 (tolerance <1e−8).
- Examine the **rank** of v7's chosen candidate under D nomination; count how often V changed the identity from D rank 1. Compare v7/v8 V accept flags to outer-test delta, described *only post hoc*.
- Extract all 30 v8 method-wise acceptance decisions including GLOBAL_CAL. Independently verify `lower = validation_mean − validation_se − 0.002` and `accepted = lower > 0.002` for each; the original code's formal policy is fixed. The lower value is a heuristic conservative gain statistic, **not a confidence interval**.
- Verify the six v7/v8 outer performance records match the frozen manifests. No changes to external data or results.

## Output

Separate branch `analysis/crit-v8-stageC-selection-audit-v1`, new worktree and new result directory `results/crit_v8_stageC_selection_v1/`:

- `selection_folds.csv` — six-fold D/V size, v7 candidate count, v7 best D rank vs v8 first D rank and validation scores, heldout summary
- `v7_candidate_rankings.csv` — all v7 candidates with D-ranking and V-confirmation metrics
- `v8_confirmation_arms.csv` — 30 V decisions, raw mean / SE / conservative statistic and outcome
- `stageC_R.log`, `independent_audit.log`, `status.txt`, locked SHA files.

## Interpretation gates

- If D first candidate and v7 candidate 1 match exactly, then the major v7/v8 difference is **which candidate is accepted**, not changing the RLR baseline or D/V split. When the V-selected candidate has D rank >1, the v7 advantage is plausibly driven by V-based selection. **It does not prove selection bias**; only independent untouched test cohorts and/or controlled simulations can quantify that.
- If V-selected candidate is D rank 1 but v7/v8 decisions diverge, this signals a bug or hidden protocol difference: stop and audit immediately.
- Negative V scores even for GLOBAL_CAL suggest either a strong RLR baseline, instability between D and V, overly stringent acceptance or regression model mismatch; these are *competing hypotheses*, not evidence of a specific flaw.
- Only consider a new candidate-ranking objective or statistical confirmation policy as a **new prospectively predeclared method**, not tuned on these six previously exposed outer folds.

**Do not label any later experiments on the same six folds or 19 development cohorts independent validation.**
