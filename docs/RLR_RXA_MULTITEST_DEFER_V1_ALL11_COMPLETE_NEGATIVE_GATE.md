# RXA MultiTest rank-margin deferral cascade V1 — complete negative science gate

**2026-10-31 — TECHNICAL PASS; SCIENTIFIC NEGATIVE.** No changes to original RLR source, or to the CMPB manuscript. This exploratory design was frozen before run outcomes, but 11 CompCancer cohorts were repeatedly exposed in earlier research.

Real completed **12/12 GitHub Actions PASS**, 11 cohorts × 5 E folds × 6 methods = 330 method-fold evaluations with independent every-patient prediction and routing audit. GitHub Actions run: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41232466179 . Fitted code SHA: ab38971d07ec078e107582d25b87412783ee0f40.

**Provenance:** public Schliep source-reconstructed datasets, original publication processed ARFF SHA **NOT VERIFIED FOR THIS RUN**. The separate author's original processed ARFF verification from the earlier CompCancer study is NOT being attributed to the present GitHub Actions run.

## Design

- Each split uses the exact original outer T/E 5-fold partition and original RLR20 K20/d2/top1000 switchBox relation pool. Same six first gene-disjoint pair identities across all multi-test arms.
- Two sequential RXA nodes, each a 3-pair multi-test. A node is confident only with >=2 oriented votes for the same class and no contradictory vote; otherwise it defers all patient observations to the next node.
- Binary state m1 versus rank-margin ternary state m10 (+1,0,-1), same pair orientation learned on outer T only. No supervised E feature nomination.
- Two terminal modes: smoothed train-only class-proportion leaves (weight=10, minimum support=5) versus **global RLR20 logistic rescue ONLY for unresolved patients**, fully trained on outer T. This is not local leaf-specific glmnet.
- Crucial controls: flat original RLR20, and flat penalized elastic net using exactly the same six rank-margin pair columns as the node candidates.
- No patient duplication, no separate branch-learned RLR fits, no teacher-surrogate distillation, no GPU original EMTTree/RMCT performance claim.

## Audited cohort-equal 11-dataset means

| Model | Logloss | Brier | AUC | Positive F1 | ECE | Pair comparisons per patient | Logistic rescue rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| Original RLR20 | 0.237944 | 0.062638 | 0.961673 | 0.876180 | 0.092092 | 18.9818 | 0% |
| Binary multi-test two nodes, ordinary leaves | 0.343810 | 0.103441 | 0.908488 | 0.810722 | 0.133577 | 4.0456 | 0% |
| Ternary multi-test two nodes, ordinary leaves | 0.341689 | 0.102616 | 0.910289 | 0.821349 | 0.128666 | 4.0026 | 0% |
| Binary multi-test, original RLR fallback | 0.319945 | 0.090980 | 0.916525 | 0.836579 | 0.147287 | 7.9535 | 21.10% |
| Ternary multi-test, original RLR fallback | 0.315619 | 0.089791 | 0.921959 | 0.851601 | 0.143661 | 7.7365 | 20.25% |
| No-tree six-pair ternary elastic-net logistic | 0.313667 | 0.073401 | 0.950279 | 0.842862 | 0.107443 | 6.0000 | 0% |

Ternary cascade plus RLR rescue: 66.58% of patient cases exit at first multi-test, 13.17% exit second, and 20.25% are deferred to global RLR. Mean number of evaluated pair indicators is 7.74 versus 18.98 for full RLR (59.2% fewer evaluations), **but the full original feature universe remains required for rank computation; this does not mean fewer measured genes**.

Probability quality is substantially WORSE than full RLR: mean logloss +0.07768 (~32.6% worse), positive F1 −0.02458. The flat six-pair logistic control has lower logloss, lower Brier, higher AUC and fewer evaluated pair relations than the ternary+RLR cascade; its F1 is 0.00874 lower, but this is not preregistered evidence of clinical superiority.

## Fixed paired dataset-level exact 2^11 sign flips; exploratory BH across four tests

| Ternary minus binary | Endpoint | Mean delta | Raw exact p | BH q |
|---|---|---:|---:|---:|
| Ordinary leaves | Logloss | −0.002121 | 0.531250 | 0.708333 |
| Ordinary leaves | F1 | +0.010627 | 0.750000 | 0.750000 |
| Global RLR fallback | Logloss | −0.004326 | 0.437500 | 0.708333 |
| Global RLR fallback | F1 | +0.015022 | 0.375000 | 0.708333 |

**All differences non-significant.** This is a genuinely negative science result, not equivalence proof.

## Decision

**STOP V1**; do not adjust margin, node width/number, vote conditions, support/smoothing, pair search or child-logistic weighting after reading these E results and then portray the rerun as prospective. This preliminary cascade is NOT a full new DT induction or the historic full CUDA RMCT/EMTTree. Earlier model-tree and teacher-guided depth2 explorations were also scientific negative. It is not defensible to add a new tree superiority claim to the almost-submitted RLR CMPB manuscript or assert a new habilitation-ready contribution from this experiment.

The current CMPB manuscript should be submitted on its existing RLR results. Further genuine tree claims would require different independent biological tasks or prospective frozen evaluation, appropriately matched real RXA tree-induction comparators.

## Auditable references

Frozen protocol: docs/RLR_RXA_MULTITEST_DEFER_V1_FROZEN_PROTOCOL.md
Actual six-arm implementation: R/experimental/rlr_rxa_multitest_deferral_v1.R
Independent Python audit: scripts/audit_rlr_rxa_multitest_deferral_v1.py
GitHub workflow: .github/workflows/rlr-rxa-multitest-deferral-v1.yml
The run artifact includes all 330 fold-method rows, 55 train/test manifests, all per-patient scores/routes, per-dataset scores, exact paired comparisons and the full independent audit log.
