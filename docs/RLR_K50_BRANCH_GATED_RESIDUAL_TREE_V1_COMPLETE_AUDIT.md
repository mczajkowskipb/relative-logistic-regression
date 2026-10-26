# RLR K50 rank-gated residual tree V1 — complete scientific gate, 2026-10-26

## Task

Develop a **genuine interpretable tree extension of the original, unpublished RLR manuscript**, limited to the eleven binary CompCancer datasets already present in that manuscript, so the eventual single extended article could fit a habilitation cycle on interpretable tree induction for omics. Preserve the original RLR core and locked paper. Do NOT claim novelty or superiority unless supported by externally evaluated data.

## Exact predeclared mechanism

- The root is the **literal original RLR K=50** (m=10 ternary rank evidence, d=2 strict gene disjointness, original supervised switchBox top1000 nominations, cv.glmnet alpha=0.5/lambda.min, unchanged `R/new/07_world_protocol_utils.R`). It supplies an original K50 logistic *offset*; K10 does NOT replace it.
- Candidate ternary roots and branch-specific sets of 3 disjoint gene-pair rank relations are nominated on inner D using a pool of original top1000 plus a deterministic expansion from **the same original dataset's full gene universe**, no other cohorts.
- A single globally regularized `glmnet` correction, `alpha=0.5`, `lambda=0.05`, `intercept=FALSE`, adds only branch-gated relational effects to the fixed original RLR logit. **No independent unstable per-leaf regression.**
- Root and conditional pair sets selected using **D only**, then candidate corrections must satisfy the inner V absolute logloss improvement **≥0.01** and BA change **≥−0.02**, with a penalty for correction nonzero coefficients. On failure, return the **byte-identical original RLR K50 heldout E probabilities**.
- This is a **hierarchical branch-gated logistic interaction classifier**, not a conventional terminal-class tree. Such interaction models have substantial prior art; originality and true decision-tree qualification are not established by nomenclature.

## Scope and reproducibility

Only paper-inventory `chen-2002` (179 samples), `gordon-2002` (181), `yeoh-2002-v1` (248). These cohorts have been exposed in previous exploratory tree studies and ARE NOT a fresh independent confirmatory sample. Public Schliep source text matrices are reconstructed, original processed ARFF SHA parity **NOT VERIFIED**.

Three disjoint stratified outer folds, original RLR K10/K50, original TSP1 and kTSP9, and a two-level greedy TSPDT **proxy** (NOT the author's original CUDA TSPDT). The comparison has 3×3×6=54 method×fold evaluations and **3648 heldout patient×method predictions**, independently audited for BA, macro-F1, logloss and Brier on exactly matched patient IDs.

[GitHub Actions #40381422977 — all 3/3 SUCCESS, complete independent score/partition audit](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40381422977). Earlier run `40380452256` failed because the Python AUDITOR inherited a 1e−7 logloss probability clipping threshold while RGT50 uses 1e−8. The audit was fixed to match the model, without changing its fitting/training/evaluation. Later diagnostic changes only count valid training candidates, not model predictions.

## Scientific results

| Original-paper dataset | RLR K10 BA / LL | RLR K50 BA / LL | kTSP9 BA / LL | RGT50 BA / LL | Accepted genuinely tree-augmented outer folds |
|---|---|---|---|---|---:|
| Chen | 0.91863 / 0.18255 | 0.91291 / 0.19818 | 0.87331 / 0.32202 | 0.91291 / 0.19818 | 0/3 |
| Gordon | 0.96485 / 0.04878 | 0.98152 / 0.03275 | 0.98152 / 0.05194 | 0.98152 / 0.03275 | 0/3 |
| Yeoh | 0.93316 / 0.10098 | 0.92446 / 0.10106 | 0.86829 / 0.45874 | 0.92446 / 0.10106 | 0/3 |

RGT50 exactly returns original K50 on all 608 heldout people pooled across three complete 3-fold datasets. **BA/LL parity is BY FALLBACK, not empirical success of a new tree.** The stronger `K50` hypothesis itself is not universally the strongest baseline; e.g. Chen and Yeoh original K10 have slightly better mean BA/LL. Do not claim K50 is always better.

### New decisive train-only causal/algorithmic diagnosis (unchanged search and lambda)

Across 9 inner D/V folds, root proposals count was **8 per fold**, 72 in total.

- Chen: 8 valid triple-relation branch designs AND 8 successfully fitted penalized logistic corrections in EACH of three folds → **24 valid candidate fits**. Nonzero added elastic-net coefficients appeared in **1 of 24** only. **0/24** gave even strictly positive V logloss improvement; recorded maximum was 0 within displayed rounding. Therefore no accepted tree in Chen.
- Gordon: zero valid 3-relation branch designs over 3 outer folds (0 candidates fitted). Pair/root availability alone did not guarantee enough branch-local sample and class support under the fixed rules. This is an **algorithmic feasibility bottleneck**, not a negative E outcome from fitted trees.
- Yeoh: likewise zero eligible three-relation branch designs over 3 outer folds. This exposes the strong limitations of stringent local support/three-relation disjointness and residual feature gating. The algorithm never learned a real conditional tree on this dataset.
- Thus, 9/9 exact K50 fallbacks. The protection against severe loss works; **the feature/branch selection mechanism has not been validated**.

## Scientific interpretation and stop rule

**No publication-ready drzewiasta extension yet.** Do not characterize these three dataset results as direct evidence that trees cannot improve RLR; only this prespecified branch-gated penalty/feature availability construction failed to produce a model passing its inner training/validation gate.

The repeated failures have at least TWO distinct mechanisms: (a) no eligible branch-wise 3-relation data support on Gordon/Yeoh, and (b) successful fit but overwhelmingly zero additional coefficients / no inner validation LL gain on Chen. Without an entirely redefined and justified inductive bias, further E-driven parameter tweaking (lambda, branch-size thresholds, roots, K) would overfit the already exposed datasets.

The math-proven certified RLR early exit C-RET remains the strongest robust interpretation improvement: exact original RLR class, bounded posterior approximation, 25–43% fewer evaluated pair-relations in prior paper cohort pilots. But it remains largely a sequential/anytime evaluation representation; learned ternary policy/global two-step optimization have not established additional tree-specific utility. Neither that C-RET result nor these RGT50 results prove predictive superiority over the **original CUDA TSPDT**, because prior comparisons used a GREEDY PROXY only.

**Next logically defensible work:** stop this specific search branch, audit novelty and data provenance, and distinguish the habilitation contribution of already accepted evolution/multitest tree work from incremental unpublished RLR. A genuinely new revised TREE algorithm could be considered only with predeclared hypothesis and confirmatory evaluation on original processed ARFF data; do not claim success from these negative pilots. The original RLR manuscript and its main source remain unchanged.
