# Selective Anchor-Gated RLR (SA-RLR) V1 — implementation and verification status

**Date:** 2026-10-27. **Mechanical implementation COMPLETE; biological efficacy NOT YET MEASURED.**

- Isolated source branch: `analysis/rlr-selective-anchor-gated-v1`.
- Code: `R/experimental/rlr_selective_anchor_gated_v1.R`
- Executable deterministic synthetic test: `R/experimental/test_rlr_selective_anchor_gated_v1.R`
- **Real GitHub Actions mechanical test [#40474192037](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40474192037) SUCCESS**; R fit, original manuscript presence check, reference-gene exchange/score invariance, and exact inactive-branch fallback PASS.
- Original RLR implementation and original locked unpublished article were not edited.

## Motivating defect in prior model tree

Previous D-only out-of-fold gradient lookahead fit working *two* RLR child experts in only 2 of 6 exposed Chowdary/Singh outer partitions. That unnecessarily discards a valid branch that contains residual biological signal when the other branch is nearly homogeneous, has poor class support, or provides no gain. The user also proposed **optimizing the reference genes by their downstream children's quality**, not merely reranking global TSP-optimal anchor panels.

## Proposed model and ordinal invariance

Let original RLR-50 predict `p0(x)` on fixed-gene-universe patient x. A single decision-stump root uses ordinal comparisons:

`S(x;g,A)=sum_{a in A} 1[x_g>x_a]`, with `|A| in {3,5}`.

`b(x) = 1[S(x;g,A)>=t]`, for `t=1..|A|`.

The conditional extension is

`logit p(x)=logit p0(x)+1[b(x)=0]Delta0(x)+1[b(x)=1]Delta1(x)`

Each `Delta_b` is an elastic-net-regularized (alpha .5, fixed lambda .10) additive signed ternary m10 rank-relation model fitted on at most 8 train-selected supervised branch pairs. The original RLR is still the fixed, shared parent logit offset. **Either Delta0 or Delta1 may be exactly zero and then its patients receive exact original RLR probabilities.** Prediction includes no raw expression values or ratios; an independently chosen strictly increasing transform of each patient's measurement vector leaves the classification invariant given the same gene/probe universe and unchanged tie ordering.

## Bounded actual anchor-panel optimization

1. On inner D only, nominate <=8 candidate target genes from top supervised original RLR/TSP pairs.
2. For each, form baseline reference panels of size 3 or 5 with same D-only rank-stability/TSP-neighbor heuristic.
3. **Actively swap one reference anchor** at each panel position with at most four train-only alternatives (nearest candidate neighbors or rank-stable genes). Enumerate ordinal thresholds t=1..k and reject duplicate patient partition patterns. Cap at **240 candidate root partitions** in total.
4. Train the original RLR50 on D and **crossfit D-only RLR50** predictions from three internal stratified folds; use crossfitted residual `e=y-p_OOF`.
5. Cheap screen all roots by second-order approximate gain available from branch-specific ternary relation corrections over **the same first 180 supervised pair feature columns**, subtracting a proportionate global nonbranch correction score and penalizing larger reference panels. Importantly, the screen is a proposal score, not an honest estimator of validation performance; **it evaluates future branch-specific effects, not just the root's direct class separation**.
6. Retain <=6 root finalists. On D, refit *real* branch-specific supervised local pair selection and conditional RLR models for each candidate. **Require at least one**, not two, nonzero child conditional expert. Apply those frozen candidate structures only on D-independent inner V and choose min V logloss among eligible finalists.
7. Inner V acceptance gate: selected conditional model must improve V logloss by >=.005 relative to original RLR50, by >=.002 relative to one global non-tree model using **exactly the same selected child pair union and root gate**, and worsen BA by no more than .02. Otherwise exactly return RLR50 on outer E.
8. Refit parent RLR50 on all outer training T, keep D-selected anchor target/reference IDs, root threshold and local pair identities frozen, refit only eligible local coefficients on T. Predict heldout E once. E outcomes never enter candidate selection or threshold choices.

The model is intentionally a compact **rank-gated relational model tree**, not a multi-level EA-DT. If this novel selective-leaf mechanism is ineffective, no further EVOLUTIONARY optimization is justified on already inspected E folds.

## Mechanical tests: all actual CI PASS

- One branch has a fitted sparse relational logistic expert while the other uses **exactly** original RLR parent predictions.
- Prediction is invariant to **independent monotone exponential transform of each of 180 patient profiles**.
- Ordinal target/anchor reference one-swap candidate enumeration and conditional OOF gradient-proxy scores are invariant under the same transforms.
- Constant gene schema remains an explicit requirement; reordered input gene universe is rejected.
- Real glmnet conditional RLR fit and prediction succeed on planted synthetic branch-conditional structure.
- These are software assertions, **NOT real biological predictive-performance or novel-method scientific validation**.

## Required next scientific check — no additional E fishing

1. Mechanistic D-only diagnostic on original-paper cohorts: how many proposed root candidates have one versus two viable local child experts, how often swaps change actual partitions, number of unique conditional pairs, stability under training bootstrap, and performance on inner V.
2. **Freeze SA-RLR code and thresholds** before any more heldout E evaluation, and use the **original full 11 datasets / original 470 patient partitions**. The user's original processed ARFF checksums must be verified on Ubuntu server first if we claim an exact manuscript reproduction. GitHub's downloaded public-source reconstructions do not have identical ARFF checksums.
3. Comparison must include original RLR10/20/50, actual SWAP.KTSP 10/20/50, LR train-only MAD, gate-only, global same-pair (not tree), greedy-root child RLR50, SA-RLR, class-threshold control fitted on inner training only and ideally actual original CUDA TSPDT and 2024 ranktreeEnsemble (Lu et al.). The last latter models are external/literature contemporaries, not inherently equivalent to original TSPDT.
4. Primary target: **discrimination** (manuscript positive-class F1 and AUC, plus Macro-F1 and BA), **non-inferior logloss** vs RLR50, and explanatory compression (unique genes, number of active pairs on one patient path). Actual improvement over original RLR is welcome but not necessary if SA-RLR adds *distinct tree interpretability* while RLR and SA-RLR as a family beat corresponding kTSP on predeclared metrics.
5. Analyze dataset-level paired differences across 11 study units, not 470 correlated folds. Report uncertainty, not just arithmetic wins; avoid retuning on the repeatedly exposed Chowdary/Singh E.

## Existing paper-comparison factual anchors

Original **locked manuscript** CompCancer aggregated (11 datasets): kTSP50 AUC .963 / F1 .884 / logloss .378; RLR50 AUC .957 / F1 .869 / logloss .226. Thus the scientific gap is primarily **F1 and AUC**, not another huge reduction in probability loss.

Previous **two exposed cohort** diagnostic [real benchmark run #40473159026](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40473159026) compared real switchBox kTSP, RLR, value LR and greedy model trees on exactly matched six E folds, 2060 independently audited patient-method scores: RLR50 Macro-F1 .9400 / logloss .25337; greedy-root RLR50 model-tree Macro-F1 .9400 / logloss .24705; kTSP50 Macro-F1 .9201 / logloss .28129; value LR Macro-F1 .9556 / logloss .33387. The tree improved loss but **not** Macro-F1 or AUC over RLR50. None of these pilot metrics may be substituted for the 11-dataset manuscript metrics.

## Novelty and publication warning

Anchor groups themselves are covered in TSPG and related rank-based methods. Decision trees over relative gene-expression orderings are prior art, including the user's original TSPDT and the 2024 **ranktreeEnsemble** ensemble publication. Model-based recursive partitioning and shrinkage are also prior art. The candidate distinctive contribution is the **explicit one-swap anchor optimization driven by projected downstream conditional RLR information, with asymmetric selective sparse relational experts and exact parent fallback**, supported by rigorous full-cohort controlled evidence. Whether this is genuinely novel and useful remains to be tested and reviewed.

No claim yet that SA-RLR beats kTSP, LR, original RLR or ranktreeEnsemble, nor that it constitutes a publication-ready method. Original finished RLR article is the safe alternative if full 11/470 evidence fails.
