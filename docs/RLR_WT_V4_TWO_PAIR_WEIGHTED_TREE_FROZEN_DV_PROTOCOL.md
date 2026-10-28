# RLR-WT v4 — Two-pair ordinal logistic evidence nodes (frozen exploratory protocol)

**Date:** 2026-10-28. **Original RLR manuscript is unpublished.** Lock `manuscript_lock/RLR3_submission_final.tex`, original `R/new/07_world_protocol_utils.R`, and all previous original 470-fold benchmarks. This branch is isolated; no original method alteration.

## Testable hypothesis

The previous one-pair trees `RLR-DT v1/v2/v3` were genuinely small (usually 1–2 comparisons), but cross-fitted RLR probabilities did not confer stable F1 advantage over label-only trees. v3 first original E fold 11-dataset means label-only F1 .794970 vs RLR-guided .791507; teacher-guided AUC .876619 vs label-only .862528, unadjusted paired p=.50, logloss .359315 vs .396890 p=.125. All earlier E datasets already exposed, not independent.

**New parameterization:** consider a binary tree node operating on TWO ternary relative gene pair rank indicators rather than one:

[
r_{ab}(x)=
\begin{cases}
+1 & rank_x(b)-rank_x(a)\ge 10,\\
-1 & rank_x(b)-rank_x(a)\le -10,\\
0 & \text{otherwise}.
\end{cases}
\qquad
N(x)={\bf1}\{w_1r_{a_1b_1}(x)+w_2r_{a_2b_2}(x)>0\}.
]

The weights `w1,w2` are the signed coefficients of the **original RLR20**, learned from D and normalized by max absolute β; no new fit in a small leaf. Candidate features are the 28 unordered combinations of the **eight highest-|β| original RLR20 pairs** (after original switchBox selection, disjointness `d=2`, original `m=10`, `glmnet alpha=.5, lambda.min`). A `SIGN` ablation replaces every β by `sign(β)`. Equal contribution on opposite ternary relations becomes a tie (=negative branch), whereas nontrivial magnitudes can resolve it; the separate mathematical unit test proves the distinct behavior.

Greedy best-first induction uses train D soft cross-entropy; fixed maximum **two levels, three internal nodes, 4 leaves, 3 distinct two-pair predicates**, hence at most **4 evaluated gene-pair comparisons** per actual patient path, up to six distinct pair relations and 12 genes in the entire model. **At inference, do not compute all 8 original candidate relations**; compute full assay ranks once and check only the two relational features of visited nodes. The actual *assay probe footprint* is the whole measured rank universe and is distinct from the 2–4 path tests.

Two teacher factors define **four frozen arms**, all on exactly the same D data:

| Arm | Coefficients in node | Teacher probability in split objective | Teacher probability in leaf |
|---|---|---|---|
| WT_BETA_LABEL | true original `β` magnitude/sign | No, α=0 | No |
| WT_BETA_OOF | identical original `β` | Yes, 25% true threefold out-of-fold `p_RLR20` | No |
| WT_SIGN_LABEL | `sign(β)` | No, α=0 | No |
| WT_SIGN_OOF | `sign(β)` | Yes, α=.25 | No |

Leaves use smoothed **real-label** class proportions in all four arms. The α=.25 cross-fitted RLR teacher improves only split-score target, not in-sample labels. `WT_BETA_OOF − WT_BETA_LABEL` isolates probabilistic teacher influence *given same weighted candidate predicates*. `WT_BETA_LABEL − WT_SIGN_LABEL` isolates learned coefficient magnitudes. **Both β and SIGN four arms still depend on RLR preselection/orientation**, so only the separate TSP-label binary tree is the *RLR-independent* control.

Other controls on same V patients: `LABEL_BINARY_TSP_TOP20` (top20 switchBox pairs, same depth/max-internals/cost penalty; TSPDT-like proxy only), original `RLR20_TEACHER` fit on D, **real switchBox native** `KTSP5_NATIVE`. Native kTSP labels are not overwritten with a 0.5 threshold.

## Study scope and reproducibility

- [Synthetic unit workflow](https://github.com/mczajkowskipb/relative-logistic-regression/actions/workflows/rlr-weighted-pairpair-v4-unit.yml): explicit two-conflicting-relations β-vs-sign differential split, patientwise strictly increasing nonlinear transform invariance, max 4 actual comparisons, model standalone teacher-free; must show real CI success.
- [11-cohort D/V pilot](https://github.com/mczajkowskipb/relative-logistic-regression/actions/workflows/rlr-weighted-pairpair11-DV-v4.yml): original manuscript 11 Schliep public named datasets, original first outer T only, **inner D/V**, no outer E read for this study. Original processed paper ARFF SHA parity remains **NOT VERIFIED**. Original `rg_parent` parameters fixed; RLR OOF 3fold teacher trained solely D; tiny Pomeroy OOF rare-class cases may use *explicitly marked* ridge surrogate, never concealed as exact original. Exactly 7 arms x 11 cohorts = 77 patient-replayed heldout V cohort-model summaries expected.
- Record every patient's real probability, decision, path, actual number of relational comparisons, target and terminal leaf; every node's two pair IDs, β weights, children, training group sizes and class posterior. Check 0/2/4 path tests for weighted, ≤2 for label-only TSP proxy; no global RLR teacher needed at prediction.
- Independent Python auditor reruns each metric from heldout V patient scores and verifies all D/V sample IDs are disjoint, all 11 named cohorts and full 77 methods, no outer E usage, tree node shape and exact paths. Exact 2^11 signflip and seeded dataset bootstrap across eleven independent dataset-mean units (NOT 77 independent groups or 470 independent folds), with no p-value cherry-picking.
- Report F1-positive, MacroF1, BA, AUROC, logloss, Brier, student agreement/fidelity vs teacher, mean/max evaluated pairs, number of unique pair relations and genes in full tree, number of internal nodes and leaves. Negative results preserved.

## Publication and novelty gates

This is **NOT** a priori novel merely because it has *two* RXA tests in nodes, logistic coefficients or a greedy tree: the user's earlier TSPDT, GTSPDT, 2020 tree-based ARXA fractional comparisons and EMTTree+RX already cover extensive multirelational tree features. The potential increment is to use original globally sparse **RLR coefficients and ternary activation** in a strictly tiny standalone tree, distinguishing beta magnitudes vs signs and coefficient selection vs genuine out-of-fold teacher effect. This must be shown as a distinct practical or inferential contribution against **the original published TSPDT and EMTTree+RX**, not only our `LABEL_BINARY_TSP_TOP20` proxy.

Proceed to first outer E only after positive or informative predeclared D/V results, as additional **exploratory** evidence on already-exposed source cohorts. A generalization/publication claim requires new unexposed data or original fully verified source ARFF, calibrated probabilistic controls, and an audit of class coding, sample IDs and paired inference.

If v4 loses to simpler sign-only/label-only TSP node with equal complexity, freeze a negative science gate. Do not retune α, β or K using previously seen E.

**Actual code paths:** `R/experimental/rlr_weighted_pairpair_tree_core_v4.R`, `R/experimental/test_rlr_weighted_pairpair_tree_v4.R`, `R/experimental/rlr_weighted_pairpair_first11_DV_v4.R`, `scripts/rlr_weighted_pairpair_first11_DV_audit_v4.py`.
