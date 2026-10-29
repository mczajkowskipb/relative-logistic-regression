# RLR-GLVT v5 — global logit-variance budgeted relational tree (PRE-RESULT PROTOCOL)

**Initiated 2026-10-29. STATUS: exploratory algorithm proposal only. Not proof of novelty or clinical utility.**

## Motivation and boundary
RLR-DT v1–v3 (OOF teacher cross-entropy guiding greedy splits) and RLR-WT v4 (two signed logistic-weighted ternary features per node) did not show stable scientific benefit on 11 repeatedly exposed Schliep cohorts. V4 F1=0.758 vs label-only TSP-like F1=0.840 and 2.43 vs 1.30 path comparisons; v4 is frozen NEGATIVE. Never revisit v4 alpha, beta, K, threshold on exposed V/E.

This new method is **not** C-RET: it does not evaluate/short-circuit an existing RLR model at inference and does not claim to preserve any individual parent decision or guaranteed probability interval. It is **not** another pair-pair feature. It learns one binary rank-comparison per internal node; original RLR contributes **only** an OOF log-odds homogeneity penalty to the *global* tree induction objective. Prediction is autonomous, leaf probabilities from true D labels alone. This design is distinct mechanically but not automatically novel: globally optimized trees, teacher/student surrogates and cost-sensitive paths predate this work. Prior art review required before any manuscript claim.

## Frozen objective and exact search
For each training-only D patient obtain truly out-of-fold original RLR20 probability p_i (same train-only 3fold procedure already audited; report surrogate fallback folds in rare-class Pomeroy), and define t_i=logit(clip(p_i,1e-4,1-1e-4)). Do not substitute in-sample parent predictions. For a partition into leaves l, optimize:

J(T) = (1-w) * sum_l CE_y(l) / CE_y(root)
     + w * sum_l SSE(t_i : i in l) / SSE(t_i : all D)
     + 0.005 * number_of_internal_nodes
     + 0.005 * mean_training_path_length.

CE_y uses smoothed true-label leaf probabilities p_l=(positive_l + 4*p_D)/(n_l+4), with p_D=(positive_D+1)/(n_D+2). If OOF logit variance is zero, set w_eff=0. Fixed w in {0, 0.25}, **not tuned on already seen V/E**. w=0 is class-only control with *identical* candidates and exact search. SSE term minimizes teacher logit heterogeneity, not teacher probability cross-entropy or label blending. The normalizations establish same dimensionless root loss for both terms, not a promise that their statistics are equally informative.

One rank pair predicate at each node: I(rank(g_right)-rank(g_left)>0). Candidate set: train-D switchBox TOP20 original TSP pair shortlist, **IDENTICAL for w=0 and w=.25**. This isolates effect of original RLR OOF *conditional logit variance* on topology; RLR NOT used for node shortlist in primary ablation. Enumerate all valid roots and the optimal optional split independently on each root child. This exactly minimizes J among <=3 internal nodes and depth<=2 for the *fixed pool*, >=max(4, floor(.07*n_D)) samples per leaf; 0-node model allowed. Cap <=2 actually evaluated rank-pair comparisons per patient, <=3 distinct model pair rules and <=6 model genes. The full transcriptomic rank universe is still required for within-patient ranks: <=2 comparisons is NOT 2-pair wet-lab assay.

## Planned preregistered study — NO outcome yet
Reuse the frozen original outer-first-T -> stratified D/V partition of eleven paper-named cohorts. All pair selection, OOF cross-fitting and optimization use D ONLY, E never accessed; RLR20 full-D and real switchBox kTSP5 continue as context, greedy TSP-top20 label-only small tree as an additional reference. Report at least:
- GLOBAL_LABEL_ONLY (w=0), GLOBAL_RLR_LOGIT_OOF (w=.25), GREEDY_LABEL_TSP_TOP20, original RLR20, native kTSP5.
- Positive F1, Macro-F1, BA, AUC, logloss, Brier; soft fidelity and label agreement vs parent original RLR20, path comparisons and full-tree pair/gene counts, no-split frequency, fallback folds.
- Independent replay of all held-out V patient predictions and D/V identifiers, 11 equal-weighted cohorts, exact 2^11 sign flips and dataset bootstrap. Label-only GLOBAL must use the same branch candidates, cost terms and leaf smoothing; a greedy control alone would confound optimizer gain.
- Critical science gate: same-pool RLR-guided global optimum improves useful accuracy/fidelity at comparable mean path cost across several cohorts, not just Chen/Pomeroy. A gain only in teacher mimicry without meaningful classification/interpretability benefit is insufficient.
- Before publication: historical **actual** CUDA TSPDT/GTSPDT/EMTTree+RX comparison, original processed ARFF probe/SHA equivalence, and genuinely not previously observed biological cohort. Current Schliep data are development and heavily re-used.

## Literature novelty audit (no novel claim)
Bertsimas & Öztürk (2023), Global optimization via optimal decision trees, DOI 10.1007/s10898-023-01311-x; loss-optimal logistic classification trees, DOI 10.1007/s11750-024-00674-y; sparse optimal trees and selective methods. Distinguish from author's TSPDT/GTSPDT/EMTTree+RX and former C-RET in detail. Global search, path budget or RLR teacher alone is known/insufficient. Original unpublished RLR and manuscript lock remain unchanged.

**Code:** R/experimental/rlr_global_logit_variance_tree_core_v5.R
**Unit:** R/experimental/test_rlr_global_logit_variance_tree_v5.R (synthetic mechanics only)
**Next gate:** actual PASS in Actions first, then original11 D/V paired run + independently replayed patient-level audit. No 470 E until the above scientific gate passes.
