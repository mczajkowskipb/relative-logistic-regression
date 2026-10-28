# RLR-DT v1 — Compact OOF RLR-guided ordinal relational tree

**Date:** 2026-10-28, first experimental implementation. **Original unpublished RLR manuscript/source is locked and untouched.** New branch: [analysis/rlr-distilled-compact-rank-tree-v1](https://github.com/mczajkowskipb/relative-logistic-regression/tree/analysis/rlr-distilled-compact-rank-tree-v1).

## Scientific question

Can the original RLR20 **teach a stand-alone, shallow relational decision tree** to preserve much of RLR prediction quality while requiring only **0–2 within-patient pair tests along the path**, rather than checking ~16.63 globally retained RLR20 ternary relations per patient? Is any benefit due to the information from RLR, not merely using good pair candidates? Can the method yield smaller explanations than a TSPDT-like label-greedy control with the same tree budget?

**No claim yet that this method outperforms or is novel relative to the author's published TSPDT/GTSPDT**. TSPDT (Czajkowski & Krętowski 2011, DOI 10.1007/978-1-4419-7046-6_3) used pairwise gene-expression comparisons at tree nodes. GTSPDT (Czajkowski & Krętowski, EuroGP 2013) uses global evolutionary optimization of tree structure and multirelational tests. More recent rank-based tree ensembles also exist (J Transl Med, 2024, DOI 10.1186/s12967-024-04940-2). Those publications must be properly positioned in any article, with an implementation-level comparison before claiming true methodological novelty.

## Key methodological distinction from old TSPDT and original RLR

1. **Teacher**: literal original frozen RLR20 trained on class labels, using original switchBox top1000, `apply_disjointness(d=2)`, patient rank ternary margin `m=10`, elastic-net `alpha=.5` and `cv.glmnet(lambda.min)`.
2. **Honest soft labels**: for every **D-only** patient, obtain the original RLR20 predicted class1 probability from a model in which that patient is HELD OUT by a stratified, threefold teacher crossfit. No in-sample RLR teacher probability is used as a distillation target. Teacher models perform their own switchBox candidate screening inside respective teacher train folds. Full teacher RLR20 fitted on D separately for V comparison and for specifying the **RLR-selected candidate node panel**, which is legally learned from D only.
3. **Small independent student**: up to 20 candidate pairs from the original RLR's globally coefficient-weighted shortlist, not all 1000 pairs. Two predestination node forms:
   - Binary: `rank(gene_right)>rank(gene_left)` (2 branches).
   - Ternary: `rank(gene_right)-rank(gene_left)>=10`, `<=-10`, or inside the dead band (3 branches). This deliberately mirrors original RLR ternary activation/inactivation, not just a standard TSPDT binary test.
4. **Teacher-guided leaf evidence**: each D patient gets a mixed target `z_i=(1-α)y_i + α p_RLR20^{OOF}(x_i)`, **α=0.55 fixed**. A leaf has smoothed `p_l=(sum_{i∈l} z_i +4p_0)/(n_l+4)`, where `p_0` is smoothed D-global target prevalence. Greedy best-first node split maximizes D aggregate reduction in proper **binary cross-entropy of soft targets**, penalized by **0.003 per added internal split**. It does NOT do class entropy or TSP score only. Controls α=0 use the SAME cost-complexity objective and exactly same pair candidates, isolating whether soft RLR OOF teacher information adds value rather than just acting as supervised feature nomination.
5. **Hard complexity constraints**: `max_depth=2`, `max_internal_nodes=3`, `min_leaf=max(4,floor(.07n_D))`. At most 2 pair comparisons on any patient path, at most 3 distinct tree test pairs and at most 6 genes in tree, binary up to 4 leaves or ternary up to 7. Early termination/pruning produces smaller trees (possibly a zero-split majority/probability leaf). None of these bounds depends on V or E outcomes.
6. **Autonomous inference**: the student contains only the chosen pair gene names, tree branches and smoothed leaf posterior. It does **not** contain, load or run the RLR teacher or its original large gene-pair panel. The inference algorithm ranks the same whole measured assay gene universe within the patient and **evaluates only each visited node pair**; the underlying array still measures all assay genes. Ranked-feature tree decisions are invariant to **each sample's individual strictly increasing nonlinear transformation**, given the same assay feature IDs.
7. **What the tree is NOT**: not an RLR model attached to each leaf, not an RLR stacked with TSP voting, not a claim of reduced wet-lab assay costs, not an original TSPDT source reproduction, and not an independently validated new clinical diagnostic.

## Frozen pilot and actual controls

[11-dataset exact first outer T D/V trial](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40658027175), **original E never accessed** by this study. Pilot compares **six** methods on exactly the same V patients, with original repeat/split protocol:
- `RLR20_TEACHER` (original RLR20 trained D, V prediction);
- `KTSP5_NATIVE` (real switchBox original kTSP5, native labels);
- `LABEL_BINARY_SAME_RLR20`: α=0 label-only compact binary tree, same 20 RLR-selected candidate pairs;
- `LABEL_BINARY_TSP_TOP20`: α=0 label-only compact binary tree, 20 original top switchBox pairs; **TSPDT-like greedy proxy only**, not an exact reproduction of original TSPDT;
- `DISTILL_BINARY_RLR20`: α=.55 compact binary relational tree;
- `DISTILL_TERNARY_RLR20`: α=.55 compact original-m10 ternary relational tree.

Report positive-class F1, Macro-F1, BA, AUC, logloss, Brier, agreement with original RLR20 V class decisions, soft CE and Brier fidelity to teacher V probabilities, mean/max path comparisons and actual patient path strings, internal node count, terminal leaves, unique genes. Node table and heldout patient predictions archived for every method. Independent Python replay of `6×11=66` method/cohort V scores, strict D/V no overlap, path bound ≤2, plus unit mechanics on monotonic patient transforms, exact lack of teacher at inference. All source datasets are previously explored Schliep public reconstructions, not original processed paper ARFF SHA matched; all pilot conclusions exploratory.

### Go/no-go for a habilitation-compatible new TREE

The **teacher-guided tree** is interesting if at least one of the frozen binary/ternary variants retains reasonably close F1/AUC to RLR20 or actual kTSP5 while using <2 actual node pair checks on average, **and** improves over BOTH fixed-capacity label-only controls on more than just cherry-picked single cohorts. Must report the model's *cost and faithfulness vs teacher* AND *true accuracy against class labels*. Report equal-weighted 11 datasets and per-cohort heterogeneity.

If all teacher-guided trees lose badly to α=0 label-only and kTSP5, the specific OOF distillation architecture fails; do not force a positive spin. A novel tree will then require a different objective/conditional representation, not a decorative "RLR teacher" box. Do not promote it to 470 E folds until stable D/V evidence, and do not claim numerical superiority or generalization based on previously exposed cohorts.

## External literature pointers

- Czajkowski, Krętowski, *Top scoring pair decision tree for gene expression data analysis*, Advances in Experimental Medicine and Biology **696** (2011), 27–35, DOI: 10.1007/978-1-4419-7046-6_3.
- Czajkowski, Krętowski, *Global Top Scoring Pair Decision Tree*, EuroGP 2013, original published full text https://www.aragorn.wi.pb.edu.pl/~mczajk/eurogp2013.pdf.
- *Ensemble methods of rank-based trees for single sample classification with gene expression profiles*, Journal of Translational Medicine **22** (2024), DOI 10.1186/s12967-024-04940-2.
