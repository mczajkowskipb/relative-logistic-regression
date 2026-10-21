# TSPDT reconstruction v1: post-result methodological audit

**Audit date: 2026-10-21. Status: FAIL as faithful TSPDT-2011 reproduction;
PASS only as a documented independent, deterministic, testable tree baseline.**

A completed, independently replayed 95-fold result from this implementation
was substantially weaker than RLR on most publication datasets. This does
**not** establish that the original 2011 TSPDT was weak or that R as a
programming language cannot implement it efficiently/correctly.

## Evidence: the original paper versus our source

The primary 2011 paper:
https://wi.pb.edu.pl/p/wp-content/uploads/sites/3/2024/07/AEMB2011.pdf

The model as actually tested is
\`R/baselines/tspdt_reconstruction_v1.R\`.
No source of the actual original Weka/Java TSPDT implementation is available.

| Issue | 2011 source | Our 95-fold reconstruction | Consequence |
| --- | --- | --- | --- |
| Feature selection | **Supervised Relief-F**, 10 neighbors, top **1000** genes (§3) | **Unsupervised MAD**, top **100** genes | Different feature universe; discards class-discriminative low-MAD genes; at most 4950 gene pairs instead of 499500 |
| Root split | **Identical to original TSP** on its selected feature space (§2.2) | Best TSP-scoring **feasible** pair subject to \`min_leaf=5\` | Root can differ from TSP even on the same features |
| Score ties | Original TSP-family may use secondary ranking; 2011 paper does not explicitly pin all tie rules | First pair encountered in MAD order | May select different gene pairs and distinct subtrees |
| Growth/stopping | Recursive node-specific split, details not fully enumerated | **max_depth=3**, \`min_leaf=5\`, stop at pure nodes | Artificially limits accessible trees, especially on imbalanced samples |
| Pruning | C4.5-like pessimistic pruning (§2.2) | Sum of training-error Wilson upper bounds at cf=0.25 | Heuristic is **not** C4.5's exact subtree error/standard error mechanism; may prune viable splits |
| Node tests | TSP, k-TSP, weighted TSP-family alternatives (§3) | **Single TSP** pair per node only | Does not represent the broader TSPDT family |
| Prediction | Historical hard class decision | Leaf-majority class and \`alpha=1\` Laplace score | Laplace logloss is reconstruction-specific, not original 2011 |
| Model comparisons | Historical fixed holdout/10×10 CV, dataset-specific structures | 19 RLR publication matrices, only 1×5 | Historical accuracy table is **not** directly comparable |

The existing TSPDT code is **not necessarily computationally wrong**
in selecting its restricted node-local TSP test: on synthetic data it
reproduces its own declared score, routing, pruning and RDS predictions.
But those synthetic tests **do not test equivalence** to the 2011 model
and cannot establish that its arbitrary design shortcuts preserve
predictive accuracy. One specific conceptual failure is that the
reconstructed root need not equal unrestricted TSP on the same
selected features, contrary to the original paper.

Current results are therefore most precisely named:
**\`TSPDT2011_INSPIRED_RECON_V1_MAD100_DEPTH3_WILSON\`**.
No plot, table or paper should label them unqualified as
"TSPDT-2011" or "replicated TSPDT".

## Separate defect in the four-model REPORT (not in R model fits)

The first attempt to build a unified four-method table failed on
\`west-2001\` with \`original ref audited metric mismatch (...,'RLR','logloss')\`.

Root-cause hypothesis pending rerun: the Python comparison script
previously used \`tuple(sorted(set(y_true)))\` to assign class0
and class1. However R's \`foreign::read.arff\` retains the
\`@attribute\` nominal factor level order, which can be \`{1,0}\`
or some non-lexical order. In those datasets Python can **invert**
the meaning of saved \`prob\` (probability of the *second R level*).
This invalidates recomputed logloss **without altering RLR accuracy
or any fitted model**.

Fix: derive the **unique original two-class orientation** using
the independently verified RLR invariant
\`y_pred = class1 if prob >= 0.5 else class0\`,
check it against *every* reference prediction, then use this
order consistently for all 380 already-saved models. An explicit
synthetic reversed-order regression test is now present.

This fix must pass the real-data four-method report before the
underlying logloss mismatch is regarded as confirmed resolved.
No tolerance threshold is loosened and no original CSV modified.

## What to do next — in the correct scientific order

1. **Complete read-only diagnostics**: sizes of all 95 existing
   trees, fraction collapsed to root, depths of internal nodes,
   repeated gene usage and test class distributions. This
   identifies whether underperformance is driven primarily by
   trivial trees, weak feature selection or unstable routing,
   without opportunistic refitting.
2. Freeze the 95-fold negative comparison under v1; **do not
   retrospectively optimize MAD-100, min_leaf=5, pruning,
   max_depth or tie rules on these same outer-test results**.
3. Specify **new** reconstruction scientifically from source:
   a train-fold-only Relief-F 1000 implementation, TSP with
   primary and secondary tie ranking, root TSP invariance,
   and a correctly defined pessimistic C4.5-like pruning
   algorithm. Write unit tests comparing with trusted
   TSP/Relief-F reference implementations where possible.
   Do not invent undocumented historical constants.
4. Solve high-dimensional pair-ranking efficiently:
   chunked evaluation / compiled array operations;
   no dense \`n×choose(p,2)\` representation in RAM.
   Reuse mature libraries for Relief-F, trees, fold generation
   and metrics when semantically appropriate.
5. Call the improved algorithm **\`TSPDT2011_RECON_V2\`**,
   not a recovery of the 2011 original source. Run synthetic,
   leakage and complexity gates **before** real-cohort fitting.
6. The existing 19 publication datasets are now part of model
   development. If a v2 design is compared there, label the
   analysis as exploratory and use a separate fixed v2 output
   root. Preserve v1 and the original RLR publication unchanged.
   A publication-facing superiority claim requires independent,
   never-used biomedical data.

**Decision:** V1 model code and its negative benchmark are frozen.
New reporting/audits may improve; changes to model semantics
belong in a separately specified V2, never in the existing v1 source.
