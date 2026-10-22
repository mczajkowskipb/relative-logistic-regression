# RLMT-v3 — tree redesign after honest RLMT-v2 negative evidence

Research checkpoint: 2026-10-22.

## Why model family was redesigned

The previous fully paired retrospective benchmark was 19 previously
exposed cancer datasets, 5 original folds per cohort (95 fold positions),
four methods and 380 real R fits, independent saved-prediction audit PASS:

| Model | Accepted/95 | Equal-dataset mean logloss vs original RLR |
| --- | ---: | ---: |
| H1 RLR leaves | 9 | +0.014324918 |
| H2 ridge offset | 25 | +0.000992198 |
| TSPDT_GINI, REAL RLR leaves | 7 | −0.001087952 |
| RLMT_LA, REAL RLR leaves | 0 | +0.001236397 |

Negative logloss difference is better. TSPDT_GINI dataset-level
7 win, 10 loss and 2 tie vs published archived RLR. Leave-one-
dataset-out estimate changed sign for 4 of 19 datasets. A mean
benefit of ~0.0011 dominated by a few cohorts is NOT a credible
superiority claim. H1 was harmful; RLMT_LA is inert and costly.
Across all 380 fits glmnet flagged 4209 small-class warnings.
All were already development-exposed, not independent tests.

Do NOT represent this new work as an improved proven method,
or as an unbiased independent method comparison on these data.
No cherry-picking a positive method/endpoint and retroactively
labelling it prespecified validation.

## Methodological change: separate routing from leaf estimation

Previous prototype used a fixed 4096-random-pair proposal, imposed
min 16 samples **and 6 of each class in every leaf**, used a newly
fitted full RLR in small leaves, and then required separate
validation logloss improvement >1SE before accepting a split.
This can create zero-split models even when a simpler small
tree would make sense.

The new methods ARE DIFFERENT scientific hypotheses, declared before
running the new 18-fit pilot. They are intentionally not parameter
fine-tunes of the v2 methods. They do not guarantee better accuracy
or logloss than the frozen published RLR.

### Variant A — STUMP_TSP_BETA

A one-rule tree with a split \`x[a] > x[b]\`. Candidate selection on
outer training only, using the magnitude of the TSP class contrast
among a deterministic 4096-**sampled** relational pair inventory.
At least 8 observations are required on each side, but individual
sides may contain only one class: no glmnet is fit inside a leaf.
Prediction is a hierarchical Beta-binomial shrinkage of leaf
class prevalence toward its parent probability:
\`p_leaf = (n_positive + tau * p_parent)/(n + tau)\`,
with \`tau=16\`. A stump is trained whenever a feasible
informative route exists. No 1SE rejection relative to RLR.

### Variant B — STUMP_RIDGE_OFFSET

Same one-rule TSP-scored routing and 8-sample leaf minimum,
but the output is \`logit(p_RLR(x)) + b_false*I(false) +
b_true*I(true)\`. Only **two coefficients**, jointly ridge
penalized with \`lambda=1\`, are fit. Aligned parent training
out-of-fold predictions are used to construct the training
logit offsets (no in-sample offset leakage). No per-child
supervised pair reranking, no internal child cv.glmnet.
One valid fitted decision stump need not improve over original
RLR; the heldout outcome will measure that truth.

### Variant C — TREE_D2_PRUNED_BETA

Hierarchical rank-relation tree of maximum depth **2** and at
most **3 internal nodes**. Each node searches its own TSP
contrast on the training subset, using 4096 deterministic pair
candidates and no per-class minimum in the children. Each
terminal node's Bernoulli probability shrinks hierarchically.
In the outer-training 2/3 discovery and 1/3 confirmation:
build depth2 on discovery, compare heldout internal
logloss at depth1 versus depth2, and retain depth2 only if
mean validation loss improvement exceeds its SE.
Otherwise **prune to a stump**. Refit the selected depth
on all outer-training samples. There is intentionally no
1SE requirement against RLR: this is a pure tree family.
The outer-test cohort/fold is untouched until final scoring.

### Not faithful replication of historical TSPDT

The earlier TSPDT (Czajkowski & Kretowski, *Top Scoring Pair
Decision Tree for Gene Expression Data Analysis*, AEMB 696,
2011, DOI 10.1007/978-1-4419-7046-6_3) may use one or more
TSP/k-TSP pairs in internal nodes and C4.5-like pessimistic
pruning. This prototype currently uses one TSP-scored
*sampled candidate pair* per node. It must NOT be called
a faithful reproduction of the original TSPDT or a new
method established by that prior paper.

The NEXT exploratory family, contingent on actual v3 pilot
results, is a **faithful old TSPDT reimplementation**:
a) exact train-only top pair (full candidate search or
a justified provably safe approximation), b) k-TSP
multi-pair voting node, c) C4.5-like pessimistic
postpruning. Benchmark as historical comparator, not
new contribution. If its old implementation is found,
prefer reusing and auditing it.

## Frozen small gate / larger research study

Gate 0: parse R sources; real synthetic 180/60 training/test
for all three modes; exact RDS prediction replay;
nondegenerate smoothed probabilities and train-only
D/V pruning.

Gate 1: two exact previously exposed cohorts \`chen-2002\`
and \`gordon-2002\`, first 3 original publication fold
positions, 3 methods = 18 NEW real fits, archived original
RLR (285 total original benchmark records) reused read-only.
Independent audit verifies paired sample IDs, original
RLR scores, aligned global parent, logloss, Brier,
model route sizes and captured warnings.

Only after Gate 0+1 PASS and inspecting numerical
failures, not just favorable scores, consider an exact
19x5 cohort benchmark in a **new untouched output root**.
Nothing yet run, no superiority results asserted.

The needed real validation is a newly sourced patient
cohort or laboratory batch held aside **before method
selection**. Exposed test fold performance is not that.

Main RLR paper remains intact and publication-ready.
New model-tree results could motivate a separate paper,
not retroactive rewriting of already reported RLR effects.
