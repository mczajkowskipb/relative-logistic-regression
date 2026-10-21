# TSPDT-inspired reconstruction v1 — model specification and DEV contract

**Stage:** synthetic DEV only; NOT an original historical 2011
implementation and NOT yet benchmark-approved. This file is stored
on a separate GitHub branch so the 285-fit RLMT-DC process remains
completely untouched.

## Evidence and limitation

Primary method description: M. Czajkowski, M. Krętowski,
*Top Scoring Pair Decision Tree for Gene Expression Data Analysis*
(2011), https://www.aragorn.wi.pb.edu.pl/~mczajk/AEMB2011.pdf,
Sect. 2.2 and 3, Fig. 1, Table 2.

The paper states recursive, greedy node-specific TSP/k-TSP ranking
and C4.5-like pessimistic pruning; its root is the TSP decision rule
on the full training set. It used Relief-F with 10 neighbors to
select 1,000 features for historical classification accuracy
comparisons. It does **not** completely specify tie-breaking,
exact C4.5 pruning confidence constants, leaf posterior
probabilities, or all original Java/Weka implementation details.
Consequently the original 2011 model cannot be claimed reproduced.

## Implemented approximation

Source: `R/baselines/tspdt_reconstruction_v1.R`.

- API: `tspdt_recon_fit(X_train,y_train)` and
  `tspdt_recon_predict(model,X_test)`.
  `X_test` must match the training numeric feature schema;
  the prediction function accepts **no y_test**.
- Only training data chooses **top 100 MAD-ranked genes**.
  In contrast, 2011 used **supervised Relief-F, f1000**.
  The reconstruction is therefore a *fair, low-complexity
  comparator*, **not** a recreation of the paper's specific
  feature selection.
- For each node containing labels `C0,C1` define a new TSP
  score `S(a,b) = |Pr(x_a<x_b | C0) -
                   Pr(x_a<x_b | C1)|`.
  Search all unordered gene pairs in the selected feature set,
  **recomputing probabilities on the samples routed to that node**.
  Training comparisons use strict less-than; equality routes
  to the `FALSE` branch. No labels from the outer holdout appear.
- Guardrails: require **min_leaf=5** samples in each proposed
  child, `S>0` and `max_depth=3` (root depth=0).
  For score ties within 1e-12, use the first pair in the
  training MAD-ranked gene order. A root split may therefore
  differ from the unrestricted one-pair TSP when the highest
  scoring relation would violate the reconstruction's
  child-size gate. **This is an explicit non-equivalence**.
- Terminal class = training-sample majority, with ties assigned
  to the training factor's first level. Training-only
  Laplace class-1 score `(n1+1)/(n+2)` is stored
  separately from the original 2011 class rule. It is
  an *uncalibrated reconstructed leaf score*, not a
  historically published TSPDT probability.
- Optional post-pruning uses a **Wilson upper estimate** for
  the node training misclassification fraction, confidence
  `cf=0.25`. Collapse a subtree when the parent estimate
  does not exceed the sum of child estimates. This is a
  transparent training-only **pessimistic-error heuristic**,
  **not the source-equivalent C4.5 implementation**.
  For comparability, `pruning="wilson"` is the registered
  default; disabling pruning is reserved for synthetic tests.
- Pair screening costs O(n*p²) at each visited node but holds
  only an O(n*p) temporary boolean matrix by iterating the
  first gene; no `n × choose(p,2)` global matrix is created.
  With p=100, there are at most 4,950 comparisons per node.
- Recursive tree nodes retain the TSP rule, node-local score,
  training-class counts, leaf estimate, child routes and
  leaf probabilities. Model is standard `saveRDS` compatible.

## Separate BigTSP baseline

Source: `R/baselines/bigtsp_pair_cart_v1.R`.

- Uses only the **reviewed BigTSP version 1.0**, if preinstalled,
  with its exact source `BigTSP::tsp.tree` and
  `tree.control(mincut=2,minsize=5,mindev=0.001)`.
- Enforces a **hard 100-gene cap** before BigTSP constructs
  the full 4,950 pair-feature binary columns; no full
  high-dimensional all-pairs expansion.
- Supplies a **dummy response column** to BigTSP's unusual
  `predict.tsp.tree` interface, never test labels.
- Reports **class predictions only** and NA probability
  rather than pretending to have calibrated logloss values.
- Identity: `BIGTSP_1_0_PAIR_CART_F100`, not TSPDT 2011.
  See `docs/BIGTSP_1_0_SOURCE_AUDIT.md`.
- R package license GPL >=2; dependency verification and
  source/installation provenance must be recorded before
  publication-facing redistribution.

## Synthetic gate (not yet executed on the Ubuntu compute host)

Run `bash scripts/tspdt_reconstruction_devcheck.sh`, which checks:

1. Parse all three R source/test files without touching omics.
2. A controlled two-class synthetic data design forces a
   TSP root and distinct **fresh TSP rules in both children**.
3. Deterministic model object, selected genes, node counts,
   held-out routing, equality convention and RDS replay.
4. Strict schema rejection and input feature count cap.
5. Training-only feature filtering is insensitive to label
   inversion. Prediction functions have no access to test labels.
6. Finite nonextreme Laplace estimates; pruning never
   expands a tree.
7. If BigTSP is already installed, smoke the separate
   class-label adapter; otherwise explicitly skip it.
   The main TSPDT-inspired DEV result can still pass.

This DEV gate requires a real `Rscript` run. A GitHub commit or
simple static bracket check **does not** mean the DEV gate passed.

## Future fair publication-data evaluation — NOT launched

Once the currently running 19-dataset × five-fold RLR/RLMT/RLMT-DC
benchmark exits successfully AND all four independent audits
pass, design a new frozen experiment on **the same 95
train/test assignments** for `TSPDT2011_INSPIRED_RECONSTRUCTION_V1`,
and optionally a distinct `BIGTSP_1_0_PAIR_CART_F100`
if source/dependencies validate.

Include all previously examined Gordon, Singh, leukemia, breast2
and the remaining datasets. Do not select 100/200/1000
features, depth, pruning threshold or probability smoothing
based on their test performance. Current parameters
(f100,minleaf5,depth3,cf0.25,Laplace1) are registered BEFORE
any TSPDT-inspired real-cohort scores. Treat logloss as a
**reconstructed leaf-score** comparison; balanced accuracy,
macro-F1, accuracy, tree nodes/leaves and gene relations
are primary descriptive model-comparison measures.
The original 2011 Table 2 is not numerically comparable
because data splits and preprocessing differ.

No old publication/RLMT model code, branch or result directory
is modified by these baseline implementations.
