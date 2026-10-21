# BigTSP 1.0 `tsp.tree`: actual-source verification

**2026-10-21 — inspected read-only archival R source.**
This is a **model identity assessment**, not an empirical benchmark run.
**Verdict: BigTSP `tsp.tree` is NOT an implementation of the
Czajkowski–Krętowski 2011 node-local greedy TSPDT algorithm.**

## Reproducible provenance

- User-supplied manual:
  https://www.rdocumentation.org/packages/BigTSP/versions/1.0/topics/tsp.tree
- CRAN source archive:
  https://stat.ethz.ch/CRAN/src/contrib/Archive/BigTSP/
- Directly inspected public CRAN mirror:
  https://github.com/cran/BigTSP/blob/af8320b2788076657556728b25c53ed04681c7d5/R/LDCA.R
- Mirrored commit:
  `af8320b2788076657556728b25c53ed04681c7d5`;
  `R/LDCA.R` Git blob
  `efeed33f102e9864e37341da53d50273fed0f3f3`.
- Package `DESCRIPTION`: version **1.0**, published August 2012,
  authors **Xiaolin Yang, Han Liu**, license **GPL (>=2)**,
  dependencies `glmnet, tree, randomForest, gbm`.
- CRAN delisted the package in 2022 (maintainer email reportedly
  undeliverable); this does not itself establish an algorithm defect.
- **Constraint:** this audit reads the GitHub source mirror. The raw
  archive tarball was *not downloaded* in this assistant runtime,
  so no archive SHA256 or successful local R installation is claimed.

## What the actual implementation does

Source `R/LDCA.R`:

1. Lines **6–21**, `tspfeature(X)`, creates a binary column for
   **every unordered gene pair**, via the strict rule
   `as.numeric(X[, i] < X[, j])`. Feature transformation is done
   **globally once**, before tree induction. It is independent of
   the response labels but quadratic in input feature count.
2. Lines **99–107**, `tsp.tree(X,response,...)`:
   `newx <- tspfeature(X)`,
   `d <- data.frame(response,newx)`,
   `rt <- tree(response ~ ., data=d, control=control,
              method=method, split=split,...)`.
   This uses the standard R **`tree::tree` impurity criterion**
   (`deviance` or `gini`) over binary pair features.
   **No TSP pair score or recursive node-local TSP candidate
   reranking is implemented by this wrapper.**
3. The resulting model has class
   `c("tsp.tree","tree")`; its pruning/stopping rules are those
   of `tree.control` and `tree`, **not** the
   source-described C4.5-like pessimistic pruning of TSPDT-2011.
4. Lines **110–117**, `predict.tsp.tree`, **expects the first
   column of `newdata` to contain a response placeholder** and
   uses `tspfeature(newdata[,2:ncol(newdata)])`.
   Passing just an unlabeled predictor matrix would
   accidentally discard its first predictor. This is a legacy
   API hazard. The value of the response column is not used for
   constructing the pair comparison matrix.
5. Since this wrapper materializes `p(p-1)/2` pair columns,
   it is not directly feasible as an unconstrained baseline on
   every original high-dimensional publication ARFF. At
   `p=2000`, this is **1,999,000 binary columns** before R
   data-frame duplication and tree working memory.
   At `p=200`, **19,900**; at `p=100`, **4,950**.
   These are column counts, not verified actual memory usage.
6. The standard `tree` object may include fitted node
   `yprob` values, but these are empirical tree probabilities,
   **not calibrated** and not a result from the original
   TSPDT 2011 implementation. To compare logloss, probability
   extraction and smoothing must be prespecified inside
   training data without test tuning.

## Comparison against the 2011 TSPDT design

| Property | 2011 TSPDT original | BigTSP 1.0 `tsp.tree` |
| --- | --- | --- |
| Gene-pair comparison test | Yes | Yes |
| Pair proposals reevaluated on node-local subset | Yes | **No explicit TSP reranking** |
| Internal-node selection | TSP-family score / accuracy | `tree()` deviance or Gini |
| Pruning | C4.5-like pessimistic pruning | `tree()` controls/pruning |
| Binary pair matrix | Not required as one global dense design | Full `p(p-1)/2` data frame |
| Output class | Historical custom TSPDT classifier | R `tree` wrapped as `tsp.tree` |

The title and citation in package documentation **are not evidence
of source-level methodological equivalence**. Appropriate distinct
label: `BigTSP_pair_Cart` or
`BigTSP_1.0_tree_on_all_pair_features`, not
`TSPDT_2011_REPRODUCTION`.

## Additional discovery: the PI's own AIR benchmark already uses BigTSP

The existing repository
https://github.com/mczajkowskipb/AIR-relational-benchmark
contains a working-looking benchmark wrapper at

https://github.com/mczajkowskipb/AIR-relational-benchmark/blob/main/R/methods/method_bigtsp.R

Its code:

- Calls `BigTSP::tsp.tree` with `tree.control(mincut=2,
  minsize=5,mindev=0.001)`.
- Correctly works around BigTSP's legacy prediction interface
  by passing a first `response` column. However it supplies
  the test label `y_test` to that column. The BigTSP source
  does **not** use this value in pair construction, but a
  dummy response column of the correct type would be safer
  for a future adapter.
- Uses `type="class"`, so **does not claim calibrated
  probabilities**, and records score as NA.
- Supports bounded prior filtering: the AIR benchmark main
  protocol is train-only MAD feature selection with f200;
  it also has a `scripts/23_make_tspdt_f100_grid.R`
  specifying `n_features=100` for `tspdt_bigtsp`.
- The AIR repository's `METHODS.md` currently calls this
  method `TSPDT` and `config/methods.csv` assigns ID
  `tspdt_bigtsp`; that **label risks overstating fidelity**
  to the paper. When updating that separate repository,
  distinguish `BigTSP pair-CART` from original TSPDT.
- The AIR wrapper's existing fit/predict plumbing can be
  **reused**, rather than reimplementing BigTSP or using a
  full `p=1626+` pair matrix on the RLR paper datasets.
  No source from the other research repository has been
  copied or modified in this check.

AIR is a distinct publication companion project with **different
omics datasets and CV protocol**, not evidence of RLMT benchmark
performance on the 19 RLR-paper matrices.

## Recommended next steps without touching running RLMT jobs

**Option A (recommended):** implement a fresh, clearly identified
`TSPDT2011_APPROX_V1` algorithm from the original paper,
using node-local TSP scoring, deterministic class/tie conventions,
train-only preprocessing and documented pruning substitutes
where original source details are unavailable. Never call it an
exact reproduction absent source-level equivalence evidence.

**Option B (additional comparator):** reuse the AIR
`BigTSP_pair_Cart` wrapper with train-only feature filtering
to at most 100 or 200 features, fixed **before any test scores
are viewed**; compare classification metrics on exactly the
same outer splits. Treat node counts and accuracy
as comparators; probabilities are optional and method-specific.

**Do not** merge either alternative into the currently executing
285-fit RLR/RLMT run, change its branch, or infer model rankings
from incompatible historical metric tables.

No extra baseline model has been fitted in this audit.
