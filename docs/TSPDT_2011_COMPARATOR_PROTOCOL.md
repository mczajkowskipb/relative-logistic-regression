# TSPDT: historical baseline for the RLR / RLMT / RLMT-DC program

**Status: research comparator registered, implementation identity NOT yet
verified. No runs authorized on biomedical datasets.**

**Isolation:** this document is deliberately committed to
`feature/tspdt-baseline-evidence`, not the branch executing the
pre-registered 19-dataset × 5-fold, 285-fit `RLMT-DC` comparison.
Do not switch branches, update running R files, cancel, restart,
or run post-hoc TSPDT fits in
`results/rlmt_dc_published_data_5fold_v1`.

## Historical primary evidence

Marcin Czajkowski and Marek Krętowski (2011),
*Top Scoring Pair Decision Tree for Gene Expression Data Analysis*,
Advances in Experimental Medicine and Biology 696, 27–35.
DOI [10.1007/978-1-4419-7046-6_3](https://doi.org/10.1007/978-1-4419-7046-6_3).

Accessible author-provided chapter:
https://www.aragorn.wi.pb.edu.pl/~mczajk/AEMB2011.pdf

- §2.2, p. 29–30: greedy, recursive top-down trees whose
  node tests are newly selected TSP, k-TSP or weighted TSP-family
  relations on the current node's observations. The **root test
  coincides with the analogous original TSP method**. Subsequent
  nodes recompute their own rules; the tree is simplified by
  C4.5-like pessimistic pruning.
- §3, p. 30: **four separately reported variants**:
  TSPDT, k-TSPDT, Weight-TSPDT and Weight-k-TSPDT.
  Historical maximum k=9, not a hyperparameter we may choose
  based on our new RLR holdout scores.
- §3, p. 30–31: 11 **Kent Ridge microarray** studies, with
  original fixed train/test partitions for certain cohorts and
  10 repetitions of tenfold CV for cohorts without predefined
  partition. Feature selection is stated as **Relief-F, 10
  neighbors, 1000 attributes**, prior to classification.
  The historical wording does **not** itself establish that
  preprocessing was performed within each CV training fold;
  our prospective baseline MUST do all supervised Relief-F
  selection on each training portion only.
- Table 2, p. 31: comparison primarily uses **classification
  accuracy** and an algorithm-specific model-size count.
  It is **not** a table of probabilistic logloss or Brier scores.
  Historical Table 2 cannot be directly combined with new
  5-fold RLR results, because experimental partitions,
  preprocessing and sample representations may differ.
- The **GTSPDT** (2013) uses global evolutionary search and
  is a **different algorithm**. It is not a substitute for
  this 2011 greedy TSPDT baseline.

## Comparator selection — before inspecting results

The main priority is the original **TSPDT (one pair per
internal node)**, rather than GTSPDT or a generic CART trained
on precomputed gene pairs. Secondary comparator:
**k-TSPDT (maximum k=9)**, if its precise node voting,
tie-breaking and pruning behavior can be independently verified.
Weighted variants are optional additional *separately labeled*
methods, not interchangeable with unweighted TSPDT.

**Model identity gates:**

1. Retrieve the original 2011 implementation/source from PI
   archives, or identify a verifiable implementation that
   demonstrably matches node-local TSP pair ranking, tie-breaking,
   classification routing and pessimistic pruning. Capture
   license and SHA256 of the implementation. The paper alone
   lacks enough detail to infer all source-level decisions,
   pruning constants and tie rules. Do NOT call an
   approximate reimplementation `TSPDT2011_REPRODUCED`.
2. One candidate external implementation is R **BigTSP
   1.0** (Xiaolin Yang and H. Liu, 2012), function
   `tsp.tree`, with documented citation to the 2011
   chapter and possible class-probability output.
   Historic archive:
   https://cran.r-project.org/src/contrib/Archive/BigTSP/
   Package license: **GPL >= 2**; it was removed from CRAN
   in 2022. No direct code-level equivalence to the original
   paper has been established. Source and dependencies must
   be reviewed *before* deciding whether to use it.
   If semantics differ, label the comparator `BigTSP_tsp_tree`,
   **not** original TSPDT.
3. A generic `rpart` or `tree` classifier trained on a
   static matrix of globally pre-ranked comparison features
   is only a `pair-feature CART` comparator, not TSPDT
   unless it provably implements node-wise fresh TSP ranking
   and original prune semantics.
4. Check 1-TSP versus k-TSP node orientation, repeated
   genes, ties, zero-variance expression profiles,
   stopping criteria, pruning, default leaf decision
   and lack of full source certainty. Provide unit tests
   on synthetic data and *only training-partition* behavior.
5. Never use held-out outcomes of known Gordon, Singh,
   leukemia or breast2 to set depth, prune penalty,
   candidate screening or other parameters.

## Prospective apples-to-apples experimental contract

**Current ongoing benchmark remains unchanged:**
`RLR / RLMT-v1 / RLMT-DC-v0`, 19 historical datasets,
5 original outer folds each, **285 total** fits.

When that benchmark has finished AND passed
`PUBLISHED_5FOLD_COMPLETE_AUDIT=PASS`, we may run a
**separate frozen TSPDT comparison**:

- Use exactly the **same 95** outer train/test partitions
  and the same checksummed 19 original ARFF files.
  Both previously seen improvements and degradations must
  be included; no 'best-looking' dataset subset.
- Baseline target **TSPDT-1** (95 fits).
  Secondary, if reference identity is resolved:
  **k-TSPDT(k<=9)** (another 95 fits).
  Both add **0 fits to the currently running job**.
- Make input feature selection strictly train-only.
  Explicitly label `ReliefF-1000` as a reconstruction of
  the published pipeline; ensure its training folds are not
  leaked into test. Do not silently re-use the **RLR global
  supervised pair top1000 pool** as TSPDT's node search space:
  those are algorithmically different selections.
- **Primary metric for historically comparable classifier
  performance: accuracy / balanced accuracy**, with class
  counts, per-class recall and macro-F1. Also report tree
  internal nodes, gene pairs and leaves. Keep the RLR/other
  methods in the same metric table on matching holdouts.
- **Probabilistic comparison:** RLR/RLMT logloss and Brier
  are well-defined from their fitted probabilities.
  TSPDT hard class labels alone **do not yield calibrated
  probabilities**. Define a separately identified,
  train-only leaf-probability model (e.g., prespecified
  Laplace smoothing) if comparing logloss/Brier, and **do
  not** portray it as original 2011 TSPDT.
  Its smoothing hyperparameter must be frozen before testing,
  not selected based on outer holdout scores.
- Report paired within-dataset 5-fold differences,
  then summarize 19 dataset-level means. Five folds
  are correlated; datasets are not independent newly
  recruited cohorts. Avoid a claim of prospectively
  independent validation or cross-method significance
  based on these already-used inputs.
- Reuse the completed 285-fit RLR/RLMT output as a
  **read-only comparator**; no rerun or overwrite of
  trained models is required if fold identity/prediction
  audits pass.

## Immediate noncomputational tasks

- Resolve original TSPDT source and executable/license status.
- Audit the BigTSP `tsp.tree` source archive only as a
  *candidate*, not a definitive reproduction.
- Check overlap between 2011 study IDs and the
  19 RLR-paper ARFFs; identical labels are not evidence
  of identical sample IDs or original study preprocessing.
- Create an immutable TSPDT implementation/specification
  with tests and freeze a new output root only after
  model identity and preprocessing are accepted.

**No TSPDT performance on the 19 new paired five-fold splits
has been computed or claimed in this checkpoint.**
