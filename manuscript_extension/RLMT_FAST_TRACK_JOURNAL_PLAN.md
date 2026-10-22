# RAPID JOURNAL REVISION — RLR → RELATIONAL LOGISTIC MODEL TREES (RLMT)

Research program: RLR/RLMT, 2026-10-22. This is a new, unpublished
manuscript revision branch based on PAPER_FROZEN_MAIN
\`4079a5ac400f1d2cfded523c7000f655470925ce\`.
Never modify \`main\` or \`manuscript_lock/\` or overwrite benchmark outputs.

## Why this paper, rather than another RLR-only conference submission

The original agreed scientific aim (2026-10-20) was **a methodologically
new interpretable TREE algorithm**, suited to the user's tree-methods
habilitation: TSPDT-style within-patient relational internal splits
\`x_i>x_j\` and fit RLR logistic models as locally interpretable leaves.
A single leaf yields the RLR method family, but numerical parity with
the archived RLR version must be demonstrated under the exact fitting
protocol; current \`RLR_ALIGNED_V2\` has a separately named CV algorithm.

Keep the user's existing ten-page IEEE-style paper as the *base
scientific study*, not as a second simultaneous conference publication.
Target one coherent journal paper, with the added RLMT method central:
**Relational Logistic Model Trees for Interpretable Probabilistic
Omics Classification** (PROVISIONAL title; revise only once empirical
claims are verified).

The original PDF and associated 20-page supplementary document are
available as LaTeX at:
- \`manuscript_lock/RLR3_submission_final.tex\`
- \`manuscript_lock/Supplement_RLR_submission_final.tex\`
These files are locked; DO NOT change or overwrite them.

## Model taxonomy, mandatory controls and order

1. **Original published-form RLR**, exactly as in the existing paper:
   mandatory preserved reference; K=10 plus K=20/50 original results.
2. **RLR_ALIGNED_V2 (single-leaf root)**: distinct exact code/label,
   explicit stratified 5-fold \`glmnet\` CV; a fair new aligned-parent
   control, never mislabeled as original RLR without verified parity.
3. **TSPDT-style relational routing + RLR leaves**, simple Gini/
   class-impurity routing *control*; same leaf learner and budget.
4. **RLMT-R**: TSPDT-style relation candidate routed by RLR
   parent-residual criterion, train-only independent D/V confirmation;
   current H1 (aligned estimator procedure) is an experimental
   implementation/starting point, not proof of improvement.
5. **RLMT-LA**: important **original planned contribution**;
   nomination requires (computationally bounded) model-aware
   child-RLR look-ahead with FULL local fit and proper-score check,
   selected solely inside discovery; independent confirmation.
   Explicit top-pair candidate pre-screen followed by small
   child-look-ahead shortlist; freeze exact budgets from a
   synthetic-only complexity gate, not old test outcomes.
6. **RLMT-H2 offset**: a separately named *probabilistic ablation*,
   using aligned-parent OOF offsets and fixed-lambda ridge
   corrections, NOT mathematically identical to independently
   fitted RLR leaf. Keep as an ablation, not a semantic substitute
   for RLR-leaf model tree.
7. **Original RLMT-DC-v0 and RLMT-CA-PS-v1**: historical attempts/
   negative ablations in supplement; never silently cherry-pick them.

**Binary classification ONLY**; original agreed small-structure scope:
prepruned root-to-two-child leaves, at most 2 realized leaves
(bounded-depth 2 with depth convention explicitly stated), no
multiclass, no growing a large tree then post-hoc pruning.
Respect minimum total/class-specific child support.

## Mandatory fast gates in order

### P0 — Already completed; preserve, do not rerun

- Main original RLR paper baseline and 20-page Supplement.
- Frozen RLR/RLMT 19×5 original comparator 285 fits.
- RLMT-DC-v0 and CA-PS-v1 95-fold each; CA-PS 13/95
  naturally accepted, 4 improved and 9 worsened outer
  logloss in accepted cases; dataset-equal ΔLL vs RLR
  +0.006177806. **Negative/mixed**: do not relabel as victory.
- v2 Stage1: real H1/H2 synthetic SWAP/glmnet fits, forced
  test-only two-leaf RDS replay, test-label independence PASS;
  downloaded GEO GSE25055+GSE25065 source matrices
  310+198 GSM. No clinical outcome test has been run.

### P1 — Algorithm finalization (NOW)

- Wait for Stage2 H1/H2 synthetic regression and GEO
  schema/ID audit outputs; do not start clinical testing
  until clinical labels/group gate is cleared.
- Finish the **missing TSPDT+RLR / RLMT-LA controls**, not
  only H2: identical depth/leaf budget, comparable
  train-only relation selection, final consistent
  D/T estimator, deterministic CV, complete RDS/replay,
  tested accepted and rejected cases, data leakage checks.
- Synthetic-only methodological check: one-leaf identity,
  candidate freeze, no usage of held-out test labels,
  adverse probability tails, no split exact fallback.
- Pin one unchanging experiment SHA and a sample/metric
  manifest before any original biomedical benchmark.
- Scientific decision: accept negative/neutral outcomes.
  No retuning on old 19 result outliers.

### P2 — Fast, directly comparable *exploratory* testing

- First 2 cohorts × 3 identical publication folds; complete
  all planned methods/controls, detect errors, measure CPU.
- Then 19 original datasets ×5 identical folds (95 fold
  records per method) with frozen original RLR predictions
  imported READ ONLY and \`RLR_ALIGNED_V2\` separately fit;
  compare AUC, macro-F1, BA, accuracy, proper logloss,
  Brier, acceptance, leaf sparsity and min support.
  If serious harm/no leaf acceptance, STOP escalation
  and honestly report failures.
- **For paper-compatible primary results**: original
  CompCancer *11 datasets* with the EXACT manuscript's
  repetition schedule: 9 datasets × (10 repetitions ×5
  folds) + 2 datasets × (2 repetitions ×5 folds) = **470
  paired outer fold positions per method**. 19×5 results
  are a new *development screen*, not numerically identical
  to the manuscript's primary 470-fold analysis.
- Keep original PTB fixed-split and PTB LOCO as structural
  and cross-cohort contexts; don't recompute all paper
  figures until a credible RLMT result exists.
- Formal inference only at dataset level, with
  prespecified family of comparisons and appropriate
  correction; 470 correlated folds are not 470 independent
  scientific units. Separate exploratory from confirmatory
  analysis transparently.

### P3 — Manuscript preparation IN PARALLEL with P1/P2

- Preserve original Introduction and methods as the
  baseline; rewrite contribution as a new conditional
  probabilistic TREE algorithm instead of claiming
  unequal RLR relation weights are new by themselves.
- Add Methods: relational route \`x_i>x_j\`, model-aware
  nomination on discovery D, confirm on independent V,
  root-to-two-leaf execution, S-only supervised relation
  selection, RLR leaf logistic fitting, no-split fallback
  and the D/T estimator alignment protocol.
- Add explicit equations, algorithm pseudo-code,
  honest comparison with original 2011 TSPDT and 2008
  Weight-kTSP by Czajkowski & Kretowski.
  Historical Weight-kTSP is NOT the same as cost-sensitive
  WTSP of other authors; avoid conflation.
- Add figure: 1 RLR-only tree, 2 leaves with relation
  split and gene-pair logistic terms; separate optional
  nested H2 offset ablation.
- Main table: original RLR vs aligned root, TSPDT+RLR,
  RLMT-R, RLMT-LA and (if applicable) H2; paired
  proper score + structural footprint. Table 2: per-
  dataset wins/losses/ties and accepted/safe fallback.
- Discussion: probabilities can worsen despite no
  change in 0.5-threshold accuracy. Directly report
  the negative RLMT-DC/CA-PS findings; do not
  present old exposure as new heldout discovery.
- Update abstract, conclusion and source/code
  availability only after actual results complete.
- For the journal format, migrate frozen IEEEtran
  structure to \`elsarticle\` in the derived new
  manuscript (not the locked original); check
  Computer Methods and Programs in Biomedicine
  current official author instructions before sending.

### P4 — Publication go/no-go

**Journal RLMT story is substantiated** only when a
methodologically real TSPDT-like tree with RLR
leaves and explicit look-ahead has passed code gates,
is reproducible on frozen benchmark folds, and
adds something scientifically interesting: reliable
probability/accuracy gains, an important structural/
interpretability trade-off, or a carefully evidenced
negative design lesson, transparently quantified.

If RLMT remains predominantly identical to RLR
or consistently worsens proper score, do not delay
indefinitely or fabricate an improvement. Submit
the already evidence-backed **RLR-only original paper**
after its separate final presentation/author review,
while treating RLMT as a distinct later study.

## Immediate writer work

The supplied 10-page original and 20-page supplement
were reviewed. Start with **the new METHODS block now**,
then placeholders for new figures/tables and
preserve a clear \`PENDING RESULTS\` boundary. Do NOT
alter numbers/tables for CompCancer, PTB or legacy8
without a regenerated source-verified analysis.
