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


## Executed Stage 2 + pilot checkpoint — verified Ubuntu outputs (2026-10-22)

User provided actual full terminal logs after executing both isolated
background tasks on the Ubuntu system (no long-running processes remain).

**First real-paper exploratory 2×3 matched fold pilot (12 fits)**
- \`RLMT_V2_PILOT_EXIT=0\`, \`RLMT_V2_PILOT_COMPLETE=PASS\`,
  \`RLMT_V2_PILOT_INDEPENDENT_AUDIT=PASS\`.
- H1 residual-screen RLR-leaf: 1/6 accepted, 1/0/5
  better/worse/tie vs **RLR_ALIGNED** on heldout logloss;
  unweighted mean fold ΔLL = **−0.0001580**.
- H2 fixed-ridge logit-offset ablation: 6/6 accepted,
  5/1/0 better/worse/tie vs **RLR_ALIGNED**; mean fold
  ΔLL **−0.0016300**. These two numbers are small and
  observed on only TWO exposed cohorts, three outer
  folds each. NO significance or generalization claim.
- **Different control identity:** the original published
  \`RLR\` and experimental \`RLR_ALIGNED_V2\` cannot be
  treated as mathematically the same fitted model;
  comparisons must show both.

**New tree-selector synthetic model fits:**
- \`RLMT_LA_SYNTHETIC_GATE_EXIT=0\`,
  \`RLMT_V2_GINI_LA_SYNTHETIC_REAL_GATE=PASS\`,
  \`RLMT_LA_SYNTHETIC_BATCH=PASS\`.
- Both TSPDT_GINI and RLMT_LA completed real R/glmnet
  two-leaf RDS replay and separate outer-label
  prediction insensitivity.
- Both NATURAL synthetic tests fell back to 1 leaf:
  two-leaf verification was intentionally FORCED
  in a synthetic TEST harness. Do not call this a
  naturally accepted scientific split.
- The actual R runtime emitted >=50 warnings in the
  H1/H2 biological pilot and 15 warnings in parts of
  LA synthetic fitting. Their **types are not yet
  independently triaged**; successful exit does not
  imply warning-free fits.

**Next fixed experiment**, under independent branch
\`feature/rlmt-paper-fourway-pilot\`:
- Preserve the completed 12 H1/H2 models, train ONLY 12
  Gini/LA models on same 2 exposed CompCancer cohorts ×
  same 3 original folds; then independently recompute
  24 fold-method metrics with exactly matched original
  RLR and common aligned parent.
- No ability to retrospectively change the already
  viewed splits or design settings. If this pilot fails
  or remains close to 1-leaf, report it and do not
  silently add favorable cohorts.

**Manuscript boundary:** The derived
\`manuscript_extension/RLMT_journal_draft.tex\`
is a writing draft with frozen old RLR figures/tables
and an imported proposed RLMT Methods section.
No v2 results tables or conclusion sentences are
ready until all four methods and the original RLR
are actually compared.

The full 470-fold article primary design remains
9 datasets × 10 repetitions × 5 folds +
2 datasets × 2 repetitions × 5 folds.
The earlier 19 datasets × 5-fold retrospective
development screen is NOT interchangeable with
the publication's 470-fold benchmark.


## Four-method matched real biomedical pilot — operator status (2026-10-22)

The new two-dataset, three-original-fold paired experiment
completed successfully on Ubuntu. User-supplied console
evidence:

\`\`\`
RLMT_V2_FOURWAY_PILOT_EXIT=0
RLMT_V2_FOURWAY_PAPER_PILOT_COMPLETE=PASS
SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY
\`\`\`

All 12 additional TSPDT_GINI and RLMT_LA model-fitting
records were produced and the runner's independent audit
was executed before return code 0. Prior H1/H2 results were
reused, not trained again.

The observed last **five** fit records were no-split RLR_ALIGNED
fallbacks (delta logloss vs aligned parent exactly zero).
This subset alone does **not** establish the acceptance total
for either new method across all six splits. The independent
fourway summary log at
\`results/rlmt_v2_fourway_pilot_2x3_v1/independent_audit.log\`
contains the definitive per-method totals but those numbers
were NOT included in the operator's pasted excerpt.
Do not invent them in the paper.

Multiple \`glmnet::lognet\` warning headers occurred.
The snapshot from \`fourway_background.log\` only shows the
header: the full, actual warning **message bodies** reside in
\`results/rlmt_v2_fourway_pilot_2x3_v1/real_biomedical.log\`.
Warning category, frequency and implications for deviance
and probability convergence must be evaluated before
expanding to full publication folds.

**Next research decision** is contingent on reading the
already-computed four-method audit in full and the first
12 complete warning messages. If LA/Gini have few/no
accepted splits and H2 is only a minor proper-score
improvement, do not represent this as a superior
relational decision tree. Preserve negative controls and
prioritize credible journal claims and the existing
strong RLR-only fallback.


## VERIFIED four-method paper-fold pilot results — 2026-10-22

User supplied the complete independent audit and every new Gini/LA
per-fold result (6 exact paired original CompCancer folds:
\`chen-2002\` and \`gordon-2002\`, folds 00–02).
The independent metrics audit passed: **24/24 method-fold** records;
run exit **0**.

All Δlogloss values below are versus a separately fit
\`RLR_ALIGNED\`, **NOT** the frozen original RLR.

| Method | Accepted / 6 | Fold LL win/loss/tie | Mean per-fold Δlogloss |
|---|---:|---:|---:|
| H1 residual RLR leaves | 1 | 1/0/5 | −0.000158015 |
| H2 fixed ridge log-odds offset | 6 | 5/1/0 | −0.001629974 |
| TSPDT_GINI / RLR leaves | 1 | 1/0/5 | −0.000158015 |
| RLMT_LA / RLR look-ahead leaves | 0 | 0/0/6 | +0.000000000 |

For Gini, the *only* accepted fold was
\`gordon-2002/rep00_fold00\`, Δlogloss −0.0009481;
all other Gini and all six LA folds fell back exactly
to the aligned global parent. The numerical equality
of H1 and Gini means does NOT demonstrate identical
models unless accepted fold identities/probabilities
are compared.

**Warning root identified:** genuine \`glmnet::lognet\`
messages repeatedly state \`one multinomial or binomial
class has fewer than 8 observations; dangerous ground\`.
These appear during nested, small, imbalanced node CV
fits. Successful program termination and RDS replay do
not establish stable leaf coefficient estimates or
reliable probabilities. Report class-support diagnostics
and method-specific warning counts. Do not merely
suppress this warning or retrospectively relax the
independent confirmation gate to force split acceptance.

**Publication impact:** H2 showed a tiny positive
proper-score signal in only two EXPOSED development
cohorts and is an offset model, NOT a conventional
two-leaf RLR tree. The proposed central RLMT_LA
algorithm has currently **zero natural accepted splits**.
The expanded RLMT-centered paper is therefore
**NOT ready for submission** on this evidence.
Do not present these numbers as generalization,
statistical significance, or better than original RLR.

Next: independently report delta vs ORIGINAL archived
RLR from the same pilot CSV; implement all-19 matched
full-fold model evaluations with measured warning
counts and no output rewriting. Keep clean original
RLR manuscript as an immediate fallback route to
submission rather than adding weak RLMT claims.


## GO / NO-GO after full 19×5 retrospective experiment — 2026-10-22

**FINAL ACTUAL OPERATOR-SUPPLIED RESULTS, NOT SPECULATION:**
Experiment \`feature/rlmt-v2-19x5-all-methods\` completed all
380/380 actual fits. Shell status:
\`RLMT_V2_19X5_BATCH_EXIT=0\`;
\`RLMT_V2_19X5_PATIENT_LEVEL_PREDICTION_AUDIT=PASS\`;
\`RLMT_V2_19X5_ALL_380_NEW_MODELS_AUDITED=PASS\`.

19 formerly examined biomedical datasets; each 5 original
stratified outer folds, total 95 paired outer fold positions per
model. All means are **equal-dataset means** of held-out
logloss, never pooled-as-independent fold statistics:

| Model | Accepted/95 | Dataset W/L/T vs RLR_ALIGNED | mean ΔLL vs RLR_ALIGNED | mean ΔLL vs ORIGINAL archived RLR | glmnet <8-class warnings |
|---|---:|---:|---:|---:|---:|
| H1 RLR-leaf | 9 | 1 / 5 / 13 | +0.013088521 | +0.014324918 | 618 |
| H2 fixed-ridge offset | 25 | 5 / 3 / 11 | −0.000244199 | +0.000992198 | 496 |
| TSPDT_GINI with RLR leaves | 7 | 3 / 1 / 15 | −0.002324349 | −0.001087952 | 562 |
| RLMT_LA lookahead RLR leaves | **0** | 0 / 0 / 19 | +0.000000000 | +0.001236397 | 2533 |

Total captured small-class warnings: **4209**, messages, not
distinct failure counts. All fits have status OK, but warning
volume remains scientific concern. Test examples include
\`prostate rep00_fold02\`: H1 accepted but worsened heldout
LL +0.0113941 relative to aligned parent; Gini accepted and
improved LL −0.0050370. Holdout outcomes were not permitted
to influence inner selection.

**Important baseline accounting:** \`RLR_ALIGNED\` averages
+0.001236397 LL relative to archived ORIGINAL RLR.
Thus a method that is tied with the aligned root is
slightly worse than published RLR in this retrospective
benchmark. Do not call aligned root the original RLR.

**Scientific interpretation and stop/go:**
- RLMT_LA is *not* an empirically supported improvement:
  natural splitting 0/95. Stop optimizing it on exposed
  retrospective fold outcomes and do not conduct a 470-fold
  main-paper experiment merely to chase a favorable result.
- H1 is distinctly detrimental (+0.0143 ΔLL vs original).
  Do not integrate into a superiority claim.
- H2 is a small offset ablation, not an RLR-leaf model tree,
  and is worse than the original RLR (+0.00099).
- TSPDT_GINI is the only experimental method with average
  directionally positive LL vs original RLR (−0.00109),
  but with 7/95 accepted folds and 3/1/15 dataset
  W/L/T relative to aligned parent; its improvement could
  be dominated by a dataset or a few route decisions.
  Its stability, outlier contribution, clinical utility and
  replication have NOT been independently established.
  No statistical-significance, publication-quality
  superiority or novel RLMT-LA evidence follows.

**NEXT FIXED ACTION:**
1. Read \`independent_19x5_audit/method_19_dataset_means.csv\`
   and \`independent_19x5_audit/four_method_summary.csv\`
   with an **independent read-only** dataset-contribution
   diagnostic. Show positive/negative contribution and
   leave-one-dataset-out robustness for Gini, compare
   both baseline definitions. Preserve all adverse folds,
   no parameter tuning.
2. Prioritize completing **original RLR manuscript**,
   reproducing original core tables and supplementary
   results, checking venue format and statistics, and
   preparing truthful submission. Keep RLMT separate
   until new independently frozen cohort or methodological
   evidence justifies its own article.
3. If publication authors require 470-fold evaluation
   of a fixed Gini tree model despite weak exploratory
   evidence, use a prospectively frozen protocol,
   budget estimates and no posthoc screening; clearly
   label reused CompCancer datasets as exposed, not
   external validation. A new independent
   cohort/locked validation design is preferable.

**NO REWRITE of original paper, no altering archived
results, no algorithm changes, no rerunning 380 fits.**
