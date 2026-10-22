# RLMT-v2 — 19×5 FOUR-METHOD LARGE RETROSPECTIVE SCREEN

Implementation/checkpoint: 2026-10-22. Distinct branch
\`feature/rlmt-v2-19x5-all-methods\`. This is a **predeclared
computational study**; no results are claimed until the operator
executes and the independent patient-level audit reports PASS.

## Goal and protocol

The next large execution compares all four **already implemented**
algorithm versions without changing design constants after examining
their previous 6-fold pilot:

- H1: residual-route nomination, independent RLR leaves, aligned root;
- H2: fixed-ridge local logit-offset correction, not genuine independent
  RLR leaves; retained as explicit probabilistic ablation;
- TSPDT_GINI: Gini-selected relational route, then actual two-leaf RLR;
- RLMT_LA: top-20 parent-residual-screened candidate relations with
  genuinely fitted RLR leaves in independent inner OOF folds,
  model-aware proper-score nomination, then D/V confirmation.

Each method freezes exactly the **same** algorithm and split-selection
settings from the pilot; user-exposed pilot LL results do not change
per-class support, proposal budget, ridge lambda or confirmation cutoff.

- Dataset inventory: 19 original exact SHA-verified cohorts:
  9 CompCancer, 2 CompCancer distinct repetition setting, legacy8.
- Outer folds: the *original* first repetition, stratified
  5fold \`seed=42\`. 95 paired outer fold positions.
- Planned 380 NEW fit attempts, each with own artifact directory.
- Input SHA256 checks repeated per cohort; all 95 original outer
  training/test sample sets cross-checked with immutable RLR archives.
- Original RLR 285-published-model reference suite is strictly
  read-only. No reruns of original models.
- Per-model: save RDS and test predictions, acceptance, logloss,
  Brier, aligned and original root logloss, class support,
  training/runtime warnings (including \`glmnet\` class <8), and
  complete failure messages. Metrics are regenerated independently
  from the saved per-patient predictions.
- Independent audit recomputes 380 expected paired model/fold rows,
  exact fallback identity and method-aligned root parity; returns
  dataset-equal 19-cohort means and 19-dataset win/loss/tie counts.
- If any one fit fails, the runner still writes a \`FIT_ERROR\` record
  and preserves all previously successful fits. The final job fails;
  it CANNOT claim 380 clean fits or launch independent PASS audits.

## Critical scientific and computing constraints

This is a **retrospective exposed development screen**. Old training
and test cohorts have already informed experimental RLMT design.
No superiority claims, no unbiased external validation, no
posthoc parameter tuning. There are 95 related fold positions, NOT 95
independent cohorts. Primary external validation requires a
genuinely independent source with frozen data/label/grouping protocol.

Known numerical hazard: glmnet binomial warning
\`one multinomial or binomial class has fewer than 8 observations;
dangerous ground\`. It was observed during previous nested
Gini/LA tests. Do NOT suppress or treat it as inconsequential;
current batch stores every captured warning message and
per-model count in \`warning_inventory.csv\`.
Console repetition may be suppressed AFTER being recorded in the
inventory to avoid enormous logs; no model or fold is hidden.
Model success alone is not scientific validity.

**Performance caveat:** RLMT_LA may require up to 20 candidates
times 3 nested folds times 2 leaf fits per original fold.
For 95 large gene-expression folds this can be CPU-expensive.
Only start with stable disk, adequate quota, and a clear plan to
leave a long-lived \`nohup\` process running. Never start a duplicate
batch, nor delete a partial output root to force a rerun.

**Journal gate:** A robust original published RLR manuscript is
already complete. If this four-method benchmark shows mostly
single-leaf fallbacks and near-zero or adverse proper-score
changes, STOP RLMT augmentation as a claim of performance
improvement and proceed with author-reviewed RLR publication.
The model-tree branch remains a separate methodological
investigation. Never omit the negative results or alter
the original archived paper.

**Paper-number distinction:** 19×5 = 95 exposures per method
is **NOT** the manuscript's primary CompCancer 11 dataset
9×10×5 + 2×2×5 = 470 paired outer fold protocol. Even if
95 passes, a 470-fold paper-aligned run (plus statistical
reporting by dataset, not pseudoreplicated folds) is a
separate, downstream stage, not interchangeable evidence.
