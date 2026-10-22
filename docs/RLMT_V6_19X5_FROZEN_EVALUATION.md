# RLMT-v6 — frozen 19×5 original-RLR-aligned retrospective experiment

**SCIENTIFIC STAGE: EXPOSED RETROSPECTIVE DEVELOPMENT ONLY.**
**STATUS: prepared source/runner/auditor; not executed on Ubuntu.**
The complete original-paper RLR results are retained READ ONLY.

## Completed original rank-margin v6 pilot (operator log 2026-10-23)

The 2-cohort ×3 frozen original folds ×5 model modes run:
- \`RLMT6_PILOT_BATCH_EXIT=0\`;
- \`RLMT6_ORIGINAL_RLR_REAL_SYNTHETIC_PARITY=PASS\`;
- \`RLMT6_30_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS\`;
- \`RLMT6_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS\`;
- \`ORIGINAL_RLR_UNMODIFIED=TRUE\`.

| Mode | 6-fold mean paired Δ heldout LL vs exact original RLR | Fold wins/losses/ties vs original |
|---|---:|---:|
| ORIGINAL_RLR | +0.000000000 | 0/0/6 |
| GLOBAL_CAL | -0.000489915 | 5/1/0 |
| **ONEPAIR_MARGIN** | **-0.003264461** | **6/0/0** |
| CONSENSUS_MARGIN | -0.000462556 | 5/1/0 |
| ALLPAIR_MARGIN | -0.002772675 | 4/2/0 |

ONEPAIR_MARGIN mean ΔLL vs GLOBAL_CAL:
**−0.002774546** (paired mean difference).
In contrast to v4/v5, these use NO new random gene-pair
candidate search: only the frozen original RLR K<=10 pairs.
This is still **highly exposed development**; a favorable
six-fold result must not be interpreted as proof.

Note: the pilot log field \`used_orig_pair_count=59/6\`
for CONSENSUS_MARGIN is a sum of existing selected pair
features across six folds, NOT 59 binary split events.
The 19×5 auditor explicitly separates
\`folds_with_rank_margin\` from
\`total_original_pair_feature_uses\`.

## Exact unchanged algorithm freeze

Original pilot algorithm commit:
\`a57c86efbff99fc8dbf8025b6db1e9f4fc6aa510\`.
The new benchmark's preflight compares Git blob hashes of:
- \`R/experimental/rlmt_v6_original_rank_margin.R\`,
- \`R/experimental/rlmt_v4_rank_repair.R\`,
- its same real synthetic test,
- original RLR utilities / parent dependencies.

A changed algorithm aborts before any real biomedical run.
**No changes to gamma=.15, global beta=.02, maximum
rank-margin |z|=3, log-odds cap ±2, supervised 1-of-K
selection, OOF seed, feature pool K10, or parent fitting.**
Only orchestration, exact original cohort matching, saving
results, and independent audit were added.

All 19 original publication cohorts are matched using
the exact source SHA256, same original seed42 5-fold patient
ID split, same 285-model original archive RLR heldout IDs
and labels. Each outer fold fits one source-parent model and
the fixed v6 alternatives. Original root predictions must
match archived RLR within 1e-9 or abort.

- 19 cohorts ×5 folds=95 source training bundles
- 5 modes per fold = 475 full method-fold prediction
  vectors, of which one is an exact original RLR control
  and four are proposed correction models.
- Outputs under the **new** immutable result directory
  \`results/rlmt_v6_original_rank_margin_19x5_v1\`.
- Same real synthetic parity+affine rank invariance before
  real fits, and independent per-patient post-run logloss/
  Brier/BA/MacroF1 reproduction from saved source rows.

## Primary analysis

Equal-cohort mean paired heldout **logloss change**
against BOTH original RLR and GLOBAL_CAL, separately
for each frozen predeclared mode. Lower ΔLL better.
We report all 19 dataset results, win/loss/tie counts,
Brier score, balanced accuracy, macro-F1, leave-one-
dataset-out sensitivity and selected original gene-pair
repeatability/stability across five folds.
Do not select/omit results based on favorable cohort
membership. Summaries are descriptive, not externally
confirmatory p-values.

**Key falsification checks:** a single bad dataset should
not determine the entire global mean improvement;
real heldout patient probability predictions should be
finite and matched to original, even for zero-correction
and high-confidence RLR cases; ONEPAIR should not receive
a gain solely by performing a different global calibration.

## Important causal/study limitations

All nineteen datasets and fold scores have previously
affected model development. Even if v6 19×5 yields an
improvement, it remains exploratory and **does not
establish generalization to unseen populations**.

The shortlist for up to 10 selected original RLR gene
pairs is supervised on outer training. Although parent
probabilities for fitting the rank-margin correction
are generated using strictly OOF original-parent fits,
the chosen full-training shortlist and its category
means/residual features are not themselves completely
reselected in each OOF fold. This *can* overfit correction
fitting. Outer-test labels are never used to train or
choose an individual fold's model, but prospective
validation is still needed.

For publication: external untouched cohort and
prospectively fixed model/hyperparameters, or a nested
out-of-source design with train-only model selection,
are preferred to another favorable exposed dataset
screen. Do not change the already-ready ORIGINAL_RLR
paper on the strength of this retrospective trial alone.
