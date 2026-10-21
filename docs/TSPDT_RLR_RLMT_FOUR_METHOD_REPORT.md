# Four-method comparison — matched 19 RLR-publication datasets, five folds

**Status: reporting code prepared; not executed on the server at
this commit.** The training results themselves were already audited:
`PUBLISHED_5FOLD_COMPLETE_AUDIT=PASS` (285 prior methods)
and `TSPDT_PUBLISHED_COMPLETE_AUDIT=PASS` (95 TSPDT-inspired
methods), as reported by the operator.

## Four fully separated method identities

1. `RLR` — original published relative logistic regression.
2. `RLMT_v1` — earlier two-leaf model-tree prototype.
3. `RLMT_DC_v0` — discovery/confirmation model tree.
4. `TSPDT2011_INSPIRED_RECON_V1` — newly implemented and
   explicitly **not historically source-equivalent** TSPDT
   reconstruction (outer-training-only MAD f100, max_depth3,
   min_leaf5, Wilson cf0.25 pruning, Laplace1 leaf scores).

The archived `BigTSP::tsp.tree` is a **different pair-feature
CART** that was not installed or run. It must not be mislabeled
as a fifth model in this report.

## Read-only reporting design

Script:
`scripts/tspdt_four_method_published_report.py`.

Inputs:

- 285 already-saved RLR/RLMT model records and
  `predictions_test.csv` from the existing
  `results/rlmt_dc_published_data_5fold_v1/` tree.
- 95 saved TSPDT-inspired model records and predictions from
  `results/tspdt_recon_19x5_v1/` in a separate worktree.
- Independent, previously completed RLR/RLMT five-fold audit,
  accepted-leaf replay, first-fold parity and summary, and
  independent TSPDT 95-model replay audit.
- Identical checksummed source inventories, exact paired
  `world / dataset / split_id / sample_id / y_true` across
  all four methods on every outer fold.

No `fit` calls, no Rscript, no new hyperparameters or
tree construction occur in this report.

Independently recompute from each saved probability/label
vector with one **common** definition:

- Accuracy: total correct / test-set count.
- Balanced accuracy: unweighted mean of two class recalls.
- Macro-F1: unweighted mean of two class F1 values.
- Logloss: held-out binary negative log-likelihood,
  numerical clipping only to [1e-15, 1-1e-15].
- Brier: held-out binary squared probability error.

All four models must have exactly the same test identities,
class labels and paired outcomes within an outer fold.
Compare the recomputed RLR/RLMT logloss and Brier to the
independent 285-run audit; compare TSPDT logloss/Brier,
accuracy, balanced accuracy and macro-F1 with the 95-run
TSPDT audit. Fail closed if any model or fold is missing,
mislabeled, incomplete or mismatched.

Expected counts:
`19 datasets × 5 outer folds × 4 models = 380 model-fold
metrics`, summarized in 95 paired fold comparisons and
19 equal-weight dataset-level means.

Outputs, **not created until actual successful server execution**:

`results/tspdt_recon_19x5_v1/four_method_comparison_v1/`

- `all_380_model_fold_metrics.csv`
- `matched_95_folds.csv`
- `matched_19_dataset_means.csv`
- `dataset_level_contrasts.json`
- `README.txt`

The script refuses to overwrite existing report output.

## Scientific interpretation

Report a dataset-level paired contrast (each method minus RLR),
its cross-dataset mean and median and wins/losses/ties; do not
mistake the 95 fold-level observations for 95 independent
biomedical cohorts. Report both successes and degradations.
Do not compare old manuscript Table 2 accuracy with
present results as if original study partition/preprocessing
were identical. Treat TSPDT's **Laplace leaf-score** logloss
as a reconstruction-specific secondary metric, not
probabilistic output of the 2011 original.

No subsequent tuning is authorized on the examined holdouts.
This unified descriptive report is exploratory and cannot
on its own establish superiority, independent validation,
or model equivalence.
