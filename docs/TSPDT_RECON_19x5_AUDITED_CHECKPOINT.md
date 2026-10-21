# Matched TSPDT-inspired 19×5 — actual 95-fit checkpoint

**Recorded 2026-10-21 from the PI's compute-server transcript.**
Original RLR publication `main` remains
`4079a5ac400f1d2cfded523c7000f655470925ce`.
TSPDT model code was separately proven to pass the synthetic gate on
`aa4dae241ba583394501f2f11efd757c8d74ad46`.

## Completed run — observed, not speculative

The existing TSPDT output root was present when the launcher was
called again; its overwrite guard stopped with
`TSPDT_LAUNCH_STOP: refusing to overwrite existing independent TSPDT result root`.
**This is intentional no-overwrite behavior, not evidence of a
failed fitting algorithm.** The original launch log is not included
in this checkpoint.

The user then ran the read-only original-model replay auditor in
`/opt/2026/RLR/relative-logistic-regression-tspdt-dev`.
That produced:

```
AUDITED_TSPDT_FOLDS=95
TSPDT_95FOLD_AUDIT=PASS
TSPDT_PUBLISHED_COMPLETE_AUDIT=PASS
```

The auditor required an existing, successful full TSPDT run,
all 95 fit rows with `status=OK`, a 285-fit successful RLR/RLMT
reference manifest, original sample/label matching, serialized
`model.rds` replay, agreement of test probabilities and class
decisions and recomputation of per-fold metrics. The TSPDT metrics
remain **reconstructed**, not historical 2011 scores.

## Complete descriptive dataset-level findings

The audit compared 19 original publication datasets, each averaging
five matched held-out folds. The numbers below are **TSPDT-inspired
minus RLR**; positive accuracy/F1 is better, negative is worse.
The unit of these summaries is the **dataset**, not any individual fold.

| Metric | Datasets improved | Datasets worsened | Mean dataset-level difference | Median difference |
| --- | ---: | ---: | ---: | ---: |
| Accuracy | 1/19 | 18/19 | -0.08573 | -0.10857 |
| Balanced accuracy | 3/19 | 16/19 | -0.08375 | -0.11000 |
| Macro-F1 | 2/19 | 17/19 | -0.09566 | -0.11942 |

Mean accuracy difference represents approximately **-8.57 percentage
points**; balanced accuracy **-8.37 p.p.**, macro-F1 **-9.57 p.p.**
These aggregate values have been recalculated from the **rounded
19 per-dataset differences supplied in the console transcript**;
slight last-digit variation is expected compared with raw CSV.

The main accuracy outlier is `cns`, for which TSPDT-inspired
accuracy was `0.718531` against RLR `0.554196`
(positive difference `+0.164336`). Other 18 datasets had
lower accuracy than RLR.

Balanced accuracy improved on only
`cns` (+0.182500), `adenocarcinoma` (+0.029487),
and `colon` (+0.007500).
Macro-F1 improved on `cns` (+0.187841) and
`adenocarcinoma` (+0.052780).

These are not estimates of external-domain generalization; all
19 cohorts were included in the already exposed publication/method
development data. There is no justification to retune
TSPDT-inspired parameters or to select CNS as the only
subsequent benchmark.

## Scope and interpretation

- **Strong observed conclusion:** The *specified reconstruction*
  (training-only MAD f100, node-local TSP ranking, depth3,
  min-leaf5, Wilson cf0.25) underperforms the historical RLR
  baseline on the majority of these matched datasets.
- **Invalid extrapolation:** This does not imply that the
  original 2011 TSPDT algorithm performs identically, because
  its exact source is unavailable and its reported experiments
  used other feature selection and split protocols.
- **BigTSP:** A distinct pair-feature CART implementation,
  not installed on the compute host and not run as part of these
  95 experiments. It must not be represented as tested.
- **Important limitation:** The audit verifies saved model replay
  and same-split RLR comparison; for a complete four-method
  head-to-head table we should **read** saved RLR, RLMT-v1,
  RLMT-DC and TSPDT outputs and independently recompute a
  common accuracy, balanced accuracy and macro-F1 definition.
  Existing 285-fit and 95-fit results must **never be retrained**
  just to generate that table.

**Decision:** Keep TSPDT-inspired v1 fixed, preserve its negative
comparison as evidence, and do not optimize parameters against the
already examined holdout folds. Prepare common read-only reporting
across four algorithms instead. To revisit an original TSPDT source
reproduction would require recovering more implementation detail or
explicitly preregistering another historical reconstruction.


## Four-way common-metric audit completed — observed PASS (2026-10-21)

The operator executed the corrected, read-only four-model report in
the *separate TSPDT worktree*. All 380 model-fold records, 95
pair-matched folds and 19 dataset-level summaries passed.

```
FOUR_METHOD_CLASS_ORDER_REGRESSION=PASS
TSPDT_STRUCTURAL_POSTMORTEM=PASS
RLR_RLMT_TSPDT_FOUR_METHOD_REPORT=PASS
COMPARISON_MODELS=380
COMPARISON_PAIRED_OUTER_FOLDS=95
COMPARISON_DATASETS=19
```

### Descriptive dataset-level comparisons to original RLR

W/L/T refers to gains/losses/exact ties among **19 five-fold
dataset means**. Positive logloss difference means **worse**.

| Method | Accuracy W/L/T | Balanced accuracy W/L/T | Mean delta logloss (method minus RLR) |
| --- | --- | --- | ---: |
| RLMT-v1 | 0/3/16 | 1/3/15 | +0.024252 |
| RLMT-DC-v0 | 2/0/17 | 2/0/17 | +0.017218 |
| TSPDT-inspired v1 | 1/18/0 | 3/16/0 | +0.038045 |

**Science-critical qualification:** RLMT-DC improves accuracy on
two datasets without dataset-level accuracy losses, but **mean
dataset-level logloss is worse**. Classification and probability
performance point in different directions. The aggregate does not
identify which datasets/folds or probability tails dominate the
logloss harm; these must be examined in the read-only
`matched_19_dataset_means.csv` and `matched_95_folds.csv`.
Do not call this superiority.

The previous four-way reporter's stop at west-2001 was a **report
class-order mismatch**, not a change in trained model outputs.
The regression test passed with original ARFF level order `1,0`.
The completed report identified one nonlexical level sequence:
`compcancer_2x2/west-2001: ER -> ER+`.
No model refits or tolerance relaxations were made.

### All 95 TSPDT-inspired v1 model structures audited

```
TSPDT_MODELS=95
ROOT_ONLY_TREES=0
MEAN_INTERNAL_NODES=2.2947
MEAN_LEAVES=3.2947
MAX_OBSERVED_DEPTH=3
MAX_SELECTED_GENES=100
TSPDT_STRUCTURAL_POSTMORTEM=PASS
```

These observations **falsify the simple hypothesis** that
underperformance stems from Wilson post-pruning almost every
model to a root-only leaf. They **do not** prove Wilson
pruning is harmless: trimming higher-order branches and
changing effective depth are still possible. The strict
`MAD100` rather than original `ReliefF1000`, extra root
child-support feasibility requirement and `max_depth=3`
remain large non-equivalences to TSPDT 2011.

### Freeze decision and next analysis

Preserve original TSPDT-inspired v1 code / all 95 model RDS
artifacts and the original publication RLR 285-run results.
Do not retrospectively tune the paper methods using the 95
already observed holdouts.

The next technical task is a **read-only** decomposition of
RLMT-DC's positive mean delta logloss into complete per-dataset
and per-fold contributions, explicitly listing harmful
accepted trees, unchanged-root ties and prediction uncertainty
for rare held-out patients. A separate, source-grounded TSPDT
v2 method specification can then be developed without claiming
the flawed v1 to be the original TSPDT algorithm.


## RLMT-DC accepted-split outcome decomposition — independently checked (2026-10-21)

Operator ran `scripts/rlmt_dc_5fold_readonly_outcome_diagnostic.py`
against previously audited 285 RLR/RLMT + 95 TSPDT models.
It finished `RLMT_DC_READONLY_OUTCOME_DIAG=PASS`. The resulting
**19-dataset mean** DC-minus-RLR test logloss was **+0.017217934**.
These are **dataset-level averages over five paired outer folds**:

| Dataset | Mean Δlogloss | Mean Δaccuracy | Accepted DC folds |
| --- | ---: | ---: | ---: |
| breast2 | +0.364955 | +0.03750 | 2 |
| armstrong-2002-v1 | +0.002365 | 0 | 2 |
| prostate | +0.001826 | 0 | 1 |
| gordon-2002 | +0.000130 | 0 | 2 |
| leukemia | -0.042135 | +0.02667 | 2 |
| other 14 datasets | 0 | 0 | 0 |

Out of 95 outer folds, **exactly 9** accepted two-leaf models:
six harmed external holdout logloss and three improved it.
The **only** accepted folds with observed gain in held-out
accuracy were leukemia folds 00/01 and breast2 fold00;
breast2 fold01 lost accuracy (−0.06250), but averaged across
five folds the dataset still showed +0.03750 accuracy.

**Severe failure:** `breast2/rep00_fold01` increased held-out
logloss **+1.698460** and Brier **+0.173929**, with external
accuracy **−0.06250**. Some wrongly classified patients received
probabilities within 0.0003 of 0 or 1. Example individual
increase of logloss **+6.270339**: true class 1,
RLR probability 0.162552, DC probability 0.000307.
These are **observed post-hoc diagnostic probabilities**,
not inputs to method selection.

The breast2 dataset-wide delta +0.364955, distributed equally
among the 19 dataset means, dominates the total +0.017218;
other dataset contributions partly compensate. This disproves
a hypothesis of uniform loss degradation across cohorts.

### Mechanistic source-code inspection (hypothesis, not yet proven root cause)

DC acceptance is based on an independent confirmation subset
within each **outer training fold**:
`rlmt_dc_confirm` in `R/new/40_rlmt_dc_utils.R`
requires confirmation `mean(parent_loss − child_loss)>0`
and `gain − SE(gain)>0` on at least 20 observations.
This is a **one-standard-error threshold**, not a formal
distribution-free generalization guarantee.

There is an important **procedure difference** between the
confirmation predictors and final outer-test predictors:
during confirmation, parent/child leaf models are trained
on discovery subset `D` with a supervised pair pool
constructed from `D`, while final accepted leaf models
are **refitted on full outer-training T** using the common
outer-training pair pool. Both are training-only and neither
uses outer-test labels, but this means positive confirmation
gain is not a direct evaluation of the identical fitted
final leaf models. High-dimensional small-leaf fitted
probabilities may become extreme on unseen patients.

This suggests examining the saved `dc/decision.json`
confirmation gains, standard errors, and training/leaf class
counts against their subsequent holdout deltas. A NEW,
**strictly read-only** script was added:
`scripts/rlmt_dc_5fold_confirmation_vs_holdout_audit.py`.
It has **not yet been executed** on the server.
It must not affect any frozen acceptance decision.

**Decision:** Preserve all nine accepted trees and all unfavorable
results. Do not tune confirmation threshold, prune probabilities,
alter the validation fold or perform `breast2`-specific adjustments
after looking at these holdout errors. Future candidate methods
require a distinct, preregistered V1/V2 specification and new
external validation before claiming gains.


## RLMT-DC: original 1-SE confirmation vs independent outer test — audited PASS (2026-10-21)

Operator executed `scripts/rlmt_dc_5fold_confirmation_vs_holdout_audit.py`
against the saved nine accepted models and matched four-method
reference. Verified output:

```
DC_CONFIRMATION_ACCEPTED=9/95
DC_ACCEPTED_BUT_HARMFUL_ON_OUTER_TEST=6/9
DC_ACCEPTED_AND_BENEFICIAL_ON_OUTER_TEST=3/9
DC_CONFIRMATION_VS_TEST_AUDIT=PASS
NO_RETRAINING_NO_THRESHOLD_CHANGE=TRUE
```

All nine confirmation gates originally passed
`nV>=20, gain>0, gain-SE>0`. The contrast below is the *original
training-only confirmation gain* (positive is favorable) against
outer-test **DC-minus-RLR loss** (positive is harmful).

| Dataset / fold | nV | Confirmation gain | 1-SE lower bound | Outer-test Δlogloss |
| --- | ---: | ---: | ---: | ---: |
| breast2 / 01 | 21 | +1.604980 | +0.766044 | **+1.698460** |
| breast2 / 00 | 21 | +0.398981 | +0.035286 | +0.126315 |
| prostate / 02 | 28 | +0.549441 | +0.283100 | +0.009130 |
| armstrong / 02 | 20 | +0.063650 | +0.022337 | +0.006691 |
| armstrong / 04 | 20 | +0.029328 | +0.009111 | +0.005133 |
| gordon / 04 | 49 | +0.021361 | +0.002446 | +0.001597 |
| gordon / 00 | 48 | +0.003619 | +0.001262 | -0.000948 |
| leukemia / 00 | 20 | +0.002497 | +0.000630 | -0.060396 |
| leukemia / 01 | 20 | +0.032551 | +0.011464 | -0.150278 |

**Interpretation constraint:** The severe breast2 fold01 was
*not* a marginal `gain-SE≈0` gate acceptance. It had a
substantial observed confirmation margin **+0.766044** with
only 21 confirmation observations and a large **SE=0.838936**.
This is a descriptive training-validation result; one SE is
neither a 95% confidence interval nor a bound on unseen-test
generalization. This does **not** by itself isolate whether
failure was caused by confirmation sample variance, leaf
model instability, refitting from discovery D to full T, or
feature-pool changes between confirmation and final fit.

Source review confirms `R/new/42_rlmt_dc_run.R` validates
models trained only on D and subsequently refits final
leaf models on all outer training T using T's own shared
supervised pair pool, with local child reranking. This is
proper train/test separation, but is **not the identical
fitted estimator tested during confirmation**.

New follow-on *read-only* file
`scripts/rlmt_dc_5fold_final_leaf_route_audit.py`
compares all nine accepted models' final leaf training
class counts, outer-test route/class distributions,
active leaf features and leaf-wise per-sample logloss
contributions. It has **not yet been executed** on Ubuntu.
No model is refitted and no score-derived threshold is chosen.

**Scientific decision remains unchanged:** keep all original
artifacts and accepted decisions immutable; do not 'repair'
breast2 or choose a stronger acceptance criterion based on
these exposed holdout scores. Before any novel method is
claimed, preregister a separate V2 and evaluate on never-used
external datasets.
