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
