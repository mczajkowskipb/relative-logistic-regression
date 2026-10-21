# RLMT-CA-PS-v1 — observed complete 19×5 biomedical retrospective checkpoint (2026-10-21)

**Evidence source:** Operator-provided Ubuntu console output from the
full immutable source-matched study; actual per-fold archives reside
under the separate Ubuntu worktree rather than inside the Git branch.
The measurements are **retrospective/exposed development evidence**.

## Completed verification, not planned future results

Frozen experimental source **\`e049fe20530ae7c5378db73f4e9186b5f83f789e\`**.
On Ubuntu, preflight and full model run returned:

\`\`\`
RLMT_CA_PS_RETRO_PREFLIGHT=PASS
ORIGINAL_REFERENCE_FITS=285
PUBLISHED_INPUT_SHA256_VERIFIED=19
OUTER_FOLDS_MATCHED=95
NEW_MODEL_FITS_PLANNED=95
RETRO_CA_MODELS=95
RETRO_CA_ACCEPTED=13/95
RLMT_CA_PS_RETRO_19X5_RUN=PASS
RLMT_CA_PS_RETRO_INDEPENDENT_19X5_AUDIT=PASS
RETRO_FOUR_METHOD_FOLD_MODELS=380
RETRO_FOUR_METHOD_DATASETS=19
RLMT_CA_PS_RETRO_19X5_COMPLETE=PASS
RLMT_CA_PS_95FOLD_BATCH_COMPLETE=PASS
\`\`\`

All 95 CA-PS models fitted on all 19 exact old biomedical
cohorts x the original 5 stratified held-out folds. SHA256 checks of
all input files and paired sample-level source/RLR predictor parity
passed. The 285 original model archives were read, not retrained
as comparators. All original result roots remain untouched.

**Report directory on Ubuntu:** \`results/rlmt_ca_ps_v1_retrospective_19x5/four_method_retro_audit_v1\`.

## Published 19-dataset-equal-weight descriptive results vs original RLR

| Model | Accepted splits | Accuracy dataset wins/losses/ties | Dataset-mean delta heldout logloss |
|---|---:|---:|---:|
| RLMT-v1 | not a comparable confirmation gate | 0/3/16 | +0.024251703 |
| RLMT-DC-v0 | 9/95 | 2/0/17 | +0.017217934 |
| **RLMT-CA-PS-v1** | **13/95** | **1/1/17** | **+0.006177806** |

**Interpretation:** CA-PS decreases the adverse **excess logloss**
of DC-v0 by approximately 64.1%, but its mean is still
**worse than original RLR** by +0.006177806.
This does NOT establish superiority, calibration improvement,
statistical significance, or generalization to unobserved
biomedical cohorts. The 19 datasets were involved in
method development, and 95 folds are not 95 independent
scientific studies.

The run emitted the R message \`There were 50 or more warnings\`.
It reached PASS, but this warning class/count is **not
independently triaged** from the logs. Do not imply that
absence of fatal errors means absence of model warnings.

## Frozen methodological distinction

- v0 DC: split nomination on training discovery D,
  1SE confirmation on independent train-only V,
  refit on full outer train T; local leaf model
  may generate overconfident wrong probabilities.
- v1 CA-PS: class-support-aware nomination score,
  original min class/leaf support and parent-shrunk
  probability mixture with \`B=4n0n1/(n0+n1)^2\`,
  \`g=B*n/(n+20)\`. The original D->T estimator
  alignment problem is **not completely solved**.
- The shrink weight (and split screen multiplier)
  is a heuristic proposed after exposure to old
  datasets; **do not optimize it on these 95 outcomes**.

## Next posthoc analysis — NEW read-only code, result not yet observed

\`scripts/rlmt_ca_ps_v1_accepted_13_posthoc_audit.py\` performs
a separate read-only analysis of **all 13 accepted trees
and both fitted leaves per tree (26 total)**, extracting:

1. Paired heldout CA vs RLR and old DC logloss,
   Brier/accuracy/balanced-accuracy deltas and
   actual harmful/beneficial accepted-fold counts.
2. Final two-leaf train class0/class1 support,
   node balance indices, true support-based shrink weights.
3. Correct independently reconstructed per-sample
   final probability mixtures, including zero-test-
   support leaves without inventing observations.
4. Both leaves' contributions to complete paired
   outer-fold excess logloss; number of very
   confident wrong-class predictions.
5. Best/worst accepted five folds, dataset means
   and reasons for nonacceptance on all 95 folds.

No training and no threshold changes. Results written
only under a fresh
\`results/rlmt_ca_ps_v1_retrospective_19x5/accepted_13_posthoc_diagnostic_v1\`
subdirectory, NEVER overwriting the independent baseline
audit, fitting manifest or original 285+95 model artifacts.

**Scientific decision until diagnostic is run:**
freeze RLMT-CA-PS-v1, preserve its mixed/negative
retrospective results, and prioritize a new
estimator-aligned and stable probabilistic method
with a *separate, genuinely unseen cohort* for
prospective evaluation. Avoid choosing changes
post-hoc from breast2/leukemia/prostate fold tails.
