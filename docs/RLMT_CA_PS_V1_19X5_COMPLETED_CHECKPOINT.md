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

## Independent posthoc accepted-tree audit — actually completed on Ubuntu

The operator ran the new read-only audit successfully after the full 95-fold run,
with **`CA_PS_POSTHOC_REPLAY=PASS`**, **13/95 accepted trees** and **26
fitted leaves**. Validation compared saved per-sample mixed probabilities,
original RLR parents, old DC results, final leaf support and weights and
recomputed every accepted fold's paired test loss without refitting.

Observed precise outcomes on the **13 accepted outer test folds**:

- **Benefit/harm/tie by unweighted logloss: 4/9/0.** Thus confirmed
  inner-train gains do not reliably predict final outer test benefits.
- **One** CA-PS true-class probability below 0.01 on accepted-test routes;
  the shrink rule reduces extreme confident errors but does not prevent all.
- Dataset-equal weighted Δlogloss CA-PS vs RLR: **+0.006177806**;
  original RLMT-DC-v0: **+0.017217934**. Both are adverse vs RLR.
- Worst accepted CA-PS vs RLR (test-set mean ΔLL):
  - breast2 rep00_fold04 **+0.430246**, wrong CA/RLR **7/7**;
    old DC **+0.000000**;
  - breast2 rep00_fold01 **+0.157154**, wrong **9/7**;
    old DC **+1.698460**;
  - chen-2002 rep00_fold00 **+0.066263**, wrong **3/2**;
    old DC **+0.000000**;
  - prostate rep00_fold02 **+0.065793**, wrong **0/0**;
    old DC **+0.009130**;
  - chen-2002 rep00_fold04 **+0.060272**, wrong **2/0**;
    old DC **+0.000000**.
- Four favorable accepted folds:
  - leukemia rep00_fold01 **−0.102860**; original DC **−0.150278**;
  - breast2 rep00_fold00 **−0.076087**; original DC **+0.126315**;
  - leukemia rep00_fold00 **−0.036364**; original DC **−0.060396**;
  - gordon-2002 rep00_fold02 **−0.013010**; original DC **+0.000000**.

**Interpretation:** Acc gains and probability quality diverge. The
breast2 fold04 counterexample worsens logloss substantially despite no
additional hard-label classification errors. Prostate fold02 worsens
logloss although both RLR and CA made **zero** classification errors.
These are descriptive *observations*, not evidence of disease cohort
shift or an identified causal mechanism for DC/CA instability.

**Freeze/close RLMT-CA-PS-v1:** Do not retroactively change split
balance factors, shrinkage strength, confirmation cutoff or rerank on
the exposed fold outcomes. No publication-grade generalization or
improved probabilistic modeling claim is supported by these data.
The code can remain as a research ablation with fully recorded 95-fit
artifacts, 380-model independent comparison and 13-tree posthoc
reconstruction. The next scientific question must be separately posed
and tested on new independent cohorts.

The completed posthoc audit is stored only on operator Ubuntu:
`results/rlmt_ca_ps_v1_retrospective_19x5/accepted_13_posthoc_diagnostic_v1/`.
This repository documentation transcribes console evidence; it does
not claim to contain raw, locally stored per-sample records.
