# CRIT-v7 — isolated **synthetic-only** repair after first failed launch

## Failure evidence from operator Ubuntu stdout

First attempt ran from pinned commit
\`6989cf2f2a1f4af75018c74fba5f901304115839\`
in \`/opt/2026/RLR/relative-logistic-regression-crit-v7\`.
Code and Python-auditor parsing passed. The first real
synthetic test started and terminated with:

    from glmnet C++ code (error code 7777);
    All used predictors have zero variance
    CRIT7_PILOT_BATCH_EXIT=1

The FIRST attempt therefore did **NOT** execute any
biomedical CRIT7 models. This was an incorrect synthetic
test fixture, **NOT a negative biomedical benchmark**.

## Exact root cause

The synthetic R test generated only **10 measured genes**.
The original publication uses \`m=10\` for the within-patient
rank-difference coding:

    delta = rank(gene_b)-rank(gene_a)
    code  = if delta>=10 then +1
            else if delta<=-10 then -1
            else 0

For ten genes ranks range from 1 through 10; the
maximum absolute difference is 9. Consequently **all**
RLR predictor codes equal 0 for every synthetic patient.
The unchanged original \`glmnet\` fitting code properly
refused its degenerate all-zero design.

## Fix, frozen boundaries

- Change synthetic generator \`n_genes\` from 10 to 32;
  XOR-planted relationship continues to use genes 1-4.
  Keep original RLR K10 / m10 / d2 / cv.glmnet
  completely unchanged; no ad hoc synthetic m change.
- In the separate positive-control synthetic challenge,
  use 480 D / 240 V artificial patients and the COMPLETE
  496-gene-pair candidate set, so both planted pairs are
  available rather than randomly omitted by a 256-pair
  random proposal. This exception affects only the
  SYNTHETIC test: real experiment always retains
  its original \`proposal_n=256\`, up to six gates,
  and the exact same class and model constraints.
- Replace the inappropriate reuse of a synthetic
  parent fitted on a different sample in the
  positive-control gate check with a parent fit on
  *its own synthetic discovery D*.
- Add a fail-closed preflight which checks Git blob
  content identity for the original 8 biomedical
  algorithm/runner/auditor/dependency files against
  the first pinned attempt (SHA 6989cf2...). Their
  contents MUST NOT be modified by this repair.
- Add a SSH-session-independent launcher with
  \`start_new_session=True\`, fully redirected
  streams and separated PID+logs; refuse any results
  directory, PID or log that already exists.

**No prior output directories are modified.**
Use a NEW git worktree
\`/opt/2026/RLR/relative-logistic-regression-crit-v7-fix1\`.
The old worktree remains as auditable failed execution.

## Validation state

As of preparation: commit sources updated; **not yet
executed in R on the operator's Ubuntu server**.
Next run must pass synthetic ORIGINAL_RLR parity, and
only then six exposed real original outer folds and
24 method-fold independent prediction replays. In case
of any other error, status MUST remain nonzero and
no benchmark result claim is made. R source parse and
true synthetic + archived patient-level parity are
checked by the guarded batch after launch.

Scientific designation remains:
\`SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY\`.
