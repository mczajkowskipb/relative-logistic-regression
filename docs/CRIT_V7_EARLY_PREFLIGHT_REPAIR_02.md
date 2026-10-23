# CRIT-v7 preflight repair 02 — early-exit without diagnostics

Observed on Ubuntu `obliczeniowy`, previous fix1 checkout:
- `crit_v7_background.log` existed with exactly **0 bytes**.
- `pilot.pid` contained PID `205114`, which no longer existed.
- `results/crit_v7_original_rlr_2x3_v1` did **not** exist.
- Therefore **no synthetic or patient-model training took place**.

Root cause in prior `scripts/crit_v7_2x3_guarded_batch.sh`:
the repair-01 source HEAD was checked against the original
`origin/feature/crit-v7-regime-switch-probe` branch SHA.
This must fail, because repair-01 was deliberately on a
different branch. The comparison was made BEFORE the first
`echo`, BEFORE results directory creation and status trap.
The child exited nonzero, leaving no error message in an
otherwise properly detached background log. This is a
**deterministic Git preflight bug**, not evidence that SSH
disconnect killed any computations.

Fix:
1. Pin the correct current repair branch for the batch guard:
   `origin/feature/crit-v7-preflight-repair-02`.
2. Print a preflight entry marker and trap `ERR` from the
   first script line, with line number, bash command, exit code.
3. Implement `CRIT7_PREFLIGHT_ONLY=1`: test exact main,
   original-reference SHA, the current remote branch SHA,
   byte-identical CRIT7 biomedical implementation/deps
   against `6989cf2...`, worktree cleanliness and
   required commands/old 285-fit status, then quit without
   touching any results directory.
4. The new detached Python launcher synchronously invokes
   this no-side-effect preflight and requires
   `CRIT7_PREFLIGHT_ONLY=PASS` before creating ANY log
   or PID. On error, it displays the exact test and line.
5. The actual child runs the full preflight independently
   a second time so no safety guard is bypassed.

No change to synthetic data correction from repair-01,
original 2011 RLR root, tree fitting, candidate budgets,
outer fold choices, experiment's methods or independent auditor.
The previous fix1 worktree and its zero-byte log are retained.

**State:** committed and source-checked; awaiting actual
Ubuntu execution. Full scientific claim remains
`EXPOSED_RETROSPECTIVE_DEV_ONLY`. This version only prevents
a process launch bug; it does not prove any CRIT-v7 performance.
