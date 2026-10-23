# CRIT-v7 fix2: FAIL-CLOSED branch guard and SSH-independent service

## Actual cause of failed fix1 launch (operator logs 2026-10-23)

\`/opt/2026/RLR/relative-logistic-regression-crit-v7-fix1/\`
showed:
- \`pilot.pid\` existed, PID 205114 not alive;
- \`crit_v7_background.log\` length 0 bytes;
- no \`results/crit_v7_original_rlr_2x3_v1\` directory, status,
  synthetic log, or biomarker results.

**Confirmed static bug:** fix1's guarded batch compared HEAD of
\`feature/crit-v7-m10-synthetic-repair-01\` against the remote SHA of
\`feature/crit-v7-regime-switch-probe\`, which is the prior version.
Both are immutable, different SHAs. The comparison unconditionally
failed BEFORE the first log line and BEFORE creating the output dir.
No evidence that SSH or logind terminated an actual long computation.

**fix2 changes**:
- Compare batch HEAD to its OWN
  \`origin/feature/crit-v7-preflight-and-systemd-fix2\` ref and print
  \`CRIT7_SOURCE_REF=PASS\` on successful match.
- On any later fail, print \`CRIT7_RUNNER_FAILURE_EXIT\` + shell line.
- Add \`scripts/crit_v7_launch_systemd.py\` so the batch runs
  as a managed \`systemd --user\` service, not an SSH shell child.
  Refuse to start if user manager unavailable or \`Linger=yes\`
  cannot be confirmed. Ask loginctl to enable the user's own linger
  if currently no, subject to polkit. Enabling linger allows
  user@UID.service to remain active after the final SSH logout.
- Journal is the first source of process/error output:
  \`journalctl --user -u crit7-fix2.service -n 80 --no-pager\`.
  The existing nested logs and nonzero exit status in result dir
  are retained by the original guarded batch.
- Preserve fix1 worktree/zero-byte log/PID and earlier initial
  failure from 10-gene all-zero RLR design.
- Keep REAL algorithm, fixed candidate budget, source-backed
  original RLR, biological test procedure and patient archive
  auditor **byte-identical** to commit
  \`6989cf2f2a1f4af75018c74fba5f901304115839\`, verified by
  runtime SHA blob gate. The 32-gene synthetic fix from fix1
  is retained, not modified.

## Execution

Create NEW worktree \`relative-logistic-regression-crit-v7-fix2\`.
Fetch this exact branch and pin its final SHA. Start with:

\`python3 scripts/crit_v7_launch_systemd.py /opt/2026/RLR/relative-logistic-regression\`

The call stops if the result directory exists, an instance of the
\`crit7-fix2.service\` is already recorded, source identity is wrong
or linger cannot be enabled. No file removal is required.

The six original exposed cohort folds, actual test time,
independent patient-level audit and any possible performance
results are **NOT YET OBSERVED**. Successful launch does not
imply successful model fitting or scientific superiority.
