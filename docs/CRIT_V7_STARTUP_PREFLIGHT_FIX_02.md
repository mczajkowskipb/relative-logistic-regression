# CRIT-v7 attempt02 — silent preflight failure repair

Observed 2026-10-23 from operator's Ubuntu terminal:

- previous fix1 attempt created `pilot.pid` (PID 205114),
  and empty `crit_v7_background.log` at 14:06;
- that process is no longer active;
- `results/crit_v7_original_rlr_2x3_v1/status.txt`,
  synthetic log and biomedical fold log are absent;
- consequently the fix1 experiment did NOT enter the R tests.

Root cause found in the guarded shell batch: following a legitimate
branch update from `feature/crit-v7-regime-switch-probe` to
`feature/crit-v7-m10-synthetic-repair-01`, the batch's
preflight still compared `HEAD` with the former branch reference.
It rejected the new commit before emitting ANY output or setting
the EXIT status trap. This is operator-independent and is **not**
evidence that the SSH detached process was killed by logout.

Fix:
1. pin the batch branch reference and launcher to
   `feature/crit-v7-m10-startup-fix-02`;
2. print `CRIT7_BATCH_PROCESS_ENTERED=PASS` as the first output
   and install an ERR trap before the Git preflight, so the next
   preflight failure will have a legible error;
3. use separate `crit_v7_attempt02_background.log` and
   `pilot_attempt02.pid` in the SAME existing isolated worktree;
   preserve fix1 logs/PID and NEVER delete previous output;
4. preserve original bootstrap commit SHA checks, SHA of
   `origin/main`, original RLR pinned SHA, 8 exact original
   CRIT algorithm/dependency Git blob hashes, and all original
   archived-patient parity requirements;
5. keep the original synthetic fix (32 genes, m=10) from fix1;
   do not change biomedical model code or the six-fold design.

Status: SOURCE PATCH APPLIED on GitHub; the operator must run the
new attempt on Ubuntu to establish R synthetic/biomedical results.
Scientific stage: EXPOSED_RETROSPECTIVE_DEV_ONLY.

If an unexpected results directory is present, STOP and inspect it.
Do NOT override files and do NOT claim PASS unless the guarded batch
and independent patient-level audit return success.
