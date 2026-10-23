# CRIT-v7 PRE-FLIGHT FIX 02: operator-observed systemd immediate EXIT=1

On 2026-10-23, the user correctly enabled:
\`sudo loginctl enable-linger "$(whoami)"\` (Linger=yes).

Then \`systemd-run --user\` successfully launched a transient
\`crit-v7-fix1-pilot.service\`, but this exited within milliseconds
with status 1. There were NO results, no status, and no
\`crit_v7_systemd.log\` output. This is a failed *early preflight*,
not an interrupted background job or a biomedical training outcome.

**Root cause identified in the frozen fix1 shell runner:**

    HEAD="$(git rev-parse HEAD)"
    test "$(git rev-parse origin/feature/crit-v7-regime-switch-probe)" = "$HEAD"

The fix1 worktree HEAD was a0aa4cb8... while
\`origin/feature/crit-v7-regime-switch-probe\` remained the
original 6989cf2... branch. The equality is unconditionally false.
Because \`set -e\` is enabled, the runner exits immediately,
before creating its result tree, log messages or status trap.

**Repair:** new isolated branch
\`feature/crit-v7-preflight-fix-02\`, new user worktree
\`relative-logistic-regression-crit-v7-fix2\`, and **new**
\`systemd-run --user\` unit \`crit-v7-fix2-pilot\`.

Only \`scripts/crit_v7_2x3_guarded_batch.sh\` preflight
was changed to check its OWN new branch and print explicit
fail-closed messages on mismatch. The optional detached
launcher has its \`BRANCH\` identity updated for consistency.

The main CRIT-v7 R algorithm, ORIGINAL RLR root wrappers,
original publication reference, synthetic m=10 repair,
all four controlled methods, archived biomedical
pilot runner, python independent audit, heldout patient
splits, and scientific status are all unchanged.

**Operator must run a fresh isolated worktree**, never
rerun over fix1 or remove its logs. Fetch this branch
from GitHub before using detached HEAD so origin branch
identity and frozen source match.

Do not claim any biomedical or synthetic PASS until
\`CRIT7_ALL_GATES=PASS\`,
\`CRIT7_PILOT_BATCH_EXIT=0\`, and independent archived
patient-level control parity are actually observed.

Lingering user services under \`systemd --user\` survive
normal SSH disconnection. A transient service is NOT
configured to restart after reboot. Do not promise
reboot survival.
