#!/usr/bin/env bash
# RLMT-CA-PS-v1 LARGE SYNTHETIC GATE v3.
# 0 published datasets; test-only forced 2-leaf real fit plus six real
# simulations and independent Python replay. Preserve all old gate roots.
set -Eeuo pipefail
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
HEAD="$(git rev-parse HEAD)"
test "$(git rev-parse origin/main)" = "$MAIN" || {
 echo "STOP: frozen publication main changed" >&2; exit 2;
}
test "$HEAD" = "$(git rev-parse origin/feature/rlmt-class-aware-ps-v1)" || {
 echo "STOP: checkout differs from latest pinned experimental branch" >&2
 exit 2
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "STOP: tracked worktree modifications" >&2; exit 2;
}
command -v Rscript >/dev/null || {
 echo "STOP: Rscript unavailable" >&2; exit 2;
}
command -v python3 >/dev/null || {
 echo "STOP: Python3 unavailable" >&2; exit 2;
}
PRIOR="results/rlmt_ca_ps_v1_synthetic_gate_v2/status.txt"
test -f "$PRIOR" && grep -qx 'RLMT_CA_PS_DEV_EXIT=0' "$PRIOR" || {
 echo "STOP: prior real synthetic v2 gate has no confirmed success" >&2
 exit 2
}
OUT="results/rlmt_ca_ps_v1_synthetic_gate_v3"
test ! -e "$OUT" || {
 echo "STOP: v3 evidence directory exists; refusing overwrite" >&2
 exit 2
}
mkdir -p "$OUT"
printf "%s\n" "$HEAD" > "$OUT/CODE_SHA.txt"
printf "%s\n" "$MAIN" > "$OUT/PUBLICATION_MAIN_SHA.txt"
on_exit() {
 rc=$?
 printf 'RLMT_CA_PS_LARGE_GATE_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== A. PARSE ALL REAL AND SIMULATION SOURCES ==="
Rscript --vanilla -e '
 files <- c(
  "R/experimental/rlmt_ca_ps_v1.R",
  "R/experimental/test_rlmt_ca_ps_v1.R",
  "R/experimental/test_rlmt_ca_ps_v1_real_synthetic.R",
  "R/experimental/test_rlmt_ca_ps_v1_real_forced_split.R",
  "R/experimental/rlmt_ca_ps_v1_six_case_synthetic_evaluation.R")
 invisible(lapply(files,parse))
' > "$OUT/parse.log" 2>&1 || {
 cat "$OUT/parse.log"; exit 1;
}
python3 -m py_compile scripts/rlmt_ca_ps_v1_synthetic_suite_audit.py
echo "RLMT_CA_PS_LARGE_R_PARSE=PASS"

echo "=== B. ORIGINAL SYNTHETIC UNIT + MOCK GATE ==="
Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1.R \
 > "$OUT/unit.log" 2>&1 || {
 cat "$OUT/unit.log"; exit 1;
}
grep -qx 'RLMT_CA_PS_SYNTHETIC_UNIT_AND_MOCK_INTEGRATION=PASS' "$OUT/unit.log"
echo "RLMT_CA_PS_ORIGINAL_UNIT_GATE=PASS"

echo "=== C. REAL ORIGINAL RLR/SWAP/GLMNET — UNFORCED ==="
Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1_real_synthetic.R \
 > "$OUT/real_unforced.log" 2>&1 || {
 cat "$OUT/real_unforced.log"; exit 1;
}
cat "$OUT/real_unforced.log"
grep -qx 'RLMT_CA_PS_REAL_SYNTHETIC=PASS' "$OUT/real_unforced.log"
echo "RLMT_CA_PS_REAL_UNFORCED=PASS"

echo "=== D. REAL FITTED TWO-LEAF BRANCH (FORCED TEST GATE ONLY) ==="
Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1_real_forced_split.R \
 > "$OUT/real_forced_two_leaf.log" 2>&1 || {
 cat "$OUT/real_forced_two_leaf.log"; exit 1;
}
cat "$OUT/real_forced_two_leaf.log"
grep -qx 'RLMT_CA_PS_REAL_FORCED_ACCEPTED_BRANCH=PASS' "$OUT/real_forced_two_leaf.log"
grep -qx 'FORCED_CONFIRMATION_BOOL_ONLY=TRUE' "$OUT/real_forced_two_leaf.log"
grep -qx 'REAL_PARENT_AND_TWO_LEAF_RDS_REPLAY=PASS' "$OUT/real_forced_two_leaf.log"
grep -qx 'HOLDOUT_LABEL_PERMUTATION_INVARIANCE=PASS' "$OUT/real_forced_two_leaf.log"

echo "=== E. SIX FIXED REAL TRAINING/TEST SYNTHETIC EXPERIMENTS ==="
Rscript --vanilla \
 R/experimental/rlmt_ca_ps_v1_six_case_synthetic_evaluation.R \
 "$ROOT/$OUT/simulations" \
 > "$OUT/simulation.log" 2>&1 || {
 cat "$OUT/simulation.log"; exit 1;
}
cat "$OUT/simulation.log"
grep -qx 'RLMT_CA_PS_SIX_CASE_SYNTHETIC_EVALUATION=PASS' "$OUT/simulation.log"

echo "=== F. INDEPENDENT SAVED-SAMPLE METRIC AUDIT ==="
python3 scripts/rlmt_ca_ps_v1_synthetic_suite_audit.py \
 --root "$ROOT/$OUT/simulations" \
 > "$OUT/independent_replay.log" 2>&1 || {
 cat "$OUT/independent_replay.log"; exit 1;
}
cat "$OUT/independent_replay.log"
grep -qx 'RLMT_CA_PS_SYNTHETIC_INDEPENDENT_METRIC_REPLAY=PASS' "$OUT/independent_replay.log"

echo "=== ALL LARGE SYNTHETIC GATES COMPLETE ==="
echo "RLMT_CA_PS_REAL_FORCED_TWO_LEAF_GATE=PASS"
echo "RLMT_CA_PS_LARGE_SYNTHETIC_GATE=PASS"
echo "REPORT=$ROOT/$OUT"
