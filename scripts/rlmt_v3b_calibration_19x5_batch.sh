#!/usr/bin/env bash
# SAFE new 19x5 calibration-controlled exploratory research batch.
# Single shared fit per original old-paper fold, 4 probability outcomes.
# Run in a NEW detached worktree through nohup, no old model writes.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "V3B_STOP: provide immutable old RLR reference absolute dir" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
SHA="$(git rev-parse HEAD)"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
REF_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-ridge-stump-calibration-control)" = "$SHA"
test "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
test -f "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for cmd in Rscript python3 sha256sum; do
 command -v "$cmd" >/dev/null
done
OUT="results/rlmt_v3b_calibration_19x5_v1"
test ! -e "$OUT" || {
 echo "V3B_STOP: previous output exists; never delete or overwrite" >&2
 exit 2
}
echo "=== 1. REAL CODE SYNTAX PREFLIGHT ==="
Rscript --vanilla -e '
 parse("R/experimental/rlmt_v3_stumps_pruned.R")
 parse("R/experimental/rlmt_v3b_calibration_control.R")
 parse("R/experimental/test_rlmt_v3b_synthetic.R")
 parse("R/experimental/rlmt_v3b_calibration_19x5_runner.R")
' >/dev/null
python3 -m py_compile scripts/rlmt_v3b_19x5_independent_audit.py
echo "V3B_SYNTAX=PASS"
mkdir -p "$OUT"
printf '%s\n' "$SHA" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REF_SHA" > "$OUT/IMMUTABLE_ORIGINAL_SHA.txt"
printf '%s\n' \
 "EXPOSED retrospective 19x5, four modes, shared 1 RLR parent and OOF, route min 8, 4096 candidates, ridge lambda 1, no tune by outer test" \
 > "$OUT/FROZEN_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'RLMT_V3B_19X5_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. SYNTHETIC REAL FOUR-MODEL FIT AND RDS REPLAY ==="
Rscript --vanilla \
 R/experimental/test_rlmt_v3b_synthetic.R \
 > "$OUT/synthetic_gate.log" 2>&1 || {
 cat "$OUT/synthetic_gate.log"; exit 1;
}
tail -n 40 "$OUT/synthetic_gate.log"
grep -qx "RLMT_V3B_SYNTHETIC_REAL_FITS=PASS" "$OUT/synthetic_gate.log"

echo "=== 3. NINETEEN EXPOSED COHORTS TIMES FIVE ORIGINAL FOLDS ==="
Rscript --vanilla \
 R/experimental/rlmt_v3b_calibration_19x5_runner.R "$REF" \
 > "$OUT/380_method_folds.log" 2>&1 || {
 tail -n 100 "$OUT/380_method_folds.log"; exit 1;
}
tail -n 35 "$OUT/380_method_folds.log"
grep -qx "RLMT_V3B_19X5_RUN=PASS" "$OUT/380_method_folds.log"

echo "=== 4. COMPLETELY INDEPENDENT PATIENT-LEVEL METRIC REPLAY ==="
python3 scripts/rlmt_v3b_19x5_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT_V3B_19X5_INDEPENDENT_PATIENT_AUDIT=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT_V3B_19X5_CONTROLLED_STUMPS_ALL_PASS"
echo "STUDY_IS_PREVIOUSLY_EXPOSED_RETROSPECTIVE=TRUE"
echo "ORIGINAL_RLR_UNTOUCHED=TRUE"
echo "REPORT=$ROOT/$OUT"
