#!/usr/bin/env bash
# CRIT-v8 Stage D: preserved six-fold candidate V reconstruction + paired-null simulation.
# Safe background batch; never modify prior results or original RLR.
set -Eeuo pipefail
trap 'rc=$?; echo "STAGED_RUNTIME_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 3 || { echo "STAGED_STOP: ORIGINAL_RLR V7_FIX3 V8 required" >&2; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
V7="$(cd "$2" && pwd -P)"
V8="$(cd "$3" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
BRANCH="analysis/crit-v8-stageD-selection-null-simulation-v1"
ORIG_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
V7_SHA="e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03"
V8_SHA="435b87ae527db82575d4b5cf800725250b53935d"
STAGEC_SHA="c5708bb55feb35ebfcc8baec89ac2d32d90afacb"
echo "STAGED_SOURCE_SHA=$HEAD"
test "$(git rev-parse "refs/remotes/origin/$BRANCH")" = "$HEAD" || {
  echo "STAGED_STOP: local HEAD does not match remote pinned branch" >&2; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$ORIG_SHA"
test "$(git -C "$V7" rev-parse HEAD)" = "$V7_SHA"
test "$(git -C "$V8" rev-parse HEAD)" = "$V8_SHA"
for tree in "$ROOT" "$REF" "$V7" "$V8"; do
 test -z "$(git -C "$tree" -c core.filemode=false status --porcelain --untracked-files=no)" || {
  echo "STAGED_STOP: tracked modifications in $tree" >&2; exit 2;
 }
done
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx "CRIT7_PILOT_BATCH_EXIT=0" \
 "$V7/results/crit_v7_original_rlr_2x3_v1/status.txt"
grep -qx "CRIT8_PILOT_BATCH_EXIT=0" \
 "$V8/results/crit_v8_one_shot_2x3_v1/status.txt"
for path in R/experimental/rlmt_v4_rank_repair.R \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/experimental/crit_v8_single_nomination.R \
 R/new/07_world_protocol_utils.R \
 R/new/28_rlmt_utils.R \
 R/new/40_rlmt_dc_utils.R \
 R/new/42_rlmt_dc_run.R; do
 test "$(git hash-object "$path")" = \
     "$(git rev-parse "$STAGEC_SHA:$path")" || {
  echo "STAGED_STOP: frozen algorithm changed: $path" >&2; exit 2;
 }
done
for ds in chen-2002 gordon-2002; do
 for f in 00 01 02; do
  test -f "$V7/results/crit_v7_original_rlr_2x3_v1/models/$ds/rep00_fold$f/complete_crit7.rds"
 done
done
for exe in Rscript python3 sha256sum; do command -v "$exe" >/dev/null; done
Rscript --vanilla -e 'parse("R/experimental/crit_v8_stageD_reconstruct_V_gains.R")' >/dev/null
python3 - <<'PY'
from pathlib import Path
for p in ("scripts/crit_v8_stageD_selection_null.py",
          "scripts/crit_v8_stageD_independent_audit.py"):
    compile(Path(p).read_text(encoding="utf-8"),p,"exec")
print("STAGED_CODE_PARSE=PASS")
PY
OUT="results/crit_v8_stageD_null_simulation_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "STAGED_STOP: output root exists; will not overwrite" >&2; exit 2;
}
if [ "${STAGED_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "STAGED_PREFLIGHT_ONLY=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$ORIG_SHA" > "$OUT/ORIGINAL_SHA.txt"
printf '%s\n' "$V7_SHA" > "$OUT/V7_SHA.txt"
printf '%s\n' "$V8_SHA" > "$OUT/V8_SHA.txt"
printf '%s\n' "EXPOSED_RETROSPECTIVE_DEV_ONLY; rerun only D-parent models for sixfold V-patient candidate score reconstruction; matched to saved v7 RDS to <1e-7; no outer test labels used; 5000 paired zero-mean centered sign-flip draws per fold; descriptive sensitivity not calibrated FWER and not a p-value for CRIT" \
 > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'STAGED_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 1. RECONSTRUCT ACTUAL V PATIENT GAIN MATRICES ==="
Rscript --vanilla R/experimental/crit_v8_stageD_reconstruct_V_gains.R \
 "$REF" "$V7" > "$OUT/reconstruction.log" 2>&1 || {
 tail -n 120 "$OUT/reconstruction.log" >&2; exit 1;
}
cat "$OUT/reconstruction.log"
grep -qx "STAGED_ARCHIVED_V7_SCORE_RECONSTRUCTION=PASS" "$OUT/reconstruction.log"

echo "=== 2. 30000 PAIRED ZERO-MEAN SIGN FLIP SENSITIVITY SIMULATIONS ==="
python3 scripts/crit_v8_stageD_selection_null.py \
 > "$OUT/null_simulation.log" 2>&1 || {
 tail -n 100 "$OUT/null_simulation.log" >&2; exit 1;
}
cat "$OUT/null_simulation.log"
grep -qx "STAGED_NULL_30000_PAIRED_SIMULATION_DRAWS=PASS" "$OUT/null_simulation.log"

echo "=== 3. INDEPENDENT ARCHIVED AND SIMULATION AUDIT ==="
python3 scripts/crit_v8_stageD_independent_audit.py "$V7" \
 > "$OUT/independent_audit.log" 2>&1 || {
 tail -n 100 "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "STAGED_ARCHIVED_V_GAIN_METRIC_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGED_30000_NULL_DRAW_SUMMARY_AUDIT=PASS" \
 "$OUT/independent_audit.log"

echo "STAGED_ALL_GATES=PASS"
echo "STAGED_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "STAGED_ARCHIVED_RESULTS_UNMODIFIED=TRUE"
echo "STAGED_REPORT=$ROOT/$OUT"
