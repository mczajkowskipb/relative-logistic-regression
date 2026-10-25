#!/usr/bin/env bash
# CRIT-v7 Stage B: 6 original folds x 6 methods, from pinned preserved fix3.
# Runs as one transaction without modifying original or fix3 output trees.
set -Eeuo pipefail
trap 'rc=$?; echo "STAGEB_EARLY_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 2 || { echo "STAGEB_STOP: ORIGINAL_REF FIX3_REF needed" >&2; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
ORIGINAL="$(cd "$1" && pwd -P)"
FIX3="$(cd "$2" && pwd -P)"
cd "$ROOT"
BASE_SHA="e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03"
PAPER_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
HEAD="$(git rev-parse HEAD)"
echo "STAGEB_RUNNER_STARTED=$HEAD"

test "$(git rev-parse origin/analysis/crit-v7-stageB-mechanistic-v1)" = "$HEAD" || {
 echo "STAGEB_STOP: local HEAD not pinned to stageB branch" >&2; exit 2;
}
test "$(git -C "$ORIGINAL" rev-parse HEAD)" = "$PAPER_SHA"
test "$(git -C "$FIX3" rev-parse HEAD)" = "$BASE_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$ORIGINAL" -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$FIX3" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$ORIGINAL/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx "CRIT7_PILOT_BATCH_EXIT=0" \
 "$FIX3/results/crit_v7_original_rlr_2x3_v1/status.txt"

# Frozen RLR + CRIT science source must remain byte-identical to fix3.
for path in \
  R/experimental/crit_v7_original_rlr_two_leaf.R \
  R/experimental/rlmt_v4_rank_repair.R \
  R/experimental/crit_v7_original_rlr_2x3_pilot.R \
  scripts/crit_v7_2x3_independent_audit.py \
  R/new/07_world_protocol_utils.R \
  R/new/28_rlmt_utils.R \
  R/new/40_rlmt_dc_utils.R \
  R/new/42_rlmt_dc_run.R
do
 test "$(git hash-object "$path")" = "$(git rev-parse "$BASE_SHA:$path")" || {
  echo "STAGEB_STOP: modified scientific dependency $path" >&2; exit 2;
 }
done
echo "STAGEB_FROZEN_SOURCES=PASS"

for bin in Rscript python3; do command -v "$bin" >/dev/null; done
Rscript --vanilla -e 'parse("R/experimental/crit_v7_stageB_mechanistic.R")' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_v7_stageB_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("STAGEB_CODE_PARSE=PASS")
PY

OUT="results/crit_v7_stageB_mechanistic_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "STAGEB_STOP: output root exists; refusing overwrite" >&2; exit 2;
}
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$BASE_SHA" > "$OUT/FIX3_SOURCE_SHA.txt"
printf '%s\n' "$PAPER_SHA" > "$OUT/ORIGINAL_RLR_SHA.txt"
printf '%s\n' "EXPOSED_RETROSPECTIVE_DEV_ONLY; mechanism ablation; original RLR; prior GLOBAL_INTERACTION; prior CRIT_TREE; same-gate GATE_ONLY; same-gate ADDITIVE_PAIRS; same-features UNCONDITIONAL_TREE; algebraically identical explicit global GLM; all new corrections fitted using strictly outer-training 3fold RLR OOF; no new selection or model hyperparameter tuning from outer test" > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf "STAGEB_BATCH_EXIT=%s\n" "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== R: stageB 6 folds x 6 methods ==="
Rscript --vanilla R/experimental/crit_v7_stageB_mechanistic.R "$ORIGINAL" "$FIX3" \
 > "$OUT/stageB_R.log" 2>&1 || {
 tail -n 100 "$OUT/stageB_R.log" >&2; exit 1;
}
grep -qx "STAGEB_36_METHOD_FOLDS=PASS" "$OUT/stageB_R.log"
grep -qx "STAGEB_GLOBAL_GLM_EXACT_REPRESENTATION=PASS" "$OUT/stageB_R.log"
tail -n 20 "$OUT/stageB_R.log"

echo "=== Python: independent archival patient-level audit ==="
python3 scripts/crit_v7_stageB_independent_audit.py "$ORIGINAL" "$FIX3" \
 > "$OUT/stageB_independent_audit.log" 2>&1 || {
 cat "$OUT/stageB_independent_audit.log" >&2; exit 1;
}
grep -qx "STAGEB_36_VECTORS_INDEPENDENTLY_AUDITED=PASS" \
 "$OUT/stageB_independent_audit.log"
grep -qx "STAGEB_EXISTING_METHOD_PATIENT_PARITY=PASS" \
 "$OUT/stageB_independent_audit.log"
cat "$OUT/stageB_independent_audit.log"
echo "STAGEB_ALL_GATES=PASS"
echo "STAGEB_ORIGINAL_RLR_AND_FIX3_UNMODIFIED=TRUE"
echo "STAGEB_REPORT=$ROOT/$OUT"
