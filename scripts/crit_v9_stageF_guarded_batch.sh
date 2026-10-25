#!/usr/bin/env bash
# CRIT-v9 Stage F: safe pinned offline mechanism isolation.
# Source cohort/test archive only for baseline replay; no prior mutation.
set -Eeuo pipefail
trap 'rc=$?; echo "STAGEF_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 1 || {
 echo "STAGEF_STOP: PINNED_CRIT9_CHECKOUT required" >&2; exit 2;
}
ROOT="$(git rev-parse --show-toplevel)"
V9="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
BRANCH="analysis/crit-v9-stageF-xor-oracle-diagnostics-v1"
V9_SHA="7b32f804b7be92158a14fd5645eb65d56e7da763"
echo "STAGEF_PINNED_SOURCE=$HEAD"
test "$(git rev-parse "refs/remotes/origin/$BRANCH")" = "$HEAD" || {
 echo "STAGEF_STOP: changed branch/HEAD" >&2; exit 2;
}
test "$(git -C "$V9" rev-parse HEAD)" = "$V9_SHA"
grep -qx "CRIT9_BATCH_EXIT=0" \
 "$V9/results/crit_v9_d_crossfit_2x3_v1/status.txt"
for tree in "$ROOT" "$V9"; do
 test -z "$(git -C "$tree" -c core.filemode=false \
           status --porcelain --untracked-files=no)" || {
   echo "STAGEF_STOP: tracked source modifications $tree" >&2; exit 2;
 }
done
for path in \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/experimental/crit_v8_single_nomination.R \
 R/experimental/crit_v9_d_crossfit_one_shot.R \
 R/experimental/rlmt_v4_rank_repair.R \
 R/experimental/crit_v9_synthetic_24panel.R \
 R/new/07_world_protocol_utils.R \
 R/new/28_rlmt_utils.R \
 R/new/40_rlmt_dc_utils.R \
 R/new/42_rlmt_dc_run.R
do
 test "$(git hash-object "$path")" = \
      "$(git rev-parse "$V9_SHA:$path")" || {
  echo "STAGEF_STOP: protected source changed: $path" >&2
  exit 2;
 }
done
for bin in Rscript python3; do command -v "$bin" >/dev/null; done
Rscript --vanilla -e '
 parse("R/experimental/crit_v9_stageF_xor_oracle.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_v9_stageF_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("STAGEF_CODE_PARSE=PASS")
PY
OUT="results/crit_v9_stageF_xor_oracle_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "STAGEF_STOP: result dir already exists; no overwrite" >&2; exit 2;
}
if [ "${STAGEF_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "STAGEF_ARCHIVED_SOURCE_PREFLIGHT=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$V9_SHA" > "$OUT/V9_SHA.txt"
printf '%s\n' "REUSED_SYNTHETIC_DEV_MECHANISM_ONLY; exact v9 12 null/XOR synthetic train/test cases; no E test labels for selection; D selected candidate vs planted gene01>gene02 gate and gene03>gene04 relation; oracle conditional tree and global interaction fit D + one V test and E heldout; unconditional oracle E diagnostic; 8640 heldout patient vectors; no threshold tuning; synthetic panel already exposed" > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'STAGEF_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 1. RECONSTRUCT PLANTED XOR AND V9 NOMINATION ==="
Rscript --vanilla R/experimental/crit_v9_stageF_xor_oracle.R \
 "$V9" > "$OUT/oracle_reconstruction.log" 2>&1 || {
 tail -n 130 "$OUT/oracle_reconstruction.log" >&2; exit 1;
}
cat "$OUT/oracle_reconstruction.log"
grep -qx "STAGEF_12_PAIRED_SYNTHETIC_STUDIES=PASS" \
 "$OUT/oracle_reconstruction.log"
grep -qx "STAGEF_8640_PATIENT_PREDICTIONS=PASS" \
 "$OUT/oracle_reconstruction.log"

echo "=== 2. INDEPENDENT PATIENT AUDIT ==="
python3 scripts/crit_v9_stageF_independent_audit.py "$V9" \
 > "$OUT/independent_audit.log" 2>&1 || {
 tail -n 150 "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "STAGEF_8640_INDEPENDENT_PATIENT_REPLAY=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEF_FROZEN_V9_PATIENT_PARITY=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEF_12_STUDY_DISCOVERY_RANK_PROVENANCE=PASS" \
 "$OUT/independent_audit.log"

echo "STAGEF_ALL_GATES=PASS"
echo "STAGEF_EXPERIMENT_LABEL=REUSED_SYNTHETIC_DEV_MECHANISM_ONLY"
echo "STAGEF_PRIOR_RLR_AND_CRIT9_UNMODIFIED=TRUE"
echo "STAGEF_REPORT=$ROOT/$OUT"
