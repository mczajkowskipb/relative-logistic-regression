#!/usr/bin/env bash
# Controlled exhaustive32 pool coverage ablation. Prior RLR/v9 untouched.
set -Eeuo pipefail
trap 'rc=$?; echo "STAGEG_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 1 || {
 echo "STAGEG_STOP: pinned CRIT9_WORKTREE argument required" >&2; exit 2;
}
ROOT="$(git rev-parse --show-toplevel)"
V9="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
BRANCH="analysis/crit-v9-stageG-exhaustive32-ablation-v1"
ORIG_V9="7b32f804b7be92158a14fd5645eb65d56e7da763"
echo "STAGEG_SOURCE_SHA=$HEAD"
test "$(git rev-parse "refs/remotes/origin/$BRANCH")" = "$HEAD" || {
 echo "STAGEG_STOP: remote pinned branch mismatches checkout" >&2; exit 2;
}
test "$(git -C "$V9" rev-parse HEAD)" = "$ORIG_V9"
grep -qx "CRIT9_BATCH_EXIT=0" \
 "$V9/results/crit_v9_d_crossfit_2x3_v1/status.txt"
for tree in "$ROOT" "$V9"; do
 test -z "$(git -C "$tree" -c core.filemode=false \
           status --porcelain --untracked-files=no)" || {
   echo "STAGEG_STOP: tracked files modified $tree" >&2; exit 2;
 }
done
for path in \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/experimental/crit_v8_single_nomination.R \
 R/experimental/crit_v9_d_crossfit_one_shot.R \
 R/experimental/rlmt_v4_rank_repair.R \
 R/new/07_world_protocol_utils.R \
 R/new/28_rlmt_utils.R \
 R/new/40_rlmt_dc_utils.R \
 R/new/42_rlmt_dc_run.R
do
 test "$(git hash-object "$path")" = \
      "$(git rev-parse "$ORIG_V9:$path")" || {
  echo "STAGEG_STOP: protected source differs from pinned original: $path" >&2
  exit 2;
 }
done
for bin in Rscript python3; do command -v "$bin" >/dev/null; done
Rscript --vanilla -e '
 parse("R/experimental/critG_exhaustive32_one_shot.R")
 parse("R/experimental/critG_48_study_panel.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/critG_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("STAGEG_SOURCE_PARSE=PASS")
PY
OUT="results/crit_v9_stageG_exhaustive32_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "STAGEG_STOP: refusing overwrite of existing results" >&2; exit 2;
}
if [ "${STAGEG_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "STAGEG_ARCHIVE_AND_CODE_PREFLIGHT=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/SOURCE_SHA.txt"
printf '%s\n' "$ORIG_V9" > "$OUT/ORIGINAL_V9_SHA.txt"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' \
 "32-GENE-ONLY exhaustive combinatorial 496-pair D proposal pool; frozen v9 6-nominee screen, coefficient CV selection and five independent V gates; 24 reused and 24 new synthetic DEVELOPMENT studies; exact original-RLR parent; fixed seeds; 48x6x120 saved heldout patient predictions; source/model archives never changed; not independent method validation and not scalable to multi-omics" \
 > "$OUT/PROTOCOL.txt"
on_exit(){
 rc=$?
 printf 'STAGEG_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 1. 48-STUDY COMPLETE-POOL INTERVENTION ==="
Rscript --vanilla R/experimental/critG_48_study_panel.R "$V9" \
 > "$OUT/48_studies.log" 2>&1 || {
 tail -n 130 "$OUT/48_studies.log" >&2; exit 1;
}
cat "$OUT/48_studies.log"
grep -qx "STAGEG_48_STUDIES_288_METHODS=PASS" "$OUT/48_studies.log"
grep -qx "STAGEG_34560_TEST_PATIENT_VECTORS=PASS" "$OUT/48_studies.log"
grep -qx "STAGEG_EXHAUSTIVE_496_RELATIONS=PASS" "$OUT/48_studies.log"
grep -qx "STAGEG_ORIGINAL_RLR_REUSED24_PATIENT_PARITY=PASS" "$OUT/48_studies.log"

echo "=== 2. INDEPENDENT 34560-PATIENT AUDIT ==="
python3 scripts/critG_independent_audit.py "$V9" \
 > "$OUT/independent_audit.log" 2>&1 || {
 tail -n 150 "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "STAGEG_34560_PATIENT_METRICS_INDEPENDENTLY_AUDITED=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEG_48_STUDY_NOMINATION_AND_FALLBACK_AUDITED=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEG_24_REUSED_STUDIES_EXACT_ORIGINAL_RLR_PARITY=PASS" \
 "$OUT/independent_audit.log"

echo "STAGEG_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=SYNTHETIC_DEVELOPMENT_ONLY"
echo "ORIGINAL_RLR_CRIT9_AND_STAGEF_UNMODIFIED=TRUE"
echo "STAGEG_REPORT=$ROOT/$OUT"
