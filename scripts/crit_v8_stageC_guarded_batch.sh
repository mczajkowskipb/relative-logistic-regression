#!/usr/bin/env bash
# Run Stage C as a strict archived-RDS diagnosis. NO training.
set -Eeuo pipefail
trap 'rc=$?; echo "STAGEC_PREP_OR_RUNTIME_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 3 || {
  echo "STAGEC_STOP: ORIGINAL_RLR CRIT7_FIX3 CRIT8 worktrees required" >&2
  exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
V7="$(cd "$2" && pwd -P)"
V8="$(cd "$3" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
BRANCH="analysis/crit-v8-stageC-selection-audit-v1"
V7_SHA="e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03"
V8_SHA="435b87ae527db82575d4b5cf800725250b53935d"
ORIG_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
echo "STAGEC_SOURCE_SHA=$HEAD"
test "$(git rev-parse "refs/remotes/origin/$BRANCH")" = "$HEAD" || {
 echo "STAGEC_STOP: remote branch moved; cannot run unpinned code" >&2; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$ORIG_SHA"
test "$(git -C "$V7" rev-parse HEAD)" = "$V7_SHA"
test "$(git -C "$V8" rev-parse HEAD)" = "$V8_SHA"
for worktree in "$ROOT" "$REF" "$V7" "$V8"; do
  test -z "$(git -C "$worktree" -c core.filemode=false status --porcelain --untracked-files=no)" || {
   echo "STAGEC_STOP: tracked modifications in $worktree" >&2; exit 2;
  }
done
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx "CRIT7_PILOT_BATCH_EXIT=0" \
 "$V7/results/crit_v7_original_rlr_2x3_v1/status.txt"
grep -qx "CRIT8_PILOT_BATCH_EXIT=0" \
 "$V8/results/crit_v8_one_shot_2x3_v1/status.txt"
for script in R/experimental/crit_v8_single_nomination.R \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/experimental/rlmt_v4_rank_repair.R; do
  test "$(git hash-object "$script")" = "$(git rev-parse "$V8_SHA:$script")"
done
for ds in chen-2002 gordon-2002; do
 for num in 00 01 02; do
  test -f "$V7/results/crit_v7_original_rlr_2x3_v1/models/$ds/rep00_fold$num/complete_crit7.rds"
  test -f "$V8/results/crit_v8_one_shot_2x3_v1/models/$ds/rep00_fold$num/complete_crit8.rds"
 done
done
command -v Rscript >/dev/null
command -v python3 >/dev/null
Rscript --vanilla -e 'parse("R/experimental/crit_v8_stageC_archived_selection.R")' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_v8_stageC_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("STAGEC_SYNTAX_PARSE=PASS")
PY
echo "STAGEC_SOURCE_AND_ARCHIVE_PREFLIGHT=PASS"

OUT="results/crit_v8_stageC_selection_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "STAGEC_STOP: result root exists; no overwrite"; exit 2;
}
if [ "${STAGEC_PREFLIGHT_ONLY:-0}" = "1" ]; then
 echo "STAGEC_PREFLIGHT_ONLY=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$ORIG_SHA" > "$OUT/ORIGINAL_SHA.txt"
printf '%s\n' "$V7_SHA" > "$OUT/V7_SHA.txt"
printf '%s\n' "$V8_SHA" > "$OUT/V8_SHA.txt"
printf '%s\n' "READ_ONLY_FROM_SAVED_MODELS; EXPOSED_RETROSPECTIVE_DEV_ONLY; no training; v7 D top6 candidate scores vs v8 D rank1; five V arms per fold, 30 total decisions; all six original folds; archived outer-test metrics only described" > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'STAGEC_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT
echo "=== 1. SAVED RDS SCIENTIFIC SELECTION DIAGNOSTIC ==="
Rscript --vanilla R/experimental/crit_v8_stageC_archived_selection.R \
 "$REF" "$V7" "$V8" > "$OUT/stageC_R.log" 2>&1 || {
  tail -n 120 "$OUT/stageC_R.log" >&2; exit 1;
}
cat "$OUT/stageC_R.log"
grep -qx "STAGEC_RDS_SELECTION_AUDIT=PASS" "$OUT/stageC_R.log"
grep -qx "STAGEC_V7_V8_D_FIRST_CANDIDATE_EXACT_PARITY=PASS" "$OUT/stageC_R.log"

echo "=== 2. INDEPENDENT 30 VALIDATION DECISION AUDIT ==="
python3 scripts/crit_v8_stageC_independent_audit.py "$V7" "$V8" \
 > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "STAGEC_ALL_30_VALIDATION_DECISIONS_AUDITED=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEC_FROZEN_MANIFEST_PARITY=PASS" \
 "$OUT/independent_audit.log"
echo "STAGEC_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "ARCHIVED_RLR_V7_V8_UNMODIFIED=TRUE"
echo "STAGEC_REPORT=$ROOT/$OUT"
