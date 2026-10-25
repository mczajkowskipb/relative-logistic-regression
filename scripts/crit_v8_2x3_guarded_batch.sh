#!/usr/bin/env bash
# CRIT-v8: frozen original RLR and 6 x 6 exposed-dev pilot.
# New worktree/output only. Never overwrite historical experiments.
set -Eeuo pipefail
trap 'rc=$?; echo "CRIT8_PREP_OR_RUNTIME_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 1 || { echo "CRIT8_STOP: provide ORIGINAL_RLR_CHECKOUT"; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
REF_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
BASE_SHA="e2b8116017f889bfdc1b633ccc3ac1beb9b5cc03"
MAIN_SHA="4079a5ac400f1d2cfded523c7000f655470925ce"
BRANCH="feature/crit-v8-single-nomination-controls-v1"
echo "CRIT8_PRECHECK_HEAD=$HEAD"
test "$(git rev-parse "origin/$BRANCH")" = "$HEAD" || {
 echo "CRIT8_STOP: source not pinned to remote branch SHA"; exit 2;
}
test "$(git rev-parse origin/main)" = "$MAIN_SHA" || {
 echo "CRIT8_STOP: publication main SHA differs"; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA" || {
 echo "CRIT8_STOP: original paper RLR SHA changed"; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"

# Byte-identical v7 scientific core and old RLR dependencies.
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
 test "$(git hash-object "$path")" = \
      "$(git rev-parse "$BASE_SHA:$path")" || {
  echo "CRIT8_STOP: previously frozen algorithm dependency changed: $path"
  exit 2
 }
done
echo "CRIT8_FROZEN_ORIGINAL_ALGORITHM=PASS"

for bin in Rscript python3 sha256sum; do
 command -v "$bin" >/dev/null
done
Rscript --vanilla -e '
 parse("R/experimental/crit_v8_single_nomination.R")
 parse("R/experimental/test_crit_v8_single_nomination.R")
 parse("R/experimental/crit_v8_2x3_pilot.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_v8_2x3_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("CRIT8_R_AND_PYTHON_PARSE=PASS")
PY

OUT="results/crit_v8_one_shot_2x3_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "CRIT8_STOP: output already exists; will not overwrite"; exit 2;
}
if [ "${CRIT8_PREFLIGHT_ONLY:-0}" = "1" ]; then
 echo "CRIT8_PREFLIGHT_ONLY=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REF_SHA" > "$OUT/ARCHIVED_RLR_SHA.txt"
printf '%s\n' "EXPOSED_RETROSPECTIVE_DEV_ONLY; six methods; one D-nominated gate and relational pairs; five separate one-shot V confirmation tests, no V-best-method selection; fixed source fold seed42; exact paper RLR root; fixed ridge and 1SE+complexity fallback; new outputs and patient-level audit; not statistically independent confirmation" \
 > "$OUT/FROZEN_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf "CRIT8_PILOT_BATCH_EXIT=%s\n" "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 1. REAL ORIGINAL RLR SYNTHETIC AND ONE-SHOT GATES ==="
Rscript --vanilla R/experimental/test_crit_v8_single_nomination.R \
 > "$OUT/synthetic.log" 2>&1 || {
 tail -n 110 "$OUT/synthetic.log" >&2; exit 1;
}
tail -n 25 "$OUT/synthetic.log"
grep -qx "CRIT8_SYNTHETIC_ALL_GATES=PASS" "$OUT/synthetic.log"

echo "=== 2. SIX EXPOSED FOLDS BY SIX METHODS ==="
Rscript --vanilla R/experimental/crit_v8_2x3_pilot.R "$REF" \
 > "$OUT/36_method_folds.log" 2>&1 || {
 tail -n 100 "$OUT/36_method_folds.log" >&2; exit 1;
}
tail -n 20 "$OUT/36_method_folds.log"
grep -qx "CRIT8_36_METHOD_FOLDS=PASS" "$OUT/36_method_folds.log"
grep -qx "CRIT8_EXACT_ORIGINAL_RLR_PARITY=PASS" "$OUT/36_method_folds.log"

echo "=== 3. INDEPENDENT ORIGINAL PATIENT AUDIT ==="
python3 scripts/crit_v8_2x3_independent_audit.py "$REF" \
 > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "CRIT8_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" "$OUT/independent_audit.log"
grep -qx "CRIT8_ONE_SHOT_36_VECTORS_AUDITED=PASS" "$OUT/independent_audit.log"
echo "CRIT8_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "ORIGINAL_RLR_UNMODIFIED=TRUE"
echo "REPORT=$ROOT/$OUT"
