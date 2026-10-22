#!/usr/bin/env bash
# EXPERIMENTAL, NO RETUNING. 2 predeclared CompCancer datasets x first 3
# identical original outer folds x H1/H2 = 12 REAL fits; original 285
# RLR/RLMT models read-only. All outcome scores are development-exposed.
#
# Args: ORIGINAL_REF_REPO STAGE2_V2_WORKTREE
set -Eeuo pipefail
test "$#" -eq 2 || {
 echo "USAGE: $0 /ABS/ORIGINAL_REF_REPO /ABS/STAGE2_V2_WORKTREE" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
STAGE2="$(cd "$2" && pwd -P)"
cd "$ROOT"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
OLD="4873878f57f806bf5b52144756f8087b17fa540c"
HEAD="$(git rev-parse HEAD)"
test "$(git rev-parse origin/main)" = "$PAPER" || {
 echo "V2_PILOT_STOP: original paper SHA changed" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-journal-benchmark-pilot)" = "$HEAD" || {
 echo "V2_PILOT_STOP: mismatched benchmark code SHA" >&2; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$OLD" || {
 echo "V2_PILOT_STOP: archived original 285 model source SHA drifted" >&2
 exit 2
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "V2_PILOT_STOP: local tracked changes" >&2; exit 2;
}
STAGE2_STATUS="$STAGE2/results/rlmt_v2_stage2_metadata_v1/status.txt"
test -f "$STAGE2_STATUS" && \
 grep -qx "RLMT_V2_STAGE2_EXIT=0" "$STAGE2_STATUS" || {
 echo "V2_PILOT_STOP: v2 stage2 has not passed" >&2; exit 2;
}
test -f "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt" && \
 grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt" || {
 echo "V2_PILOT_STOP: original reference experiment not completed" >&2
 exit 2
}
OUT="results/rlmt_v2_exposed_pilot_2x3_v1"
test ! -e "$OUT" || {
 echo "V2_PILOT_STOP: existing pilot evidence; refusing overwrite" >&2
 exit 2
}
for tool in Rscript python3 sha256sum; do
 command -v "$tool" >/dev/null || {
  echo "V2_PILOT_STOP: missing tool $tool" >&2
  exit 2
 }
done
echo "=== A. PARSE EXPERIMENTS ==="
Rscript --vanilla -e '
 parse(file="R/experimental/rlmt_v2_exposed_pilot_2x3.R")
' >/dev/null
python3 -m py_compile scripts/rlmt_v2_exposed_pilot_2x3_audit.py
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$PAPER" > "$OUT/PAPER_SHA.txt"
printf '%s\n' "$OLD" > "$OUT/REFERENCE_285_FIT_SHA.txt"
printf '%s\n' "chen-2002,gordon-2002;rep00_fold00-02;H1/H2" \
 > "$OUT/PREREGISTERED_PILOT.txt"
trap 'rc=$?; printf "RLMT_V2_PILOT_EXIT=%s\n" "$rc" > "'"$OUT"'/status.txt"' EXIT

echo "=== B. ALL 12 REAL V2 PAPER-COHORT FITS ==="
Rscript --vanilla R/experimental/rlmt_v2_exposed_pilot_2x3.R "$REF" \
 > "$OUT/fits.log" 2>&1 || {
 tail -n 90 "$OUT/fits.log"; exit 1;
}
tail -n 25 "$OUT/fits.log"
grep -qx "RLMT_V2_PILOT_COMPLETE=PASS" "$OUT/fits.log"

echo "=== C. INDEPENDENT HELDOUT PREDICTION AUDIT ==="
python3 scripts/rlmt_v2_exposed_pilot_2x3_audit.py \
 --reference-repo "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT_V2_PILOT_INDEPENDENT_AUDIT=PASS" \
 "$OUT/independent_audit.log"

echo "RLMT_V2_EXPOSED_PAPER_PILOT_COMPLETE=PASS"
echo "REPORT=$ROOT/$OUT"
