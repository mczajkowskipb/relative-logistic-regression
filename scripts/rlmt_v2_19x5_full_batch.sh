#!/usr/bin/env bash
# 19 previously EXPOSED real biomedical cohorts x original 5 outer folds,
# four PREDECLARED algorithms (H1/H2/Gini/LA), 380 NEW fit attempts.
# Run through nohup to survive SSH/VPN; completely separate result dir.
# All 285 previously archived RLR/RLMT models read-only.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "USAGE: bash $0 /absolute/frozen-original-ref" >&2
 exit 2
}
REF="$(cd "$1" && pwd -P)"
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
SHA="$(git rev-parse HEAD)"
for cmd in python3 Rscript sha256sum; do
 command -v "$cmd" >/dev/null || {
  echo "RLMT_V2_19X5_STOP: missing $cmd" >&2
  exit 2
 }
done
echo "=== 1. CODE PARSE ==="
Rscript --vanilla -e '
 parse("R/experimental/rlmt_v2_19x5_four_method_runner.R")
 parse("R/experimental/rlmt_v2_gini_lookahead.R")
 parse("R/experimental/rlmt_v2_aligned_offset.R")
' >/dev/null
python3 -m py_compile \
 scripts/rlmt_v2_19x5_preflight.py \
 scripts/rlmt_v2_19x5_independent_audit.py
echo "RLMT_V2_19X5_PARSE=PASS"

echo "=== 2. NO-WRITE SOURCE AND FOLD PREFLIGHT ==="
python3 scripts/rlmt_v2_19x5_preflight.py \
 --original-ref "$REF" --new-sha "$SHA"

OUT="results/rlmt_v2_19x5_four_method_v1"
test ! -e "$OUT"
mkdir -p "$OUT"
printf '%s\n' "$SHA" > "$OUT/CODE_SHA.txt"
printf '%s\n' "4079a5ac400f1d2cfded523c7000f655470925ce" \
 > "$OUT/PAPER_SHA.txt"
printf '%s\n' "4873878f57f806bf5b52144756f8087b17fa540c" \
 > "$OUT/ORIGINAL_SOURCE_SHA.txt"
printf '%s\n' \
 "EXPOSED 19x5; H1,H2,TSPDT_GINI,RLMT_LA; 4096 route screen; 20 LA; unchanged method constants" \
 > "$OUT/FROZEN_PREDECLARED_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'RLMT_V2_19X5_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT
echo "=== 3. 380 REAL MODELS / NESTED DISCOVERY AND CONFIRMATION ==="
Rscript --vanilla R/experimental/rlmt_v2_19x5_four_method_runner.R \
 "$REF" > "$OUT/380_fits.log" 2>&1 || {
 tail -n 100 "$OUT/380_fits.log"
 echo "RLMT_V2_19X5_INCOMPLETE: preserved all completed model artifacts and failure manifest"
 exit 1
}
tail -n 45 "$OUT/380_fits.log"
grep -qx 'RLMT_V2_19X5_RUN=PASS' "$OUT/380_fits.log"

echo "=== 4. INDEPENDENT 380-MODEL PATIENT-LEVEL REPLAY ==="
python3 scripts/rlmt_v2_19x5_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx 'RLMT_V2_19X5_ALL_380_NEW_MODELS_AUDITED=PASS' \
 "$OUT/independent_audit.log"
echo "RLMT_V2_19X5_ALL_METHODS_COMPLETE=PASS"
echo "EXPOSED_RETROSPECTIVE_ONLY=TRUE"
echo "PUBLISHED_MAIN_AND_285_ORIGINAL_MODELS_UNTOUCHED=TRUE"
echo "REPORT=$ROOT/$OUT"
