#!/usr/bin/env bash
# Safe isolated gate of RLMT-v4, NEVER writes old original RLR.
# Only after real synthetic fitted models passes does it touch six
# previous-paper, already EXPOSED outer training folds.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "RLMT4_STOP: provide original frozen RLR worktree" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
SHA="$(git rev-parse HEAD)"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
REFSHA="4873878f57f806bf5b52144756f8087b17fa540c"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-v4-rank-repair)" = "$SHA"
test "$(git -C "$REF" rev-parse HEAD)" = "$REFSHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for c in Rscript python3 sha256sum; do command -v "$c" >/dev/null; done
OUT="results/rlmt_v4_rank_repair_2x3_v1"
test ! -e "$OUT" || {
 echo "RLMT4_STOP: prior output already exists. Never overwrite." >&2
 exit 2
}
echo "=== 1. PRE-FLIGHT REAL R PARSE + AUDITOR SYNTAX ==="
Rscript --vanilla -e '
  parse("R/experimental/rlmt_v4_rank_repair.R")
  parse("R/experimental/test_rlmt_v4_rank_repair_synthetic.R")
  parse("R/experimental/rlmt_v4_rank_repair_2x3_pilot.R")
' >/dev/null
python3 -m py_compile scripts/rlmt_v4_2x3_independent_audit.py
echo "RLMT4_PARSE=PASS"
mkdir -p "$OUT"
printf '%s\n' "$SHA" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REFSHA" > "$OUT/OLD_ORIGINAL_RLR_SHA.txt"
printf '%s\n' \
 "EXPOSED development pilot, exact original RLR parent, crossfit OOF 3, 4096 random rank pairs + original K10, hard min side 8, ridge beta 0.02 gamma 0.1, no outer-test tuning" \
 > "$OUT/FROZEN_EXPLORATORY_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'RLMT4_PILOT_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. SYNTHETIC ORIGINAL RLR LITERAL PARITY GATE ==="
Rscript --vanilla R/experimental/test_rlmt_v4_rank_repair_synthetic.R \
 > "$OUT/synthetic_real.log" 2>&1 || {
 tail -n 100 "$OUT/synthetic_real.log"
 exit 1
}
tail -n 35 "$OUT/synthetic_real.log"
grep -qx "RLMT4_SYNTHETIC_REAL_RLR_PARITY_GATE=PASS" \
 "$OUT/synthetic_real.log"

echo "=== 3. REAL 2 DATASETS x 3 FOLDS x 5 VARIANTS ==="
Rscript --vanilla R/experimental/rlmt_v4_rank_repair_2x3_pilot.R "$REF" \
 > "$OUT/real_30_fits.log" 2>&1 || {
 tail -n 120 "$OUT/real_30_fits.log"
 exit 1
}
tail -n 40 "$OUT/real_30_fits.log"
grep -qx "RLMT4_PILOT_30_FOLD_METHODS=PASS" "$OUT/real_30_fits.log"
grep -qx "RLMT4_ORIGINAL_RLR_NUMERICAL_PARITY=PASS" "$OUT/real_30_fits.log"

echo "=== 4. INDEPENDENT PATIENT-LEVEL RESULTS AUDIT ==="
python3 scripts/rlmt_v4_2x3_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"
 exit 1
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT4_30_REAL_PILOT_PREDICTIONS_AUDITED=PASS" \
 "$OUT/independent_audit.log"
grep -qx "RLMT4_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT4_RANK_REPAIR_PILOT_ALL_CHECKS_PASS"
echo "EXPOSED_RETROSPECTIVE_NO_EXTERNAL_VALIDATION=TRUE"
echo "ORIGINAL_RLR_RESULTS_UNTOUCHED=TRUE"
echo "REPORT=$ROOT/$OUT"
