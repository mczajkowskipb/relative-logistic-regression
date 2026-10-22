#!/usr/bin/env bash
# Safe original-RLR exact-parity full 19x5 comparative benchmark.
# 95 training bundles × 5 outputs, 19 prior-development datasets.
# No writes to frozen source or earlier 2x3 pilot.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "RLMT4_19X5_STOP: one argument, old original RLR reference" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
REF_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-v4-19x5)" = "$HEAD"
test "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for cmd in Rscript python3 sha256sum; do
 command -v "$cmd" >/dev/null
done
OUT="results/rlmt_v4_rank_repair_19x5_v1"
test ! -e "$OUT" || {
 echo "RLMT4_19X5_STOP: previous result root exists, never overwrite" >&2
 exit 2
}
echo "=== 1. ORIGINAL RLR AND RLMT-v4 FROZEN-CODE PREFLIGHT ==="
Rscript --vanilla -e '
 parse("R/experimental/rlmt_v4_rank_repair.R")
 parse("R/experimental/test_rlmt_v4_rank_repair_synthetic.R")
 parse("R/experimental/rlmt_v4_rank_repair_19x5_runner.R")
' >/dev/null
# Python syntax check without writing __pycache__ in the source worktree.
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/rlmt_v4_19x5_independent_audit.py")
compile(p.read_text(),str(p),"exec")
print("RLMT4_19X5_PY_AUDITOR_SYNTAX=PASS")
PY
echo "RLMT4_19X5_R_PARSE=PASS"

mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REF_SHA" > "$OUT/IMMUTABLE_ORIGINAL_REF_SHA.txt"
printf '%s\n' \
 "RLMT-V4 same original original-RLR parent; exact parity 1e-9; unmodified seed42 4096 pair candidate pool; original 19x5; all 5 variants; beta penalty .02 and gamma penalty .10; nested parent 3fold OOF; previous exposure flagged; no outcome-directed tuning" \
 > "$OUT/FROZEN_PROTOCOL.txt"
on_exit() {
 code=$?
 printf 'RLMT4_19X5_BATCH_EXIT=%s\n' "$code" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. SYNTHETIC RLR-LITERAL PARITY ==="
Rscript --vanilla \
 R/experimental/test_rlmt_v4_rank_repair_synthetic.R \
 > "$OUT/synthetic_real.log" 2>&1 || {
 tail -n 90 "$OUT/synthetic_real.log"; exit 1;
}
tail -n 24 "$OUT/synthetic_real.log"
grep -qx "RLMT4_SYNTHETIC_REAL_RLR_PARITY_GATE=PASS" \
 "$OUT/synthetic_real.log"

echo "=== 3. MATCHED ALL 19 ORIGINAL DATASETS × FIVE FOLDS ==="
Rscript --vanilla \
 R/experimental/rlmt_v4_rank_repair_19x5_runner.R "$REF" \
 > "$OUT/475_method_folds.log" 2>&1 || {
 tail -n 110 "$OUT/475_method_folds.log"; exit 1;
}
tail -n 30 "$OUT/475_method_folds.log"
grep -qx "RLMT4_19X5_475_FOLD_METHODS=PASS" \
 "$OUT/475_method_folds.log"
grep -qx "RLMT4_19X5_ORIGINAL_RLR_PARITY=PASS" \
 "$OUT/475_method_folds.log"

echo "=== 4. INDEPENDENT PATIENT PARITY + 19 DATASET ROBUSTNESS ==="
python3 scripts/rlmt_v4_19x5_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT4_19X5_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "RLMT4_19X5_475_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT4_19X5_FULL_AUDIT=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "NO_ORIGINAL_RLR_RESULTS_WRITTEN=TRUE"
echo "REPORT=$ROOT/$OUT"
