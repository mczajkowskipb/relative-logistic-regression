#!/usr/bin/env bash
# Safe CRIT-v7 true conditional gene-order tree plus global interaction 2x3 exposed pilot.
# Refuses existing result directory and any change to frozen original.
set -Eeuo pipefail
test "$#" -eq 1 || {
  echo "CRIT7_STOP: provide immutable ORIGINAL RLR checkout" >&2
  exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
REF_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/crit-v7-regime-switch-probe)" = "$HEAD"
test "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for bin in Rscript python3 sha256sum; do
 command -v "$bin" >/dev/null
done

OUT="results/crit_v7_original_rlr_2x3_v1"
test ! -e "$OUT" || {
 echo "CRIT7_STOP: result already exists, do not overwrite" >&2
 exit 2
}

echo "=== 1. Parse all code and independent auditor ==="
Rscript --vanilla -e '
  parse("R/experimental/rlmt_v4_rank_repair.R")
  parse("R/experimental/crit_v7_original_rlr_two_leaf.R")
  parse("R/experimental/test_crit_v7_original_rlr_synthetic.R")
  parse("R/experimental/crit_v7_original_rlr_2x3_pilot.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_v7_2x3_independent_audit.py")
compile(p.read_text(),str(p),"exec")
print("CRIT7_PY_AUDITOR_SYNTAX_PASS")
PY
echo "CRIT7_R_PARSE=PASS"

mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REF_SHA" > "$OUT/ORIGINAL_ARCHIVED_RLR_SHA.txt"
printf '%s\n' \
 "EXPOSED_RETROSPECTIVE_DEV_ONLY; locked 4 variants ORIGINAL_RLR GLOBAL_CAL GLOBAL_INTERACTION CRIT_TREE; exact original RLR root; train-only inner D/V split confirmation; unsupervised balance candidates and random proposals; genuine distinct leaf one-pair slope and leaf intercept; global non-tree pair interaction control; original 3fold inner OOF; fixed ridge; fail-closed fallback; unchanged original root; only chen-2002 and gordon-2002 first 3 seed42 original outer folds" \
 > "$OUT/FROZEN_RESEARCH_PROTOCOL.txt"
on_exit() {
  rc=$?
  printf 'CRIT7_PILOT_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. Real synthetic parity + RDS + rank affine tests ==="
Rscript --vanilla \
 R/experimental/test_crit_v7_original_rlr_synthetic.R \
 > "$OUT/synthetic.log" 2>&1 || {
 tail -n 120 "$OUT/synthetic.log"; exit 1;
}
tail -n 35 "$OUT/synthetic.log"
grep -qx "CRIT7_REAL_SYNTHETIC_ORIGINAL_ROOT_AND_TREE=PASS" \
 "$OUT/synthetic.log"

echo "=== 3. Original matched 2 datasets x 3 folds x 4 models ==="
Rscript --vanilla \
 R/experimental/crit_v7_original_rlr_2x3_pilot.R "$REF" \
 > "$OUT/24_method_folds.log" 2>&1 || {
 tail -n 120 "$OUT/24_method_folds.log"; exit 1;
}
tail -n 35 "$OUT/24_method_folds.log"
grep -qx "CRIT7_PILOT_ALL_24_METHOD_FOLDS=PASS" \
 "$OUT/24_method_folds.log"
grep -qx "CRIT7_PILOT_EXACT_ARCHIVED_ORIGINAL_RLR=PASS" \
 "$OUT/24_method_folds.log"

echo "=== 4. Independent original-patient paired audit ==="
python3 scripts/crit_v7_2x3_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "CRIT7_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "CRIT7_24_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS" \
 "$OUT/independent_audit.log"
echo "CRIT7_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "ORIGINAL_RLR_UNMODIFIED=TRUE"
echo "REPORT=$ROOT/$OUT"
