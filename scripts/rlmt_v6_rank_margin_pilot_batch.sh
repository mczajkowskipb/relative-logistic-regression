#!/usr/bin/env bash
# Safe RLMT-v6 original-RLR-anchored rank margin 2x3 exposed pilot.
# Refuses existing result directory and any change to frozen original.
set -Eeuo pipefail
test "$#" -eq 1 || {
  echo "RLMT6_STOP: provide immutable ORIGINAL RLR checkout" >&2
  exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
REF_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-v6-rank-margin)" = "$HEAD"
test "$(git -C "$REF" rev-parse HEAD)" = "$REF_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for bin in Rscript python3 sha256sum; do
 command -v "$bin" >/dev/null
done

OUT="results/rlmt_v6_original_rank_margin_2x3_v1"
test ! -e "$OUT" || {
 echo "RLMT6_STOP: result already exists, do not overwrite" >&2
 exit 2
}

echo "=== 1. Parse all code and independent auditor ==="
Rscript --vanilla -e '
  parse("R/experimental/rlmt_v4_rank_repair.R")
  parse("R/experimental/rlmt_v6_original_rank_margin.R")
  parse("R/experimental/test_rlmt_v6_original_rank_margin_synthetic.R")
  parse("R/experimental/rlmt_v6_original_rank_margin_2x3_pilot.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/rlmt_v6_2x3_independent_audit.py")
compile(p.read_text(),str(p),"exec")
print("RLMT6_PY_AUDITOR_SYNTAX_PASS")
PY
echo "RLMT6_R_PARSE=PASS"

mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REF_SHA" > "$OUT/ORIGINAL_ARCHIVED_RLR_SHA.txt"
printf '%s\n' \
 "EXPOSED_RETROSPECTIVE_DEV_ONLY; locked 5 variants ORIGINAL_RLR GLOBAL_CAL ONEPAIR_MARGIN CONSENSUS_MARGIN ALLPAIR_MARGIN; original K<=10 genes only; 3-fold original RLR OOF; within ternary category normalized rank distances; train-only centering and scale; ridge beta=.02 gamma=.15; clip 3 feature 2 logodds; unchanged original root; only chen-2002 and gordon-2002 first 3 seed42 original outer folds" \
 > "$OUT/FROZEN_RESEARCH_PROTOCOL.txt"
on_exit() {
  rc=$?
  printf 'RLMT6_PILOT_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. Real synthetic parity + RDS + rank affine tests ==="
Rscript --vanilla \
 R/experimental/test_rlmt_v6_original_rank_margin_synthetic.R \
 > "$OUT/synthetic.log" 2>&1 || {
 tail -n 120 "$OUT/synthetic.log"; exit 1;
}
tail -n 35 "$OUT/synthetic.log"
grep -qx "RLMT6_ORIGINAL_RLR_REAL_SYNTHETIC_PARITY=PASS" \
 "$OUT/synthetic.log"

echo "=== 3. Original matched 2 datasets x 3 folds x 5 models ==="
Rscript --vanilla \
 R/experimental/rlmt_v6_original_rank_margin_2x3_pilot.R "$REF" \
 > "$OUT/30_method_folds.log" 2>&1 || {
 tail -n 120 "$OUT/30_method_folds.log"; exit 1;
}
tail -n 35 "$OUT/30_method_folds.log"
grep -qx "RLMT6_PILOT_ALL_30_METHOD_FOLDS=PASS" \
 "$OUT/30_method_folds.log"
grep -qx "RLMT6_PILOT_EXACT_ARCHIVED_ORIGINAL_RLR=PASS" \
 "$OUT/30_method_folds.log"

echo "=== 4. Independent original-patient paired audit ==="
python3 scripts/rlmt_v6_2x3_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT6_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "RLMT6_30_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT6_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "ORIGINAL_RLR_UNMODIFIED=TRUE"
echo "REPORT=$ROOT/$OUT"
