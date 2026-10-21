#!/usr/bin/env bash
# One isolated, NO-OLD-OMICS development gate for H1/H2.
# Operator's NEW detached v2 worktree. Previously published results untouched.
set -Eeuo pipefail
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
HEAD="$(git rev-parse HEAD)"
test "$(git rev-parse origin/main)" = "$PAPER" || {
 echo "V2_STOP: original published main changed" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-aligned-probability-v2)" = "$HEAD" || {
 echo "V2_STOP: experimental branch differs from HEAD" >&2; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "V2_STOP: tracked worktree changed" >&2; exit 2;
}
command -v Rscript >/dev/null || {
 echo "V2_STOP: missing Rscript" >&2; exit 2;
}
command -v python3 >/dev/null || {
 echo "V2_STOP: missing Python3" >&2; exit 2;
}
OUT="results/rlmt_v2_aligned_offset_synthetic_gate_v1"
test ! -e "$OUT" || {
 echo "V2_STOP: previous gate exists, preserving evidence" >&2; exit 2;
}
mkdir -p "$OUT"
printf "%s\n" "$HEAD" > "$OUT/CODE_SHA.txt"
printf "%s\n" "$PAPER" > "$OUT/PAPER_MAIN_SHA.txt"
trap 'rc=$?; printf "RLMT_V2_SYNTH_GATE_EXIT=%s\n" "$rc" > "'"$OUT"'/status.txt"' EXIT
echo "=== 1. R PARSE / PYTHON STATIC VALIDATION ==="
Rscript --vanilla -e '
 parse(file="R/experimental/rlmt_v2_aligned_offset.R");
 parse(file="R/experimental/test_rlmt_v2_real_synthetic.R")
' > "$OUT/parse.log" 2>&1 || {
 cat "$OUT/parse.log"; exit 1;
}
python3 -m py_compile scripts/rlmt_v2_gse25066_source_gate.py \
 > "$OUT/python_parse.log" 2>&1 || {
 cat "$OUT/python_parse.log"; exit 1;
}
echo "RLMT_V2_R_PARSER=PASS"

echo "=== 2. H1 PROCEDURE + H2 FIXED-OFFSET REAL SYNTHETIC ==="
Rscript --vanilla R/experimental/test_rlmt_v2_real_synthetic.R \
 > "$OUT/real_synthetic.log" 2>&1 || {
 cat "$OUT/real_synthetic.log"; exit 1;
}
cat "$OUT/real_synthetic.log"
for marker in \
 RLMT_V2_REAL_SYNTHETIC_H1_H2_GATE=PASS \
 H1_FIT_PROCEDURE_ALIGNED=TRUE \
 H2_FIXED_LAMBDA_OFFSET_OOF=TRUE \
 TEST_ONLY_FORCED_CONFIRMATION_DISCLOSED=TRUE \
 NO_BIOMEDICAL_DATA=TRUE; do
 grep -qx "$marker" "$OUT/real_synthetic.log" || {
  echo "V2_STOP: missing marker $marker" >&2; exit 1;
 }
done
echo "RLMT_V2_H1_H2_SYNTHETIC_GATE=PASS"
echo "V2_SYNTHETIC_REPORT=$ROOT/$OUT"
