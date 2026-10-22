#!/usr/bin/env bash
# Synthetics ONLY. REAL glmnet/RLR leaf fits for capped lookahead and Gini.
# No old RLR/RLMT models, no biomedical data. Never overwrite results.
set -Eeuo pipefail
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
test "$(git rev-parse origin/main)" = "$PAPER" || {
 echo "RLMT_LA_STOP: frozen RLR paper main has changed" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-lookahead-v2)" = "$HEAD" || {
 echo "RLMT_LA_STOP: experimental commit does not match origin" >&2; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "RLMT_LA_STOP: tracked source is dirty" >&2; exit 2;
}
command -v Rscript >/dev/null || {
 echo "RLMT_LA_STOP: Rscript unavailable" >&2; exit 2;
}
OUT="results/rlmt_v2_gini_lookahead_synthetic_gate_v1"
test ! -e "$OUT" || {
 echo "RLMT_LA_STOP: existing result evidence; refusing overwrite" >&2
 exit 2
}
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$PAPER" > "$OUT/PAPER_SHA.txt"
on_exit() {
 rc=$?
 printf 'RLMT_LA_SYNTHETIC_GATE_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT
echo "=== A. STRICT R PARSE ==="
Rscript --vanilla -e '
 parse("R/experimental/rlmt_v2_gini_lookahead.R")
 parse("R/experimental/test_rlmt_v2_gini_lookahead_synthetic.R")
' > "$OUT/parse.log" 2>&1 || {
 cat "$OUT/parse.log"; exit 1;
}
echo "RLMT_LA_R_PARSE=PASS"
echo "=== B. REAL GINI AND INNER-OOF RLR LOOKAHEAD (SYNTHETIC) ==="
Rscript --vanilla R/experimental/test_rlmt_v2_gini_lookahead_synthetic.R \
 > "$OUT/real_synthetic.log" 2>&1 || {
 tail -n 120 "$OUT/real_synthetic.log"; exit 1;
}
tail -n 35 "$OUT/real_synthetic.log"
grep -qx "RLMT_V2_GINI_LA_SYNTHETIC_REAL_GATE=PASS" \
 "$OUT/real_synthetic.log"
grep -qx "FORCED_NOMINATION_AND_CONFIRMATION_TEST_ONLY=TRUE" \
 "$OUT/real_synthetic.log"
grep -qx "NO_BIOMEDICAL_DATA=TRUE" "$OUT/real_synthetic.log"
echo "RLMT_LA_SYNTHETIC_BATCH=PASS"
echo "REPORT=$ROOT/$OUT"
