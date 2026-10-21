#!/usr/bin/env bash
# New experimental RLMT-CA-PS-v1: SYNTHETIC ONLY, isolated worktree.
# No real ARFF / previously published patient heldouts, no Git switching,
# no automatic R package installs and no overwrite of previous gate.
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
test "$(git rev-parse origin/main)" = \
 "4079a5ac400f1d2cfded523c7000f655470925ce" || {
 echo "STOP: RLR paper main different from frozen publication SHA" >&2
 exit 2
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "STOP: tracked worktree has content changes" >&2
 exit 2
}
command -v Rscript >/dev/null || {
 echo "STOP: Rscript is missing" >&2; exit 2;
}
OUT="results/rlmt_ca_ps_v1_synthetic_gate"
test ! -e "$OUT" || {
 echo "STOP: synthetic gate exists; refusing to overwrite" >&2
 exit 2
}
mkdir -p "$OUT"
git rev-parse HEAD > "$OUT/CODE_SHA.txt"
trap 'rc=$?; printf "RLMT_CA_PS_DEV_EXIT=%s\n" "$rc" > "'"$OUT"'/status.txt"' EXIT
echo "=== 1. SOURCE PARSING ==="
Rscript --vanilla -e '
  parse(file="R/experimental/rlmt_ca_ps_v1.R")
  parse(file="R/experimental/test_rlmt_ca_ps_v1.R")
  parse(file="R/experimental/test_rlmt_ca_ps_v1_real_synthetic.R")
' > "$OUT/parse.log" 2>&1 || {
 cat "$OUT/parse.log"
 exit 1
}
echo "RLMT_CA_PS_R_PARSE=PASS"
echo "=== 2. MATH/UNIT AND MOCK END-TO-END TEST ==="
Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1.R \
  > "$OUT/unit.log" 2>&1 || {
 cat "$OUT/unit.log"
 exit 1
}
cat "$OUT/unit.log"
grep -qx 'RLMT_CA_PS_SYNTHETIC_UNIT_AND_MOCK_INTEGRATION=PASS' "$OUT/unit.log"

echo "=== 3. REAL RLR/SWAP/GLMNET INTEGRATION — SYNTHETIC ONLY ==="
Rscript --vanilla R/experimental/test_rlmt_ca_ps_v1_real_synthetic.R \
  > "$OUT/real_synthetic.log" 2>&1 || {
 cat "$OUT/real_synthetic.log"
 exit 1
}
cat "$OUT/real_synthetic.log"
grep -qx 'RLMT_CA_PS_REAL_SYNTHETIC=PASS' "$OUT/real_synthetic.log"
grep -qx 'NO_BIOMEDICAL_DATA=TRUE' "$OUT/real_synthetic.log"
echo "RLMT_CA_PS_DEV_GATE=PASS"
echo "REPORT=$ROOT/$OUT"
