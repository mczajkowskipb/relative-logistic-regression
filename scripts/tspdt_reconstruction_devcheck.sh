#!/usr/bin/env bash
# Synthetic only; run inside a detached worktree. Does not run on RLR cohorts.
set -Eeuo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
test "$(git rev-parse origin/main)" = \
  "4079a5ac400f1d2cfded523c7000f655470925ce" || {
  echo "STOP: publication main changed" >&2; exit 2;
}
test -z "$(git status --porcelain --untracked-files=no)" || {
  echo "STOP: tracked worktree is dirty" >&2; exit 2;
}
command -v Rscript >/dev/null || {
  echo "STOP: Rscript not available" >&2; exit 2;
}
test -f docs/TSPDT_RECONSTRUCTION_V1_METHOD.md || {
  echo "STOP: research method specification missing" >&2; exit 2;
}
out="results/tspdt_reconstruction_v1_synthetic_gate"
test ! -e "$out" || {
  echo "STOP: refusing to overwrite past DEV results" >&2; exit 2;
}
mkdir -p "$out"
git rev-parse HEAD > "$out/CODE_SHA.txt"
on_exit() {
  rc=$?
  printf 'TSPDT_DEV_EXIT=%s\n' "$rc" > "$out/status.txt"
}
trap on_exit EXIT

echo "=== R SOURCE PARSE ==="
Rscript --vanilla -e '
  parse(file="R/baselines/tspdt_reconstruction_v1.R")
  parse(file="R/baselines/bigtsp_pair_cart_v1.R")
  parse(file="R/baselines/test_tspdt_reconstruction_v1.R")
' > "$out/parse.log" 2>&1 || {
  cat "$out/parse.log"
  exit 1
}
echo "TSPDT_RECON_R_PARSE=PASS"

echo "=== RECONSTRUCTION SYNTHETIC TEST ==="
Rscript --vanilla R/baselines/test_tspdt_reconstruction_v1.R \
  > "$out/test.log" 2>&1 || {
  cat "$out/test.log"
  exit 1
}
cat "$out/test.log"
grep -qx 'TSPDT_RECON_SYNTHETIC=PASS' "$out/test.log"
echo "TSPDT_RECON_DEV_GATE=PASS"
echo "REPORT=$repo/$out"
