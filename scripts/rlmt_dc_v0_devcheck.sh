#!/usr/bin/env bash
# RLMT-DC-v0 DEV GATE: only synthetic data. NEVER run a biomedical benchmark.
set -Eeuo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
branch="$(git branch --show-current)"
if [[ "$branch" != "feature/rlmt-dc-v0" ]]; then
  echo "STOP: required branch feature/rlmt-dc-v0 (actual: $branch)" >&2
  exit 2
fi
if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
  echo "STOP: tracked working tree is dirty" >&2
  exit 2
fi
main_sha="4079a5ac400f1d2cfded523c7000f655470925ce"
if [[ "$(git rev-parse origin/main)" != "$main_sha" ]]; then
  echo "STOP: publication main changed" >&2
  exit 2
fi
command -v Rscript >/dev/null || {
  echo "STOP: Rscript unavailable" >&2
  exit 2
}

out="results/rlmt_dc_v0_devcheck_v2"
if [[ -e "$out" ]]; then
  echo "STOP: refusing to overwrite prior DEV check" >&2
  exit 2
fi
mkdir -p "$out"
printf '%s\n' "$(git rev-parse HEAD)" > "$out/CODE_SHA.txt"
printf '%s\n' "$main_sha" > "$out/PUBLICATION_MAIN_SHA.txt"
on_exit() {
  local rc=$?
  printf 'DEV_CHECK_EXIT=%s\n' "$rc" > "$out/status.txt"
}
trap on_exit EXIT

echo "=== SOURCE PARSE ==="
Rscript --vanilla -e '
  parse(file="R/new/40_rlmt_dc_utils.R")
  parse(file="R/new/41_rlmt_dc_regression.R")
  parse(file="R/new/42_rlmt_dc_run.R")
  parse(file="R/new/43_rlmt_dc_synthetic_smoke.R")
' > "$out/parser.log" 2>&1
echo "RLMT_DC_R_PARSE=PASS"

echo "=== PRIMITIVE REGRESSION (SYNTHETIC ONLY) ==="
Rscript --vanilla R/new/41_rlmt_dc_regression.R \
  > "$out/primitive.log" 2>&1 || {
    cat "$out/primitive.log"
    exit 1
  }
cat "$out/primitive.log"
grep -qx 'RLMT_DC_PRIMITIVES=PASS' "$out/primitive.log"

echo "=== END-TO-END SYNTHETIC SMOKE ==="
Rscript --vanilla R/new/43_rlmt_dc_synthetic_smoke.R \
  > "$out/synthetic_pipeline.log" 2>&1 || {
    cat "$out/synthetic_pipeline.log"
    exit 1
  }
cat "$out/synthetic_pipeline.log"
grep -qx 'RLMT_DC_SYNTHETIC_PIPELINE=PASS' "$out/synthetic_pipeline.log"
grep -qx 'RLMT_DC_SYNTHETIC_ACCEPTED_BRANCH=PASS' "$out/synthetic_pipeline.log"
grep -qx 'RLMT_DC_SYNTHETIC_L1_PARITY=PASS' "$out/synthetic_pipeline.log"
grep -qx 'ACCEPTED_BRANCH_FIXTURE=FORCED_FUNCTIONAL_NOT_FOR_PAPER' "$out/synthetic_pipeline.log"

echo "RLMT_DC_DEV_GATE=PASS"
echo "REPORT=$repo/$out"
