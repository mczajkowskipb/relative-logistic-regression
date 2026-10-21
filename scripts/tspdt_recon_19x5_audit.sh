#!/usr/bin/env bash
# Read-only independent 95-model verification and paired RLR comparison.
set -Eeuo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
out="results/tspdt_recon_19x5_v1"
test -f "$out/status.txt" && \
test "$(cat "$out/status.txt")" = "RUNNER_EXIT=0" || {
  echo "STOP: TSPDT 95-fit experiment not finished successfully" >&2
  exit 2
}
test -z "$(git status --porcelain --untracked-files=no)" || {
  echo "STOP: tracked TSPDT worktree is dirty" >&2
  exit 2
}
test "$(git rev-parse origin/main)" = \
  "4079a5ac400f1d2cfded523c7000f655470925ce" || {
  echo "STOP: RLR publication main changed" >&2
  exit 2
}
sha="$(cat "$out/CODE_SHA.txt")"
test "$(git cat-file -t "$sha")" = "commit" || {
  echo "STOP: invalid code provenance" >&2
  exit 2
}
for path in \
  R/new/07_world_protocol_utils.R \
  R/baselines/tspdt_reconstruction_v1.R \
  R/baselines/50_tspdt_published_19x5_run.R; do
  git show "$sha:$path" | cmp - "$path" >/dev/null || {
    echo "STOP: trained model source differs from recorded original SHA ($path)" >&2
    exit 2
  }
done
test ! -e "$out/audit_v1" || {
  echo "STOP: audit_v1 already exists (no overwrite)" >&2
  exit 2
}
ref="$(cat "$out/REFERENCE_REPO.txt")"
test -f "$ref/results/rlmt_dc_published_data_5fold_v1/status.txt" || {
  echo "STOP: original experiment unavailable" >&2
  exit 2
}
Rscript --vanilla -e '
  parse(file="R/baselines/51_tspdt_published_19x5_audit.R")
' >/dev/null
echo "=== INDEPENDENT REAL TSPDT FIT REPLAY + RLR COMPARISON ==="
Rscript --vanilla R/baselines/51_tspdt_published_19x5_audit.R "$ref"
echo "TSPDT_PUBLISHED_COMPLETE_AUDIT=PASS"
