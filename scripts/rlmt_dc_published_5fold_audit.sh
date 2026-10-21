#!/usr/bin/env bash
# READ-ONLY five-fold audit pipeline. No training; preserves previous outputs.
set -Eeuo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo"
out="results/rlmt_dc_published_data_5fold_v1"
test "$(git branch --show-current)" = "feature/rlmt-dc-v0" || {
  echo "STOP: wrong RLR branch" >&2; exit 2;
}
test -z "$(git status --porcelain --untracked-files=no)" || {
  echo "STOP: tracked worktree is dirty" >&2; exit 2;
}
test -f "$out/status.txt" && \
test "$(cat "$out/status.txt")" = "RUNNER_EXIT=0" || {
  echo "STOP: five-fold runner not completed successfully" >&2; exit 2;
}
test "$(cat "$out/PUBLICATION_MAIN_SHA.txt")" = \
  "4079a5ac400f1d2cfded523c7000f655470925ce" || {
  echo "STOP: publication main SHA changed" >&2; exit 2;
}
test "$(git rev-parse origin/main)" = \
  "4079a5ac400f1d2cfded523c7000f655470925ce" || {
  echo "STOP: origin/main moved" >&2; exit 2;
}
sha="$(cat "$out/CODE_SHA.txt")"
[[ "$sha" =~ ^[a-f0-9]{40}$ ]] || {
  echo "STOP: invalid frozen run SHA" >&2; exit 2;
}
for source in \
  R/new/07_world_protocol_utils.R \
  R/new/28_rlmt_utils.R \
  R/new/40_rlmt_dc_utils.R \
  R/new/42_rlmt_dc_run.R \
  R/new/47_rlmt_dc_published_5fold.R; do
  git show "$sha:$source" | cmp - "$source" >/dev/null || {
    echo "STOP: model source has changed since training: $source" >&2
    exit 2
  }
done
for report in \
  "$out/published_5fold_audit_v1" \
  "$out/published_5fold_accepted_replay_v1" \
  "$out/published_5fold_firstfold_parity_v1" \
  "$out/published_5fold_summary_v1"; do
  test ! -e "$report" || {
    echo "STOP: refusing to overwrite existing audit: $report" >&2
    exit 2
  }
done
Rscript --vanilla -e '
  parse(file="R/new/48_rlmt_dc_published_5fold_audit.R")
  parse(file="R/new/49_rlmt_dc_published_5fold_replay.R")
' >/dev/null
python3 -m py_compile \
  scripts/rlmt_dc_published_5fold_pilot_parity.py \
  scripts/rlmt_dc_published_5fold_summary.py

echo "=== 1. INDEPENDENT 285-FIT TEST RECONSTRUCTION ==="
Rscript --vanilla R/new/48_rlmt_dc_published_5fold_audit.R

echo "=== 2. ALL ACCEPTED TREES: SERIALIZED FIT REPLAY ==="
Rscript --vanilla R/new/49_rlmt_dc_published_5fold_replay.R

echo "=== 3. ALL 57 FIRST-FOLD PILOT PARITY CHECKS ==="
python3 scripts/rlmt_dc_published_5fold_pilot_parity.py

echo "=== 4. DATASET-LEVEL PAIRED FIVE-FOLD RESULTS ==="
python3 scripts/rlmt_dc_published_5fold_summary.py

echo "PUBLISHED_5FOLD_COMPLETE_AUDIT=PASS"
