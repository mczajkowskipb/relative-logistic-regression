#!/usr/bin/env bash
# H1 D-only synthetic screening methods — no RLR parent fit, no V/E tests.
set -Eeuo pipefail
trap 'rc=$?; echo "H1_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 1 || { echo "H1_STOP: pinned completed H0 worktree required" >&2; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
H0="$(cd "$1" && pwd -P)"
cd "$ROOT"
BRANCH="analysis/crit-stageH1-bounded-screening-pilot-v1"
H0_SHA="21986f8ece37ae5056ea1123360e12c74bcf60d3"
HEAD="$(git rev-parse HEAD)"
echo "STAGEH1_SOURCE_SHA=$HEAD"
test "$HEAD" = "$(git rev-parse "refs/remotes/origin/$BRANCH")" || {
 echo "H1_STOP: source SHA does not match fetched origin branch" >&2; exit 2;
}
test "$(git -C "$H0" rev-parse HEAD)" = "$H0_SHA"
grep -qx "H0_BATCH_EXIT=0" \
 "$H0/results/crit_stageH0_symmetry_v1/status.txt"
for dir in "$ROOT" "$H0"; do
 test -z "$(git -C "$dir" status --porcelain --untracked-files=no)" || {
  echo "H1_STOP: tracked source modified $dir" >&2; exit 2;
 }
done
for path in \
 R/experimental/critG_exhaustive32_one_shot.R \
 R/experimental/crit_v9_d_crossfit_one_shot.R \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/new/40_rlmt_dc_utils.R
do
 test "$(git hash-object "$path")" = \
      "$(git rev-parse "$H0_SHA:$path")" || {
  echo "H1_STOP: prior scientific code modified: $path" >&2; exit 2;
 }
done
command -v Rscript >/dev/null
command -v python3 >/dev/null
Rscript --vanilla -e '
 parse("R/experimental/crit_stageH1_pool_guard_fix1.R")
 parse("R/experimental/test_crit_stageH1_fix1_guards.R")
 parse("R/experimental/crit_stageH1_screening_pilot.R")
' >/dev/null
echo "=== H1 FIX1 POOL-BUDGET REGRESSION ==="
Rscript --vanilla R/experimental/test_crit_stageH1_fix1_guards.R
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_stageH1_independent_audit.py")
compile(p.read_text(encoding="utf-8"),str(p),"exec")
print("STAGEH1_SOURCE_PARSE=PASS")
PY
OUT="results/crit_stageH1_bounded_screening_fix1_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "H1_STOP: destination already exists; never overwrite" >&2; exit 2;
}
if [ "${STAGEH1_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "STAGEH1_CODE_AND_HISTORY_PREFLIGHT=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$H0_SHA" > "$OUT/H0_SOURCE_SHA.txt"
printf '%s\n' "STAGEH1_HISTORICAL_FAILED_OUTPUT_UNTOUCHED=TRUE" > "$OUT/FAILURE_ARCHIVE_POLICY.txt"
printf '%s\n' \
 "SCREENING-ONLY D=120; true Bernoulli null offset 0.5 (not original-RLR OOF); 32/64/128/256/512 genes; XOR and null; 8 replicates; fixed random vs per-study random vs rank-covariance-weighted induced-gene graph, each 320-pair pool; exhaustive all 496 pairs only p32; shortlist six undirected disjoint pair-product interactions; all D data only; no parent fitting, V gate, holdout E, or formal type-I calibration; development only" \
 > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'STAGEH1_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT
echo "=== 1. D-ONLY 256 STRATEGY/RANKING RUNS — H1 FIX1 ==="
Rscript --vanilla R/experimental/crit_stageH1_screening_pilot.R \
 > "$OUT/256_screening_runs.log" 2>&1 || {
 tail -n 100 "$OUT/256_screening_runs.log" >&2; exit 1;
}
tail -n 35 "$OUT/256_screening_runs.log"
grep -qx "STAGEH1_80_INDEPENDENT_D_ONLY_SYNTHETIC_STUDIES=PASS" \
 "$OUT/256_screening_runs.log"
grep -qx "STAGEH1_256_SCREENING_RUNS_COMPLETE=PASS" \
 "$OUT/256_screening_runs.log"
echo "=== 2. INDEPENDENT STRATEGY AUDIT ==="
python3 scripts/crit_stageH1_independent_audit.py \
 > "$OUT/independent_audit.log" 2>&1 || {
 tail -n 100 "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "STAGEH1_256_INDEPENDENT_PROVENANCE_CHECKS=PASS" \
 "$OUT/independent_audit.log"
grep -qx "STAGEH1_D_ONLY_PAIR_BUDGETS_AND_SYMMETRIC_JOINT_AUDITED=PASS" \
 "$OUT/independent_audit.log"
echo "STAGEH1_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=D_ONLY_SYNTHETIC_SCREENING_DEVELOPMENT_ONLY"
echo "HISTORICAL_RLR_CRIT9_STAGEG_STAGEH0_UNMODIFIED=TRUE"
echo "RESULT=$ROOT/$OUT"
