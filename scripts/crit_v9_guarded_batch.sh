#!/usr/bin/env bash
# CRIT-v9: 24 new synthetic studies + original-fold 2x3 six-method audit.
# Large background-safe one-shot scientific development batch; no history overwrite.
set -Eeuo pipefail
trap 'rc=$?; echo "CRIT9_PREP_OR_RUNTIME_FAILURE rc=$rc line=$LINENO" >&2' ERR
test "$#" -eq 2 || {
 echo "CRIT9_STOP: ORIGINAL_RLR and CRIT8_PINNED_REF needed" >&2; exit 2;
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
V8="$(cd "$2" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
BRANCH="feature/crit-v9-d-crossfit-one-shot-v1"
ORIG_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
BASE_SHA="435b87ae527db82575d4b5cf800725250b53935d"
echo "CRIT9_PINNED_SOURCE_SHA=$HEAD"
test "$(git rev-parse "refs/remotes/origin/$BRANCH")" = "$HEAD" || {
 echo "CRIT9_STOP: worktree not pinned to fetched origin branch" >&2; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$ORIG_SHA"
test "$(git -C "$V8" rev-parse HEAD)" = "$BASE_SHA"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx "CRIT8_PILOT_BATCH_EXIT=0" \
 "$V8/results/crit_v8_one_shot_2x3_v1/status.txt"
for tree in "$ROOT" "$REF" "$V8"; do
 test -z "$(git -C "$tree" -c core.filemode=false status --porcelain --untracked-files=no)" || {
  echo "CRIT9_STOP: tracked source changes in $tree" >&2; exit 2;
 }
done
for path in \
 R/experimental/crit_v7_original_rlr_two_leaf.R \
 R/experimental/crit_v8_single_nomination.R \
 R/experimental/rlmt_v4_rank_repair.R \
 R/new/07_world_protocol_utils.R \
 R/new/28_rlmt_utils.R \
 R/new/40_rlmt_dc_utils.R \
 R/new/42_rlmt_dc_run.R
do
 test "$(git hash-object "$path")" = \
      "$(git rev-parse "$BASE_SHA:$path")" || {
  echo "CRIT9_STOP: frozen original RLR/v7/v8 dependency changed" >&2
  exit 2;
 }
done
for cmd in Rscript python3 sha256sum; do command -v "$cmd" >/dev/null; done
Rscript --vanilla -e '
 parse("R/experimental/crit_v9_d_crossfit_one_shot.R")
 parse("R/experimental/test_crit_v9_d_crossfit.R")
 parse("R/experimental/crit_v9_synthetic_24panel.R")
 parse("R/experimental/crit_v9_2x3_pilot.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
for filename in (
 "scripts/crit_v9_independent_audit.py",
 "scripts/crit_v9_synthetic_independent_audit.py"
):
 compile(Path(filename).read_text(encoding="utf-8"),filename,"exec")
print("CRIT9_CODE_PARSE=PASS")
PY
OUT="results/crit_v9_d_crossfit_2x3_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "CRIT9_STOP: output path exists, refusing overwrite" >&2; exit 2;
}
if [ "${CRIT9_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "CRIT9_FULL_PREFLIGHT_ONLY=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$ORIG_SHA" > "$OUT/ORIGINAL_RLR_SHA.txt"
printf '%s\n' "$BASE_SHA" > "$OUT/V8_SOURCE_SHA.txt"
printf '%s\n' "EXPOSED_BIOMEDICAL_DEVELOPMENT_ONLY; 4x6 new synthetic independent train-test studies; six model arms; D-top6 proposal then D-only cross-fitted correction rank then one nominated V tree candidate; five method-specific V tests; unchanged original RLR; 2x3 paper fold development pilot; independent patient-level auditors; synthetic scenarios and seeds fixed before running" \
 > "$OUT/PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'CRIT9_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 1. SYNTHETIC STRUCTURAL AND REPRODUCIBILITY ==="
Rscript --vanilla R/experimental/test_crit_v9_d_crossfit.R \
 > "$OUT/synthetic_structural.log" 2>&1 || {
 tail -n 100 "$OUT/synthetic_structural.log" >&2; exit 1;
}
tail -n 25 "$OUT/synthetic_structural.log"
grep -qx "CRIT9_SYNTHETIC_ALL_GATES=PASS" "$OUT/synthetic_structural.log"

echo "=== 2. 24 FRESH INDEPENDENT SYNTHETIC TRAIN-TEST STUDIES ==="
Rscript --vanilla R/experimental/crit_v9_synthetic_24panel.R \
 > "$OUT/synthetic_panel.log" 2>&1 || {
 tail -n 120 "$OUT/synthetic_panel.log" >&2; exit 1;
}
tail -n 30 "$OUT/synthetic_panel.log"
grep -qx "CRIT9_24_FRESH_SYNTHETIC_DATASETS_SIX_ARMS=PASS" \
 "$OUT/synthetic_panel.log"

echo "=== 3. SYNTHETIC PATIENT-LEVEL INDEPENDENT AUDIT ==="
python3 scripts/crit_v9_synthetic_independent_audit.py \
 > "$OUT/synthetic_independent_audit.log" 2>&1 || {
 tail -n 140 "$OUT/synthetic_independent_audit.log" >&2; exit 1;
}
cat "$OUT/synthetic_independent_audit.log"
grep -qx "CRIT9_17280_SYNTHETIC_PATIENT_PREDICTIONS_AUDITED=PASS" \
 "$OUT/synthetic_independent_audit.log"

echo "=== 4. TWO EXPOSED BIOMEDICAL COHORTS 3 ORIGINAL FOLDS x 6 ==="
Rscript --vanilla R/experimental/crit_v9_2x3_pilot.R "$REF" \
 > "$OUT/36_method_folds.log" 2>&1 || {
 tail -n 120 "$OUT/36_method_folds.log" >&2; exit 1;
}
tail -n 25 "$OUT/36_method_folds.log"
grep -qx "CRIT9_36_METHOD_FOLDS=PASS" "$OUT/36_method_folds.log"
grep -qx "CRIT9_EXACT_ORIGINAL_RLR_PARITY=PASS" "$OUT/36_method_folds.log"

echo "=== 5. ORIGINAL BIOMEDICAL PATIENT AUDIT ==="
python3 scripts/crit_v9_independent_audit.py "$REF" \
 > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log" >&2; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "CRIT9_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "CRIT9_D_CROSSFIT_36_VECTORS_AUDITED=PASS" \
 "$OUT/independent_audit.log"

echo "CRIT9_ALL_GATES=PASS"
echo "ORIGINAL_RLR_V7_V8_UNMODIFIED=TRUE"
echo "SCIENTIFIC_STAGE=DEVELOPMENT_ONLY"
echo "REPORT=$ROOT/$OUT"
