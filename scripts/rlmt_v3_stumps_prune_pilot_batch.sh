#!/usr/bin/env bash
# Small next scientific gate: REAL synthetic training + independent 18-fit
# biomedical pilot (2 previously exposed cohorts x 3 original folds x
# 3 new method families). Old 285 original RLR models always read only.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "V3_STOP: provide original frozen RLR reference worktree" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
SHA="$(git rev-parse HEAD)"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
OLD_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-stumps-pruning-v3)" = "$SHA"
test "$(git -C "$REF" rev-parse HEAD)" = "$OLD_SHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for cmd in Rscript python3 sha256sum; do command -v "$cmd" >/dev/null; done
OUT="results/rlmt_v3_stumps_2x3_pilot_v1"
test ! -e "$OUT" || {
 echo "V3_STOP: previous result root exists — no overwrite" >&2
 exit 2
}
echo "=== 1. NEW V3 CODE PARSE ==="
Rscript --vanilla -e '
parse("R/experimental/rlmt_v3_stumps_pruned.R")
parse("R/experimental/test_rlmt_v3_stumps_pruned_synthetic.R")
parse("R/experimental/rlmt_v3_stumps_2x3_real_pilot.R")
' >/dev/null
python3 -m py_compile scripts/rlmt_v3_stumps_2x3_independent_audit.py
echo "RLMT_V3_R_PARSE=PASS"
mkdir -p "$OUT"
printf '%s\n' "$SHA" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$OLD_SHA" > "$OUT/ORIGINAL_RLR_SHA.txt"
printf '%s\n' \
 "3 fixed modes; TSP sampled 4096; min leaf 8; Beta strength 16; ridge lambda 1; max tree depth 2; validation prune" \
 > "$OUT/PREDECLARED_V3_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'RLMT_V3_PILOT_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. REAL SYNTHETIC MODELS / RDS REPLAY ==="
Rscript --vanilla R/experimental/test_rlmt_v3_stumps_pruned_synthetic.R \
 > "$OUT/synthetic_real_fit.log" 2>&1 || {
 tail -n 100 "$OUT/synthetic_real_fit.log"; exit 1;
}
cat "$OUT/synthetic_real_fit.log"
grep -qx "RLMT_V3_STUMP_PRUNE_SYNTH_REAL_GATE=PASS" \
 "$OUT/synthetic_real_fit.log"

echo "=== 3. 18 NEW FITS ON ORIGINAL PAIRED BIOMEDICAL DATA ==="
Rscript --vanilla R/experimental/rlmt_v3_stumps_2x3_real_pilot.R \
 "$REF" > "$OUT/real_18_fits.log" 2>&1 || {
 tail -n 100 "$OUT/real_18_fits.log"; exit 1;
}
tail -n 40 "$OUT/real_18_fits.log"
grep -qx "RLMT_V3_PILOT_COMPLETE=PASS" "$OUT/real_18_fits.log"

echo "=== 4. INDEPENDENT REPLAY OF EXACT OUTER TEST PREDICTIONS ==="
python3 scripts/rlmt_v3_stumps_2x3_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT_V3_PILOT_INDEPENDENT_AUDIT=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT_V3_STUMP_AND_PRUNE_PILOT_COMPLETE=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "NO_MODIFICATIONS_TO_ORIGINAL_RLR=TRUE"
echo "REPORT=$ROOT/$OUT"
