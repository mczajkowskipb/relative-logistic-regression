#!/usr/bin/env bash
# Isolated RLMT-v5 guarded six-real-fold retrospective pilot.
# Strict original RLR patient-level parity, no overwrites.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "RLMT5_STOP: supply immutable original RLR checkout" >&2
 exit 2
}
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
cd "$ROOT"
SHA="$(git rev-parse HEAD)"
REFSHA="4873878f57f806bf5b52144756f8087b17fa540c"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
test "$(git rev-parse origin/main)" = "$MAIN"
test "$(git rev-parse origin/feature/rlmt-v5-honest-router)" = "$SHA"
test "$(git -C "$REF" rev-parse HEAD)" = "$REFSHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
test -z "$(git -C "$REF" -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx "RUNNER_EXIT=0" \
 "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
for bin in Rscript python3 sha256sum; do
 command -v "$bin" >/dev/null
done
OUT="results/rlmt_v5_honest_router_2x3_v1"
test ! -e "$OUT" || {
 echo "RLMT5_STOP: experiment output exists; do not overwrite" >&2
 exit 2
}
echo "=== 1. R PARSE PREFLIGHT AND AUDIT PY SYNTAX ==="
Rscript --vanilla -e '
parse("R/experimental/rlmt_v4_rank_repair.R")
parse("R/experimental/rlmt_v5_honest_router.R")
parse("R/experimental/test_rlmt_v5_honest_router_synthetic.R")
parse("R/experimental/rlmt_v5_honest_router_2x3_pilot.R")
' >/dev/null
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/rlmt_v5_2x3_independent_audit.py")
compile(p.read_text(),str(p),"exec")
print("RLMT5_PY_AUDIT_PARSE=PASS")
PY
echo "RLMT5_R_PARSE=PASS"
mkdir -p "$OUT"
printf '%s\n' "$SHA" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$REFSHA" > "$OUT/FROZEN_ORIGINAL_RLR_SHA.txt"
printf '%s\n' \
 "DEVELOPMENT_EXPOSED: 3 train-only inner stratified folds, seed1403; all frozen v4 parameter defaults; candidates GLOBAL_CAL,CONF_ROUTE,ACTIVITY_ROUTE; 4 nonzero blend weights 0.25/0.5/0.75/1; SE penalty 0.5; selection cost 0.00025; minimum conservative improvement 0.00025; strict no benefit -> original RLR. No tuning on outer test." \
 > "$OUT/FROZEN_PROTOCOL.txt"
on_exit() {
 rc=$?
 printf 'RLMT5_PILOT_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. REAL SYNTHETIC PARENT AND GATE TESTS ==="
Rscript --vanilla \
  R/experimental/test_rlmt_v5_honest_router_synthetic.R \
  > "$OUT/synthetic.log" 2>&1 || {
 tail -n 100 "$OUT/synthetic.log"; exit 1;
}
tail -n 32 "$OUT/synthetic.log"
grep -qx "RLMT5_SYNTH_REAL_GUARDED_CV=PASS" "$OUT/synthetic.log"

echo "=== 3. TWO EXPOSED DATASETS, SIX REAL HELDOUT FOLDS ==="
Rscript --vanilla \
 R/experimental/rlmt_v5_honest_router_2x3_pilot.R "$REF" \
 > "$OUT/six_real_folds.log" 2>&1 || {
 tail -n 100 "$OUT/six_real_folds.log"; exit 1;
}
tail -n 40 "$OUT/six_real_folds.log"
grep -qx "RLMT5_PILOT_ORIGINAL_RLR_PARITY=PASS" \
 "$OUT/six_real_folds.log"
grep -qx "RLMT5_PILOT_ALL_TRAINING=PASS" \
 "$OUT/six_real_folds.log"

echo "=== 4. INDEPENDENT ORIGINAL PATIENT-LEVEL REPLAY ==="
python3 scripts/rlmt_v5_2x3_independent_audit.py \
 --original-ref "$REF" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx "RLMT5_PILOT_INDEPENDENT_PATIENT_AUDIT=PASS" \
 "$OUT/independent_audit.log"
grep -qx "RLMT5_PILOT_ORIGINAL_RLR_PARENT_PARITY=PASS" \
 "$OUT/independent_audit.log"
echo "RLMT5_HONEST_ROUTER_ALL_GATES=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY"
echo "NO_MUTATIONS_TO_ORIGINAL_RLR=TRUE"
echo "REPORT=$ROOT/$OUT"
