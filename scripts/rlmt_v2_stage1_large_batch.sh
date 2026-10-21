#!/usr/bin/env bash
# Stage 1: REAL H1/H2 synthetic gate, THEN independent GEO source acquisition.
# No old 19 omics files, no RLR published main mutation, no clinical model fit.
# Designed for user-run 'nohup bash ...' to survive SSH/VPN disconnect.
set -Eeuo pipefail
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
HEAD="$(git rev-parse HEAD)"
test "$(git rev-parse origin/main)" = "$MAIN" || {
 echo "RLMT_V2_STOP: original paper main changed" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-aligned-probability-v2)" = "$HEAD" || {
 echo "RLMT_V2_STOP: experimental branch SHA changed" >&2; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "RLMT_V2_STOP: uncommitted model source changes" >&2; exit 2;
}
OUT="results/rlmt_v2_stage1_batch_v1"
COHORT="data/rlmt_v2_gse25066_metadata_only_v1"
test ! -e "$OUT" && test ! -e "$COHORT" || {
 echo "RLMT_V2_STOP: stage1 or source evidence already exists" >&2; exit 2;
}
mkdir -p "$OUT" data
printf '%s\n' "$HEAD" > "$OUT/EXPERIMENTAL_CODE_SHA.txt"
printf '%s\n' "$MAIN" > "$OUT/ORIGINAL_PAPER_SHA.txt"
trap 'rc=$?; printf "RLMT_V2_STAGE1_EXIT=%s\n" "$rc" > "'"$OUT"'/status.txt"' EXIT

echo "=== A. SYNTHETIC ONLY: H1/H2 REAL SWAP/GLMNET ==="
bash scripts/rlmt_v2_h1_h2_synthetic_gate.sh \
 > "$OUT/h1_h2_synthetic.log" 2>&1 || {
  tail -n 100 "$OUT/h1_h2_synthetic.log"; exit 1;
}
cat "$OUT/h1_h2_synthetic.log"
grep -qx 'RLMT_V2_H1_H2_SYNTHETIC_GATE=PASS' \
 "$OUT/h1_h2_synthetic.log"

echo "=== B. NCBI GEO SOURCE-ONLY MATRICES & SHA256 ==="
python3 scripts/rlmt_v2_gse25066_source_gate.py \
 --output-dir "$ROOT/$COHORT" \
 > "$OUT/geo_source.log" 2>&1 || {
  tail -n 100 "$OUT/geo_source.log"
  echo "RLMT_V2_SOURCE_GATE=NOT_PASSED" >&2
  exit 1
}
cat "$OUT/geo_source.log"
grep -qx "GEO_NEW_COHORT_METADATA_GATE=PASS" \
 "$OUT/geo_source.log"
grep -qx "GEO_CLINICAL_LABEL_AND_GROUP_GATE=PENDING_MANUAL_AUDIT" \
 "$OUT/geo_source.log"

echo "RLMT_V2_STAGE1_H1_H2_SYNTHETIC_AND_SOURCE_METADATA=PASS"
echo "CLINICAL_OUTCOME_MODEL_FITS=0"
echo "OLD_PATIENT_19X5_RESULT_ROOTS_UNTOUCHED=TRUE"
echo "CLINICAL_AND_PATIENT_OVERLAP_GATE=PENDING"
echo "REPORT=$ROOT/$OUT"
