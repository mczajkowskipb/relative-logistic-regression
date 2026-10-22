#!/usr/bin/env bash
# Stage 2: NO biomedical fit. Safe with nohup; no older v1 evidence changes.
# Proves H1/H2 regression after glmnet API-only fix and profiles GEO
# discovery observed labels + validation IDs without evaluating validation truth.
set -Eeuo pipefail
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
test "$(git rev-parse origin/main)" = "$PAPER" || {
 echo "V2_STAGE2_STOP: paper main SHA divergence" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-aligned-probability-v2)" = "$HEAD" || {
 echo "V2_STAGE2_STOP: new research branch SHA divergence" >&2; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "V2_STAGE2_STOP: tracked source dirty" >&2; exit 2;
}
V1="results/rlmt_v2_stage1_batch_v1"
test -f "$V1/status.txt" && \
 grep -qx "RLMT_V2_STAGE1_EXIT=0" "$V1/status.txt" || {
 echo "V2_STAGE2_STOP: first source/synthetic gate not complete" >&2; exit 2;
}
SOURCE="data/rlmt_v2_gse25066_metadata_only_v1"
test -f "$SOURCE/PROVENANCE_METADATA_ONLY.json" || {
 echo "V2_STAGE2_STOP: original downloaded GEO SHA manifest missing" >&2; exit 2;
}
OUT="results/rlmt_v2_stage2_metadata_v1"
AUDIT="$OUT/clinical_schema_discovery_only"
test ! -e "$OUT" && \
 test ! -e "results/rlmt_v2_aligned_offset_synthetic_gate_v2" || {
 echo "V2_STAGE2_STOP: evidence root exists; refuse overwrite" >&2; exit 2;
}
for app in Rscript python3 sha256sum; do
 command -v "$app" >/dev/null || {
  echo "V2_STAGE2_STOP: required tool missing: $app" >&2; exit 2;
 }
done
mkdir "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$PAPER" > "$OUT/PAPER_MAIN_SHA.txt"
on_exit() {
 rc=$?
 printf "RLMT_V2_STAGE2_EXIT=%s\n" "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== A. SYNTHETIC H1/H2 V2 — after glmnet control fix ==="
bash scripts/rlmt_v2_h1_h2_synthetic_gate_v2.sh \
 > "$OUT/real_synthetic_v2.log" 2>&1 || {
  tail -n 90 "$OUT/real_synthetic_v2.log"; exit 1;
}
cat "$OUT/real_synthetic_v2.log"
grep -qx "RLMT_V2_H1_H2_SYNTHETIC_GATE=PASS" \
 "$OUT/real_synthetic_v2.log"

echo "=== B. CLINICAL DISCOVERY-ONLY FIELD SCHEMA + ID AUDIT ==="
python3 -m py_compile \
 scripts/rlmt_v2_gse25066_clinical_metadata_audit.py \
 > "$OUT/python_parse.log" 2>&1 || {
 cat "$OUT/python_parse.log"; exit 1;
}
python3 scripts/rlmt_v2_gse25066_clinical_metadata_audit.py \
 --source-dir "$ROOT/$SOURCE" \
 --output-dir "$ROOT/$AUDIT" \
 > "$OUT/clinical_metadata.log" 2>&1 || {
  cat "$OUT/clinical_metadata.log"; exit 1;
}
cat "$OUT/clinical_metadata.log"
grep -qx "GSE25066_CLINICAL_SCHEMA_AUDIT=PASS" \
 "$OUT/clinical_metadata.log"
grep -qx "GSE25066_VALIDATION_RESPONSE_VALUES_NOT_EVALUATED=TRUE" \
 "$OUT/clinical_metadata.log"

echo "RLMT_V2_STAGE2_SYNTHETIC_AND_DISCOVERY_METADATA=PASS"
echo "OLD_19x5_AND_CA_PS_ROOTS_UNTOUCHED=TRUE"
echo "VALIDATION_OUTCOME_EXAMINATION=NONE"
echo "BIOMEDICAL_FITS=0"
echo "CLINICAL_PATIENT_INDEPENDENCE=UNPROVEN"
echo "REPORT=$ROOT/$OUT"
