#!/usr/bin/env bash
# Complete RETROSPECTIVE 19x5 NEW-method-only real omics experiment.
# Existing 285 original fit outputs + original biomedical inputs READ ONLY.
# One-command preflight, parse, all 95 fit, independent all-model audit.
# This is EXPOSED METHOD DEVELOPMENT: never claim prospective validation.
set -Eeuo pipefail
test "$#" -eq 1 || {
 echo "USAGE: bash scripts/rlmt_ca_ps_v1_retro_batch.sh /ABS/OLD_RLR_REPO" >&2
 exit 2
}
REF="$(cd "$1" && pwd -P)"
ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
NEW="results/rlmt_ca_ps_v1_retrospective_19x5"
MAIN="4079a5ac400f1d2cfded523c7000f655470925ce"
OLD="4873878f57f806bf5b52144756f8087b17fa540c"
HEAD="$(git rev-parse HEAD)"
test "$(git rev-parse origin/main)" = "$MAIN" || {
 echo "STOP: original publication main SHA changed" >&2; exit 2;
}
test "$(git rev-parse origin/feature/rlmt-class-aware-ps-v1)" = "$HEAD" || {
 echo "STOP: experimental branch remote HEAD mismatch" >&2; exit 2;
}
test "$(git -C "$REF" rev-parse HEAD)" = "$OLD" || {
 echo "STOP: original 285-model result source HEAD mismatch" >&2; exit 2;
}
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)" || {
 echo "STOP: new method tracked tree dirty" >&2; exit 2;
}
test ! -e "$NEW" || {
 echo "STOP: retrospective result root exists; refusing overwrite" >&2
 exit 2
}
for command_name in python3 Rscript sha256sum; do
 command -v "$command_name" >/dev/null || {
  echo "STOP: required command missing: $command_name" >&2
  exit 2
 }
done
test -f "results/rlmt_ca_ps_v1_synthetic_gate_v3/status.txt" && \
 grep -qx "RLMT_CA_PS_LARGE_GATE_EXIT=0" \
  "results/rlmt_ca_ps_v1_synthetic_gate_v3/status.txt" || {
 echo "STOP: last full synthetic test gate was not completed" >&2
 exit 2
}
echo "=== 1. PARSE R, PYTHON AUDIT ==="
Rscript --vanilla -e '
 parse(file="R/experimental/rlmt_ca_ps_v1.R")
 parse(file="R/experimental/rlmt_ca_ps_v1_retro_runner.R")
' > /dev/null
python3 -m py_compile \
 scripts/rlmt_ca_ps_v1_retrospective_preflight.py \
 scripts/rlmt_ca_ps_v1_retro_independent_audit.py

echo "=== 2. STRICT PAPER-SHA INPUT + 285/95 REFERENCE PREFLIGHT ==="
python3 scripts/rlmt_ca_ps_v1_retrospective_preflight.py \
 --reference-repo "$REF" \
 --expected-sha "$HEAD"

mkdir "$NEW"
printf "%s\n" "$HEAD" > "$NEW/CODE_SHA.txt"
printf "%s\n" "$MAIN" > "$NEW/PAPER_MAIN_SHA.txt"
printf "%s\n" "$OLD" > "$NEW/ORIGINAL_REFERENCE_SHA.txt"
printf "%s\n" "EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY" \
 > "$NEW/SCIENTIFIC_STAGE.txt"
on_exit() {
 local rc=$?
 printf "RLMT_CA_PS_RETRO_RUN_EXIT=%s\n" "$rc" > "$NEW/status.txt"
}
trap on_exit EXIT

echo "=== 3. ALL 95 NEW CA-PS FITS (NO OLD FIT RETRAINING) ==="
Rscript --vanilla \
 R/experimental/rlmt_ca_ps_v1_retro_runner.R "$REF" \
 > "$NEW/full_95_run.log" 2>&1 || {
  echo "STOP: R fitting terminated; preserving all partial evidence"
  tail -n 90 "$NEW/full_95_run.log"
  exit 1
 }
tail -n 36 "$NEW/full_95_run.log"
grep -qx 'RLMT_CA_PS_RETRO_19X5_RUN=PASS' "$NEW/full_95_run.log"
printf "%s\n" "RLMT_CA_PS_RETRO_RUN_EXIT=0" > "$NEW/status.txt"

echo "=== 4. INDEPENDENT 380-MODEL MATCHED METRIC/PATIENT AUDIT ==="
python3 scripts/rlmt_ca_ps_v1_retro_independent_audit.py \
 --reference-repo "$REF" \
 > "$NEW/independent_audit.log" 2>&1 || {
  echo "STOP: independent audit failed; model results preserved"
  tail -n 90 "$NEW/independent_audit.log"
  exit 1
 }
cat "$NEW/independent_audit.log"
grep -qx 'RLMT_CA_PS_RETRO_INDEPENDENT_19X5_AUDIT=PASS' \
 "$NEW/independent_audit.log"

echo "RLMT_CA_PS_RETRO_19X5_COMPLETE=PASS"
echo "REPORT=$ROOT/$NEW/four_method_retro_audit_v1"
