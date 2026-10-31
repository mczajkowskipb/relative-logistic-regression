#!/usr/bin/env bash
# Reproduce Supplement S7.1 on the ORIGINAL processed CompCancer matrices.
# Run as a child process: bash scripts/reproduce_compcancer_ternary_binary_ablation.sh /path/to/data
# Never write to the source data directory.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
DATA_ROOT="${1:-$ROOT/data}"
export RLR_ORIGINAL_DATA_ROOT="$DATA_ROOT"
for cmd in python3 Rscript git sha256sum; do
  command -v "$cmd" >/dev/null || { echo "MISSING_EXECUTABLE=$cmd" >&2; exit 2; }
done
Rscript --vanilla -e '
  req <- c("foreign","switchBox","glmnet","ncvreg","pROC","jsonlite")
  absent <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if(length(absent))stop("Missing R packages: ",paste(absent,collapse=", "))
  cat("S7_1_R_DEPENDENCIES=PASS\n")
'
Rscript --vanilla -e 'parse("R/new/28_compcancer_ternary_vs_binary_ablation.R");cat("S7_1_SCRIPT_PARSE=PASS\n")' >/dev/null
python3 scripts/verify_original_arff_preflight.py --data-root "$DATA_ROOT"
out="results/rlr_binary_outer5_original_arff_v2"
mkdir -p "$out"
for ds in \
  alizadeh-2000-v1 armstrong-2002-v1 bittner-2000 chen-2002 \
  golub-1999-v1 gordon-2002 pomeroy-2002-v1 shipp-2002-v1 \
  singh-2002 west-2001 yeoh-2002-v1
do
  mkdir -p "$out/$ds"
  if [[ -s "$out/$ds/fold_metrics.csv" && -s "$out/$ds/patients.csv" && -s "$out/$ds/fold_manifest.csv" ]]; then
    echo "S7_1_EXISTING_COHORT_PENDING_AUDIT=$ds"
    continue
  fi
  echo "S7_1_FITTING_COHORT=$ds"
  git rev-parse HEAD > "$out/$ds/SOURCE_SHA.txt"
  git hash-object R/new/07_world_protocol_utils.R > "$out/$ds/ORIGINAL_RLR_CORE_SHA.txt"
  RLR_DATASET_ID="$ds" Rscript --vanilla R/new/28_compcancer_ternary_vs_binary_ablation.R \
    2>&1 | tee "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_ORIGINAL_ARFF_SHA256_CONFIRMED=PASS' "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_OUTER5_FIVE_FOLDS_THREE_METHODS=PASS' "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_OUTER5_ALL_PATIENTS_EXACTLY_ONCE_PER_METHOD=PASS' "$out/$ds/fit.log"
done
audit="results/rlr_binary_outer5_original_arff_audit_v2"
mkdir -p "$audit"
python3 scripts/audit_compcancer_ternary_binary_ablation.py \
  --input "$out" --out "$audit" --original-data-root "$DATA_ROOT" |
  tee "$audit/audit.log"
grep -qx 'RLR_BINARY_OUTER5_165_FOLD_METRICS_INDEPENDENT_REPLAY=PASS' "$audit/audit.log"
grep -qx 'RLR_BINARY_OUTER5_EVERY_PATIENT_EXACTLY_ONCE_PER_METHOD=PASS' "$audit/audit.log"
grep -qx 'RLR_BINARY_OUTER5_EXACT_2048_SIGNFLIPS_AND_BH8=PASS' "$audit/audit.log"
grep -qx 'ORIGINAL_PAPER_ARFF_SHA256_PARITY=PASS' "$audit/audit.log"
grep -qx 'RLR_ORIGINAL_ARFF_11_SOURCE_AND_SHARD_SHA256_RECHECK=PASS' "$audit/audit.log"
echo "S7_1_ORIGINAL_ARFF_FULL11_EXTERNAL_AUDIT=PASS"
echo "S7_1_OUTPUT_DIR=$audit"
