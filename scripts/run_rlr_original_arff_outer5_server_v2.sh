#!/usr/bin/env bash
# READ-ONLY on original ARFF; only NEW benchmark outputs inside current checkout.
# Usage: bash scripts/run_rlr_original_arff_outer5_server_v2.sh /opt/2026/RLR/data
set -euo pipefail
cd "$(dirname "$0")/.."
DATA_ROOT="${1:-/opt/2026/RLR/data}"
export RLR_ORIGINAL_DATA_ROOT="$DATA_ROOT"
if [[ ! -d "$DATA_ROOT/CompCancer" ]]; then
  echo "BLOCKED: missing original CompCancer input: $DATA_ROOT/CompCancer" >&2; exit 2
fi
for executable in Rscript python3 git sha256sum; do
  command -v "$executable" >/dev/null || { echo "BLOCKED: missing $executable" >&2; exit 2; }
done
Rscript --vanilla -e 'pkg<-c("foreign","switchBox","glmnet","ncvreg","pROC","jsonlite");bad<-pkg[!vapply(pkg,requireNamespace,logical(1),quietly=TRUE)];if(length(bad)){cat("BLOCKED_MISSING_R_PACKAGES:",paste(bad,collapse=","),"\n");quit(status=2)}; cat("R_DEPENDENCIES=PASS\n")'
Rscript --vanilla -e 'parse("R/experimental/rlr_binary_outer5_original_arff_v2.R");cat("ORIGINAL_RUNNER_R_PARSE=PASS\n")' >/dev/null
python3 -m py_compile scripts/audit_rlr_binary_outer5_v2.py
python3 scripts/rlr_original_compcancer_sha_preflight_v2.py --data-root "$DATA_ROOT"
out="results/rlr_binary_outer5_original_arff_v2"
mkdir -p "$out"
for ds in alizadeh-2000-v1 armstrong-2002-v1 bittner-2000 chen-2002 golub-1999-v1 gordon-2002 pomeroy-2002-v1 shipp-2002-v1 singh-2002 west-2001 yeoh-2002-v1; do
  echo "RUN_ORIGINAL_ARFF_OUTER5 dataset=$ds"
  mkdir -p "$out/$ds"
  git rev-parse HEAD > "$out/$ds/SOURCE_SHA.txt"
  git hash-object R/new/07_world_protocol_utils.R > "$out/$ds/ORIGINAL_RLR_CORE_SHA.txt"
  RLR_DATASET_ID="$ds" Rscript --vanilla R/experimental/rlr_binary_outer5_original_arff_v2.R 2>&1 | tee "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_ORIGINAL_ARFF_SHA256_CONFIRMED=PASS' "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_OUTER5_FIVE_FOLDS_THREE_METHODS=PASS' "$out/$ds/fit.log"
  grep -qx 'RLR_BINARY_OUTER5_ALL_PATIENTS_EXACTLY_ONCE_PER_METHOD=PASS' "$out/$ds/fit.log"
done
audit="results/rlr_binary_outer5_original_arff_audit_v2"
python3 scripts/audit_rlr_binary_outer5_v2.py --input "$out" --out "$audit" | tee "$audit.log"
grep -qx 'RLR_BINARY_OUTER5_165_FOLD_METRICS_INDEPENDENT_REPLAY=PASS' "$audit.log"
grep -qx 'RLR_BINARY_OUTER5_EVERY_PATIENT_EXACTLY_ONCE_PER_METHOD=PASS' "$audit.log"
grep -qx 'RLR_BINARY_OUTER5_EXACT_2048_SIGNFLIPS_AND_BH8=PASS' "$audit.log"
echo "RLR_ORIGINAL_ARFF_OUTER5_ALL11_RUN_AND_INDEPENDENT_AUDIT=PASS"
echo "RESULTS=$audit"
