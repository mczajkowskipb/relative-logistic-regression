#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export RLR_PROJECT_ROOT="$ROOT"
export RLR_PTB_FIXED_ROOT="${RLR_PTB_FIXED_ROOT:-$ROOT/results/ptb_fixed_published_state}"
export RLR_PTB_TRAIN_FILE="${RLR_PTB_TRAIN_FILE:-$ROOT/data/cohort_norm/01_train_data.arff}"
export RLR_MSIGDB_GMT_DIR="${RLR_MSIGDB_GMT_DIR:-$ROOT/reference_outputs/ptb_ora/gmt}"
export RLR_OUT_ROOT="${RLR_OUT_ROOT:-$ROOT/results/supplementary/ptb_fixed_ora}"
Rscript R/new/27_ptb_fixed_biological_case_study_trainonly.R
