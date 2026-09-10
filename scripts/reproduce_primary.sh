#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export RLR_PROJECT_ROOT="$ROOT"

for f in R/new/07_world_protocol_utils.R R/new/08_run_worlds_config.R R/new/09_collect_worlds_config.R; do
  [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 2; }
done

run_world () {
  local cfg="$1"
  Rscript R/new/08_run_worlds_config.R "$cfg"
  Rscript R/new/09_collect_worlds_config.R "$cfg"
}

echo "[1/6] CompCancer 9x10 shard"
Rscript R/new/08_run_worlds_config.R configs/paper_frozen/compcancer_9x10.json

echo "[2/6] CompCancer West/Yeoh 2x2 shard"
Rscript R/new/08_run_worlds_config.R configs/paper_frozen/compcancer_2x2.json

echo "[3/6] CompCancer merged collection + reports"
Rscript R/new/09_collect_worlds_config.R configs/paper_frozen/compcancer_collect.json
Rscript R/new/10_make_paper_tables.R configs/paper_frozen/compcancer_collect.json
Rscript R/new/11_make_paper_figures.R configs/paper_frozen/compcancer_collect.json
Rscript R/new/12_make_stats_tests.R configs/paper_frozen/compcancer_collect.json

echo "[4/6] PTB fixed discovery-data split"
run_world configs/paper_frozen/ptb_fixed.json
Rscript R/new/10_make_paper_tables.R configs/paper_frozen/ptb_fixed.json

echo "[5/6] PTB LOCO"
run_world configs/paper_frozen/ptb_loco.json
Rscript R/new/10_make_paper_tables.R configs/paper_frozen/ptb_loco.json

echo "[6/6] legacy8"
run_world configs/paper_frozen/legacy8.json
Rscript R/new/10_make_paper_tables.R configs/paper_frozen/legacy8.json

echo "Primary reconstruction complete."
