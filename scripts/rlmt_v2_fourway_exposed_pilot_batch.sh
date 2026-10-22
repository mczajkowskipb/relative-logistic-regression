#!/usr/bin/env bash
# Adds only 12 NEW Gini and LA paper-cohort fits; prior H1/H2 read only.
set -Eeuo pipefail
test "$#" -eq 3 || { echo "ARGS: ORIGINAL_REF PRIOR_H1H2_WT PRIOR_LA_SYNTH_WT"; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
REF="$(cd "$1" && pwd -P)"
PREV="$(cd "$2" && pwd -P)"
LA="$(cd "$3" && pwd -P)"
cd "$ROOT"
HEAD="$(git rev-parse HEAD)"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
REFSHA="4873878f57f806bf5b52144756f8087b17fa540c"
PREVSHA="363e8c8e00718d2f834225d8434292bb408d6b3b"
LASHA="3abaf501f8f6303e14454e239af74e3fbce04bb7"
test "$(git rev-parse origin/main)" = "$PAPER"
test "$(git rev-parse origin/feature/rlmt-paper-fourway-pilot)" = "$HEAD"
test "$(git -C "$REF" rev-parse HEAD)" = "$REFSHA"
test "$(git -C "$PREV" rev-parse HEAD)" = "$PREVSHA"
test "$(git -C "$LA" rev-parse HEAD)" = "$LASHA"
test -z "$(git -c core.filemode=false status --porcelain --untracked-files=no)"
grep -qx 'RUNNER_EXIT=0' "$REF/results/rlmt_dc_published_data_5fold_v1/status.txt"
grep -qx 'RLMT_V2_PILOT_EXIT=0' "$PREV/results/rlmt_v2_exposed_pilot_2x3_v1/status.txt"
grep -qx 'RLMT_V2_PILOT_INDEPENDENT_AUDIT=PASS' "$PREV/results/rlmt_v2_exposed_pilot_2x3_v1/independent_audit.log"
grep -qx 'RLMT_LA_SYNTHETIC_GATE_EXIT=0' "$LA/results/rlmt_v2_gini_lookahead_synthetic_gate_v1/status.txt"
grep -qx 'RLMT_V2_GINI_LA_SYNTHETIC_REAL_GATE=PASS' "$LA/results/rlmt_v2_gini_lookahead_synthetic_gate_v1/real_synthetic.log"
for cmd in Rscript python3 sha256sum; do command -v "$cmd" >/dev/null; done
OUT="results/rlmt_v2_fourway_pilot_2x3_v1"
test ! -e "$OUT" || { echo "STOP: preserve prior evidence"; exit 2; }

echo "=== 1. PARSE CODE ==="
Rscript --vanilla -e 'parse("R/experimental/rlmt_v2_fourway_exposed_pilot_2x3.R");parse("R/experimental/rlmt_v2_gini_lookahead.R")' >/dev/null
python3 -m py_compile scripts/rlmt_v2_fourway_exposed_pilot_audit.py
mkdir -p "$OUT"
printf '%s\n' "$HEAD" > "$OUT/CODE_SHA.txt"
printf '%s\n' "$PAPER" > "$OUT/PAPER_SHA.txt"
printf '%s\n' "$REFSHA" > "$OUT/REFERENCE_SHA.txt"
printf '%s\n' "Two cohorts, first three exact original folds, TSPDT_GINI/RLMT_LA" > "$OUT/PREDECLARED.txt"
on_exit() {
 rc=$?
 printf 'RLMT_V2_FOURWAY_PILOT_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT

echo "=== 2. NEW REAL GINI AND LA FITS ==="
Rscript --vanilla R/experimental/rlmt_v2_fourway_exposed_pilot_2x3.R "$REF" "$PREV" > "$OUT/real_biomedical.log" 2>&1 || {
 tail -n 100 "$OUT/real_biomedical.log"; exit 1;
}
tail -n 30 "$OUT/real_biomedical.log"
grep -qx 'RLMT_V2_FOURWAY_RUN=PASS' "$OUT/real_biomedical.log"

echo "=== 3. INDEPENDENT 4-WAY PAIRED AUDIT ==="
python3 scripts/rlmt_v2_fourway_exposed_pilot_audit.py --reference-repo "$REF" --previous-pilot "$PREV" > "$OUT/independent_audit.log" 2>&1 || {
 cat "$OUT/independent_audit.log"; exit 1;
}
cat "$OUT/independent_audit.log"
grep -qx 'RLMT_V2_FOURWAY_INDEPENDENT_AUDIT=PASS' "$OUT/independent_audit.log"
echo "=== 4. WARNINGS TO REVIEW ==="
grep -Ei 'Warning|nonconver|deprecated' "$OUT/real_biomedical.log" > "$OUT/warnings.txt" || true
tail -n 50 "$OUT/warnings.txt"
echo "RLMT_V2_FOURWAY_PAPER_PILOT_COMPLETE=PASS"
echo "SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY"
echo "REPORT=$ROOT/$OUT"
