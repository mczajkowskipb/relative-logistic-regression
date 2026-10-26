#!/usr/bin/env bash
# Read-only scientific audit; create only a fresh own result directory.
set -Eeuo pipefail
test "$#" -eq 1 || { echo "H0_STOP: StageG checkout needed" >&2; exit 2; }
ROOT="$(git rev-parse --show-toplevel)"
G="$(cd "$1" && pwd -P)"
cd "$ROOT"
BRANCH="analysis/crit-stageH0-symmetry-audit-v1"
G_SHA="6ca44c30a439dbdb97b3a45f7e8f45432da3bd59"
echo "H0_AUDIT_SOURCE_SHA=$(git rev-parse HEAD)"
test "$(git rev-parse HEAD)" = \
     "$(git rev-parse "refs/remotes/origin/$BRANCH")" || {
  echo "H0_STOP: source SHA not pinned to remote" >&2; exit 2;
}
test "$(git -C "$G" rev-parse HEAD)" = "$G_SHA"
grep -qx "STAGEG_BATCH_EXIT=0" \
 "$G/results/crit_v9_stageG_exhaustive32_v1/status.txt"
for tree in "$G" "$ROOT"; do
 test -z "$(git -C "$tree" status --porcelain --untracked-files=no)" || {
  echo "H0_STOP: tracked source modified" >&2; exit 2;
 }
done
python3 - <<'PY'
from pathlib import Path
p=Path("scripts/crit_stageH0_symmetry_audit.py")
compile(p.read_text(),str(p),"exec")
print("H0_SCRIPT_PARSE=PASS")
PY
OUT="results/crit_stageH0_symmetry_v1"
test ! -e "$OUT" && test ! -L "$OUT" || {
 echo "H0_STOP: output path already exists; no overwrite" >&2; exit 2;
}
if [ "${H0_PREFLIGHT_ONLY:-0}" = 1 ]; then
 echo "H0_PREFLIGHT=PASS"
 exit 0
fi
mkdir -p "$OUT"
printf '%s\n' "$(git rev-parse HEAD)" > "$OUT/H0_SOURCE_SHA.txt"
printf '%s\n' "$G_SHA" > "$OUT/STAGEG_SOURCE_SHA.txt"
on_exit(){
 rc=$?
 printf 'H0_BATCH_EXIT=%s\n' "$rc" > "$OUT/status.txt"
}
trap on_exit EXIT
python3 scripts/crit_stageH0_symmetry_audit.py \
 --stageg "$G/results/crit_v9_stageG_exhaustive32_v1" \
 --out "$ROOT/$OUT/symmetry_report" \
 > "$OUT/audit.log" 2>&1 || {
 tail -n 80 "$OUT/audit.log" >&2; exit 1;
}
cat "$OUT/audit.log"
grep -qx "H0_SYMMETRY_EXACT_12_OF_12=PASS" "$OUT/audit.log"
grep -qx "H0_STAGEG_FROZEN_34560_PATIENT_AUDIT_VERIFIED=PASS" "$OUT/audit.log"
echo "H0_ALL_GATES=PASS"
