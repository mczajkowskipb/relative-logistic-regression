#!/usr/bin/env bash
set -euo pipefail

TARGET_URL="https://github.com/mczajkowskipb/relative-logistic-regression.git"

if [[ ! -f README.md || ! -f LICENSE || ! -f CITATION.cff ]]; then
  echo "Run this script from the root of the clean publication release." >&2
  exit 2
fi

if [[ -d .git ]]; then
  echo "A .git directory already exists. Refusing to overwrite existing Git history." >&2
  echo "If this is already the clean repository, push it manually after verification." >&2
  exit 3
fi

python3 scripts/check_release.py --quick

git init
git checkout -b main
git add .
git commit -m "Publication release: Relative Logistic Regression"
git remote add origin "$TARGET_URL"
git push -u origin main

echo
echo "Published to: https://github.com/mczajkowskipb/relative-logistic-regression"
echo "Next: fresh-clone verification, then make relative-logistic-regression-working private."
