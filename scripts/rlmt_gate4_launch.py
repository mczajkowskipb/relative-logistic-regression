#!/usr/bin/env python3
"""GATE-4 Gordon launcher: checks provenance before creating results.
Run only after fast-forwarding feature/rlmt-prototype:
  python3 scripts/rlmt_gate4_launch.py --expected-sha <40-hex-git-sha>
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
MAIN = "4079a5ac400f1d2cfded523c7000f655470925ce"
BRANCH = "feature/rlmt-prototype"
CONFIG = "configs/rlmt/gate4_gordon_frozen.json"
RESULT = "results/rlmt_gate4_gordon_frozen_5fold"
CANONICAL_PATH = "data/CompCancer/gordon-2002.arff"
SCREEN = "results/rlmt_gate4_metadata_screen_v1/metadata_inventory.csv"
HASH_MANIFEST = "checksums/DATA_SHA256_EXPECTED.txt"


def stop(message: str) -> None:
    raise SystemExit("STOP: " + message)


def git(*args: str) -> str:
    p = subprocess.run(
        ["git", *args], cwd=ROOT, check=True,
        capture_output=True, text=True
    )
    return p.stdout.strip()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def digest_from_manifest(text: str) -> str:
    hits = []
    for line in text.splitlines():
        m = re.fullmatch(r"([0-9a-fA-F]{64})\s+(.+)", line.strip())
        if m and m.group(2) == CANONICAL_PATH:
            hits.append(m.group(1).lower())
    if len(hits) != 1:
        stop("canonical publication manifest must contain one Gordon SHA-256")
    return hits[0]


def run() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--expected-sha", required=True)
    args = ap.parse_args()
    expected_sha = args.expected_sha.lower()
    if not re.fullmatch("[a-f0-9]{40}", expected_sha):
        stop("expected Git commit must be a full 40-character SHA")
    if Path.cwd().resolve() != ROOT:
        stop("run from repository root")
    if git("branch", "--show-current") != BRANCH:
        stop("incorrect branch")
    if git("status", "--porcelain", "--untracked-files=no"):
        stop("tracked working tree is dirty")
    actual_commit = git("rev-parse", "HEAD")
    if actual_commit != expected_sha:
        stop(f"wrong HEAD: {actual_commit}")
    if git("rev-parse", "origin/feature/rlmt-prototype") != expected_sha:
        stop("local HEAD differs from fetched feature branch")
    if git("rev-parse", "origin/main") != MAIN:
        stop("frozen publication main changed")
    if not (ROOT / "docs/RLMT_GATE4_GORDON_FREEZE.md").is_file():
        stop("missing pre-registered GATE-4 scientific contract")

    config_path = ROOT / CONFIG
    if not config_path.is_file():
        stop("missing frozen GATE-4 config")
    cfg = json.loads(config_path.read_text())
    world = cfg["worlds"]
    models = cfg["models"]
    if (cfg["out_root"] != RESULT or len(world) != 1 or
        world[0]["protocol"] != "compcancer_cv" or
        world[0]["datasets_include"] != ["gordon-2002"] or
        (world[0]["n_folds"], world[0]["n_repeats"],
         world[0]["max_splits_per_dataset"]) != (5, 1, 5) or
        models["seed"] != 42 or
        (models["run_baselines"], models["run_rlr"], models["run_rlmt"],
         models["skip_existing"]) != (False, True, True, False) or
        models["top_pool_values"] != [1000] or
        models["rlr_grid"] != {"K": [10], "m": [10], "d": [2]} or
        models["rlmt"] != {
            "K": 10, "m": 10, "d": 2, "pool": 1000,
            "max_leaves": [1, 2], "max_depth": 2,
            "screen_top_B": 20, "min_child_n": 16,
            "min_child_class": 6, "inner_folds": 3,
            "alpha": 0.5, "lambda_choice": "lambda.min", "seed": 42,
        }):
        stop("configuration changed from frozen GATE-4 Gordon protocol")

    screen_path = ROOT / SCREEN
    if not screen_path.is_file():
        stop("missing completed, performance-blind metadata screen")
    with screen_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if (len(rows) != 9 or
        not all(r["sha256_matches"] == "TRUE" and
                r["parse_ok"] == "TRUE" for r in rows)):
        stop("metadata screen integrity did not pass")
    eligible = [r["dataset"] for r in rows if r["candidate_eligible"] == "TRUE"]
    if eligible != ["gordon-2002"]:
        stop(f"unexpected eligible dataset set: {eligible}")

    # Both HEAD and frozen publication main must agree on this checksum.
    published_manifest = git("show", f"{MAIN}:{HASH_MANIFEST}")
    expected_digest = digest_from_manifest(published_manifest)
    current_manifest = (ROOT / HASH_MANIFEST).read_text()
    if digest_from_manifest(current_manifest) != expected_digest:
        stop("working checksum manifest differs from publication main")
    gordon_rows = [r for r in rows if r["dataset"] == "gordon-2002"]
    if len(gordon_rows) != 1:
        stop("metadata inventory must have exactly one Gordon entry")
    if gordon_rows[0]["sha256"].lower() != expected_digest:
        stop("metadata inventory Gordon digest differs from publication main")

    datafile = Path(world[0]["comp_dir"]) / "gordon-2002.arff"
    if not datafile.is_file():
        stop(f"missing Gordon input: {datafile}")
    actual_digest = sha256_file(datafile)
    if actual_digest != expected_digest:
        stop("Gordon input SHA256 mismatch versus frozen publication main")

    out = ROOT / RESULT
    if out.exists():
        stop(f"GATE-4 result directory already exists: {out}")

    print("=== R PARSER ===", flush=True)
    subprocess.run([
        "Rscript", "--vanilla", "-e",
        'parse(file="R/new/08_run_worlds_config.R");'
        'parse(file="R/new/28_rlmt_utils.R");'
        'parse(file="R/new/38_rlmt_gate4_check.R")',
    ], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    print("R_PARSE=PASS", flush=True)

    # No result directory is created until all checks above pass.
    out.mkdir(parents=True)
    (out / "CODE_SHA.txt").write_text(expected_sha + "\n")
    (out / "SOURCE_SHA256.txt").write_text(
        f"{expected_digest}  {datafile}\n"
    )
    (out / "FROZEN_CONFIG_SHA256.txt").write_text(
        sha256_file(config_path) + "\n"
    )
    script = out / "run.sh"
    script.write_text(f"""#!/usr/bin/env bash
cd '{ROOT}' || exit 2
if [ "$(git rev-parse HEAD)" != "{expected_sha}" ] || \\
   [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "STOP: GATE4 code changed after preflight" >&2
    echo "RUNNER_EXIT=90" > '{out}/status.txt'
    exit 90
fi
Rscript --vanilla R/new/08_run_worlds_config.R \\
  '{CONFIG}' > '{out}/run.log' 2>&1
rc=$?
printf 'RUNNER_EXIT=%s\\n' "$rc" > '{out}/status.txt'
exit "$rc"
""")

    with (out / "launcher.log").open("w") as logfile:
        process = subprocess.Popen(
            ["bash", str(script)],
            cwd=ROOT,
            stdin=subprocess.DEVNULL,
            stdout=logfile,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    (out / "runner.pid").write_text(str(process.pid) + "\n")
    print("GATE4_PREFLIGHT=PASS")
    print(f"CODE_SHA={expected_sha}")
    print(f"GORDON_SHA256={expected_digest}")
    print("DATASETS=gordon-2002 OUTER_FOLDS=5")
    print(f"GATE4_STARTED PID={process.pid}")
    print(f"RESULTS={out}")


if __name__ == "__main__":
    try:
        run()
    except subprocess.CalledProcessError as exc:
        stop(f"preflight command failed ({exc.cmd}): exit {exc.returncode}")
