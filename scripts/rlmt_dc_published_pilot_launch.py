#!/usr/bin/env python3
"""Launch frozen exploratory 19-publication-dataset smoke, after integrity checks.

No new GEO, no automatic cohort acquisition, no parameter tuning.
Run: python3 scripts/rlmt_dc_published_pilot_launch.py --expected-sha <sha>
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
BRANCH = "feature/rlmt-dc-v0"
PAPER_MAIN = "4079a5ac400f1d2cfded523c7000f655470925ce"
CFG = "configs/rlmt_dc/published_data_smoke_v1.json"
MANIFEST = "checksums/DATA_SHA256_EXPECTED.txt"
OUT = "results/rlmt_dc_published_data_smoke_v1"


def abort(message: str) -> None:
    raise SystemExit("PUBLICATION_PILOT_STOP: " + message)


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True,
        text=True, check=True
    ).stdout.strip()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected-sha", required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()

    if Path.cwd().resolve() != ROOT:
        abort("run from RLR repository root")
    expected_commit = args.expected_sha.lower()
    if not re.fullmatch(r"[a-f0-9]{40}", expected_commit):
        abort("expected SHA must have 40 hexadecimal digits")
    if git("branch", "--show-current") != BRANCH:
        abort("incorrect Git branch")
    if git("status", "--porcelain", "--untracked-files=no"):
        abort("tracked working tree is dirty")
    if git("rev-parse", "HEAD") != expected_commit:
        abort("current code differs from pinned SHA")
    if git("rev-parse", "origin/" + BRANCH) != expected_commit:
        abort("remote feature ref differs from pinned SHA")
    if git("rev-parse", "origin/main") != PAPER_MAIN:
        abort("publication main changed: stop")
    if (ROOT / OUT).exists():
        abort("pilot output already exists; do not overwrite results")

    cfg_path = ROOT / CFG
    cfg = json.loads(cfg_path.read_text())
    if (cfg.get("stage") != "PUBLICATION_DATA_EXPLORATORY_SMOKE_V1"
        or cfg.get("out_root") != OUT
        or (cfg.get("model_seed"), cfg.get("n_folds"),
            cfg.get("n_repeats"), cfg.get("max_splits_per_dataset"))
           != (42, 5, 1, 1)
        or cfg.get("models") != {
            "RLR": "K10_m10_d2_pool1000_alpha0.5_lambda.min",
            "RLMT_v1": "unchanged_L2",
            "RLMT_DC_v0": "discovery_confirmation_4096",
        }):
        abort("frozen exploratory pilot model/split config mismatch")
    expected_groups = [
        ("compcancer_9x10", "CompCancer_truth_9x10", "CompCancer",
         ["alizadeh-2000-v1", "armstrong-2002-v1", "bittner-2000",
          "chen-2002", "golub-1999-v1", "gordon-2002",
          "pomeroy-2002-v1", "shipp-2002-v1", "singh-2002"]),
        ("compcancer_2x2", "CompCancer_truth_2x2", "CompCancer",
         ["west-2001", "yeoh-2002-v1"]),
        ("legacy8", "comparison", "comparison",
         ["adenocarcinoma", "breast2", "cns", "colon",
          "dlbcl", "leukemia", "ovarian", "prostate"]),
    ]
    actual_groups = [
        (x["world"], x["data_subdir"], x["manifest_subdir"], x["datasets"])
        for x in cfg.get("cohorts", [])
    ]
    if actual_groups != expected_groups:
        abort("cohort inventory differs from publication preregistration")
    if not (ROOT / "docs/RLMT_DC_PUBLISHED_DATA_PROTOCOL.md").is_file():
        abort("must register scientific protocol before launch")

    paper_manifest = git("show", PAPER_MAIN + ":" + MANIFEST)
    expected_hashes: dict[str, str] = {}
    for line in paper_manifest.splitlines():
        match = re.fullmatch(r"([a-fA-F0-9]{64})\s+(\S+)", line.strip())
        if not match:
            abort("invalid publication SHA256 manifest entry")
        if match.group(2) in expected_hashes:
            abort("duplicate publication manifest path")
        expected_hashes[match.group(2)] = match.group(1).lower()
    current_manifest = (ROOT / MANIFEST).read_text()
    if current_manifest.strip() != paper_manifest.strip():
        abort("working checksum manifest differs from publication main")

    files: list[dict[str, str]] = []
    for world, subdir, reference_subdir, datasets in expected_groups:
        for dataset in datasets:
            key = f"data/{reference_subdir}/{dataset}.arff"
            sha = expected_hashes.get(key)
            if not sha:
                abort(f"missing canonical publication hash for {key}")

            # Two documented publication-era on-server data layouts are
            # allowed, but actual bytes MUST match the exact paper SHA.
            paths = [
                ROOT / "data" / subdir / (dataset + ".arff"),
                ROOT.parent / "data" / subdir / (dataset + ".arff"),
            ]
            tried: list[str] = []
            good: list[Path] = []
            for path in paths:
                if not path.is_file():
                    tried.append(f"{path}:missing")
                    continue
                actual = digest(path)
                if actual == sha:
                    resolved = path.resolve()
                    if resolved not in good:
                        good.append(resolved)
                else:
                    tried.append(f"{path}:SHA_MISMATCH:{actual}")
            if not good:
                abort(f"{dataset} not found with paper SHA256; "
                      + "; ".join(tried))
            files.append({
                "world": world, "dataset": dataset,
                "source_path": str(good[0]), "sha256": sha,
                "reference_manifest_path": key,
            })
    if len(files) != 19:
        abort("not exactly 19 verified publication datasets")

    # Verify the actual R runner/audit parse BEFORE creating or
    # executing any result directory. The Python script itself can
    # be checked separately with py_compile.
    subprocess.run(
        [
            "Rscript", "--vanilla", "-e",
            'parse(file="R/new/44_rlmt_dc_published_pilot.R");'
            'parse(file="R/new/45_rlmt_dc_published_audit.R")',
        ],
        cwd=ROOT, check=True,
        stdout=subprocess.DEVNULL,
    )
    print("PUBLICATION_PILOT_R_PARSE=PASS")
    print("PUBLICATION_PILOT_PREFLIGHT=PASS")
    print("SHARDS=CompCancer:11 legacy8:8")
    print("EXPLORATORY_FIRST_FOLD_PER_DATASET=19")
    print("MODELS=RLR;RLMT_v1;RLMT_DC_v0")
    print("PLANNED_MODEL_FITS=57")
    print("PTB=DEFERRED_SEPARATE_PUBLISHED_PROTOCOL")
    print("NEW_GEO_DATA=NONE")
    if args.preflight_only:
        print("NO_RESULT_DIRECTORY_CREATED=TRUE")
        return

    # Output created only after complete, performance-blind integrity
    # screening succeeds for ALL 19 fixed inputs.
    out = ROOT / OUT
    out.mkdir(parents=True)
    with (out / "input_manifest.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(files[0]))
        writer.writeheader()
        writer.writerows(files)
    (out / "CODE_SHA.txt").write_text(expected_commit + "\n")
    (out / "PUBLICATION_MAIN_SHA.txt").write_text(PAPER_MAIN + "\n")
    (out / "FROZEN_CONFIG_SHA256.txt").write_text(
        digest(cfg_path) + "\n"
    )
    (out / "README.txt").write_text(
        "Exploratory method-development benchmark on ALL paper ARFF datasets "
        "CompCancer11 + legacy8. First fixed outer fold, seed42. "
        "Same outer split/outer-training global top1000 pool across "
        "RLR/RLMT-v1/RLMT-DC-v0. No new external patient cohorts. "
        "No confirmatory statistical inference or test-driven tuning. "
        "See docs/RLMT_DC_PUBLISHED_DATA_PROTOCOL.md.\n"
    )
    run = out / "run.sh"
    run.write_text(f"""#!/usr/bin/env bash
cd '{ROOT}' || exit 2
if [ "$(git rev-parse HEAD)" != "{expected_commit}" ] || \\
   [ "$(git branch --show-current)" != "{BRANCH}" ] || \\
   [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "STOP: changed Git code after preflight" >&2
  echo "RUNNER_EXIT=90" > '{out}/status.txt'
  exit 90
fi
Rscript --vanilla R/new/44_rlmt_dc_published_pilot.R '{CFG}' \\
  > '{out}/run.log' 2>&1
rc=$?
printf 'RUNNER_EXIT=%s\\n' "$rc" > '{out}/status.txt'
exit "$rc"
""")
    with (out / "launcher.log").open("w") as handle:
        process = subprocess.Popen(
            ["bash", str(run)], stdin=subprocess.DEVNULL,
            stdout=handle, stderr=subprocess.STDOUT,
            cwd=ROOT, start_new_session=True,
        )
    (out / "runner.pid").write_text(str(process.pid) + "\n")
    print("PUBLICATION_PILOT_STARTED_PID=" + str(process.pid))
    print("RUN_LOG=" + str(out / "run.log"))
    print("SUMMARY_MANIFEST=" + str(out / "run_manifest.csv"))


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        abort(f"preflight subprocess failed: {error}")
