#!/usr/bin/env python3
"""Read-only SHA256 preflight for the exact processed ARFF/data snapshot.

Usage (from any directory):
  python3 scripts/verify_original_arff_preflight.py --data-root /path/to/RLR/data

Inputs are NOT downloaded and expression values are NOT exported.
A SHA256 match proves byte-for-byte identity, including original ARFF row/column
order, probe names, class attributes, and class-level ordering.
"""
import argparse
import hashlib
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
PAPER_DATASETS = (
    "alizadeh-2000-v1", "armstrong-2002-v1", "bittner-2000",
    "chen-2002", "golub-1999-v1", "gordon-2002", "pomeroy-2002-v1",
    "shipp-2002-v1", "singh-2002", "west-2001", "yeoh-2002-v1",
)


def digest(path):
    state = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            state.update(block)
    return state.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=pathlib.Path, default=REPO / "data",
                    help="Directory containing CompCancer/, CompCancer_truth_*/ and other frozen inputs")
    ap.add_argument("--manifest", type=pathlib.Path,
                    default=REPO / "checksums/DATA_SHA256_EXPECTED.txt",
                    help="Original frozen SHA256 manifest")
    args = ap.parse_args()
    root = args.data_root.expanduser().resolve()
    if not root.is_dir():
        ap.error(f"Data root is not a directory: {root}")
    if not args.manifest.is_file():
        ap.error(f"Missing expected checksum manifest: {args.manifest}")

    entries = []
    for number, line in enumerate(args.manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64})\s+(data/[^\s]+)", line.strip())
        if not match:
            ap.error(f"Invalid SHA256 manifest at line {number}")
        expected, logical_path = match.groups()
        relative = pathlib.PurePosixPath(logical_path)
        if ".." in relative.parts or relative.is_absolute():
            ap.error(f"Unsafe manifest path: {logical_path}")
        entries.append((expected.lower(), logical_path, root.joinpath(*relative.parts[1:])))
    if not entries or len({name for _, name, _ in entries}) != len(entries):
        ap.error("Empty manifest or duplicate paths")

    manifest_paths = {name for _, name, _ in entries}
    for dataset in PAPER_DATASETS:
        primary = f"data/CompCancer/{dataset}.arff"
        shard = ("CompCancer_truth_2x2" if dataset in ("west-2001", "yeoh-2002-v1")
                 else "CompCancer_truth_9x10")
        secondary = f"data/{shard}/{dataset}.arff"
        if primary not in manifest_paths or secondary not in manifest_paths:
            ap.error(f"Canonical original ARFF/shard missing from manifest: {dataset}")

    matched, missing, mismatched = 0, [], []
    for expected, name, path in entries:
        if not path.is_file():
            missing.append(name)
            print(f"MISSING {name}")
            continue
        observed = digest(path)
        if observed != expected:
            mismatched.append(name)
            print(f"MISMATCH {name} expected={expected} actual={observed}")
        else:
            matched += 1
            print(f"OK {name}")
    print(f"SHA256_MANIFEST_MATCHED={matched}/{len(entries)}")
    if missing or mismatched:
        print(f"SHA256_MANIFEST_STATUS=FAIL missing={len(missing)} mismatch={len(mismatched)}")
        return 2
    print(f"ORIGINAL_COMPCANCER_ARFF_SOURCE_AND_SHARD=PASS 11/11")
    print(f"SHA256_MANIFEST_STATUS=PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
