#!/usr/bin/env python3
"""Fail-closed preflight for NEW RLMT-CA-PS-v1 matched retrospective 19x5 run.

Old RLR/RLMT results + publication input shards are read-only.
This checks the frozen source and recorded original training-test splits.
No output files and no model fits are created.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
PAPER_SHA="4079a5ac400f1d2cfded523c7000f655470925ce"
OLD_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
BRANCH="feature/rlmt-class-aware-ps-v1"
OLD_RESULTS="results/rlmt_dc_published_data_5fold_v1"
OUT="results/rlmt_ca_ps_v1_retrospective_19x5"
EXPECTED={
 "compcancer_9x10":("alizadeh-2000-v1","armstrong-2002-v1",
   "bittner-2000","chen-2002","golub-1999-v1",
   "gordon-2002","pomeroy-2002-v1","shipp-2002-v1","singh-2002"),
 "compcancer_2x2":("west-2001","yeoh-2002-v1"),
 "legacy8":("adenocarcinoma","breast2","cns","colon","dlbcl",
            "leukemia","ovarian","prostate"),
}


def stop(msg):
    raise SystemExit("RLMT_CA_PS_RETRO_PREFLIGHT_STOP: "+msg)


def git(path,*args):
    return subprocess.run(["git",*args],cwd=path,capture_output=True,
                          text=True,check=True).stdout.strip()


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):
            h.update(chunk)
    return h.hexdigest()


def csvrows(path):
    if not path.is_file(): stop("missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference-repo",required=True)
    p.add_argument("--expected-sha",required=True)
    a=p.parse_args()
    ref=Path(a.reference_repo).resolve(strict=True)
    if Path.cwd().resolve()!=ROOT: stop("must run from experimental worktree")
    try:
        if git(ROOT,"rev-parse","HEAD")!=a.expected_sha:
            stop("new method commit differs from pinned expected SHA")
        if git(ROOT,"rev-parse","origin/"+BRANCH)!=a.expected_sha:
            stop("remote method commit differs from expected SHA")
        if git(ROOT,"rev-parse","origin/main")!=PAPER_SHA:
            stop("paper main changed")
        if git(ROOT,"-c","core.filemode=false","status","--porcelain",
               "--untracked-files=no"):
            stop("new method tracked worktree differs from checkout")
        if git(ref,"rev-parse","HEAD")!=OLD_SHA or (
            git(ref,"-c","core.filemode=false","status",
                "--porcelain","--untracked-files=no")
        ):
            stop("original 285-fit reference worktree is not frozen/clean")
        if (ROOT/OUT).exists(): stop("new outcome root already exists")
    except subprocess.CalledProcessError as exc:
        stop("Git check failed "+str(exc))
    old=ref/OLD_RESULTS
    if not old.is_dir() or not (old/"status.txt").is_file():
        stop("original 285-fit benchmark not available")
    if (old/"status.txt").read_text().strip()!="RUNNER_EXIT=0":
        stop("original reference run did not complete")
    for requirement in (
        "published_5fold_audit_v1/independent_fold_metrics.csv",
        "published_5fold_accepted_replay_v1/accepted_models.csv",
        "input_manifest.csv","run_manifest.csv"
    ):
        if not (old/requirement).is_file():
            stop("old independent audit artifact absent: "+requirement)
    inputrows=csvrows(old/"input_manifest.csv")
    oldfits=csvrows(old/"run_manifest.csv")
    if len(inputrows)!=19 or len(oldfits)!=285:
        stop("original 19 datasets or 285 fits missing")
    exp=[(world,dataset) for world,ds in EXPECTED.items() for dataset in ds]
    if [(r["world"],r["dataset"]) for r in inputrows]!=exp:
        stop("input cohorts/order diverge from frozen 19-dataset reference")
    keys=set()
    for row in oldfits:
        key=(row["world"],row["dataset"],row["split_id"],row["method"])
        if key in keys: stop("duplicate old reference model "+str(key))
        keys.add(key)
        if row["status"]!="OK": stop("old fit failed "+str(key))
    all_fold_ids=tuple(f"rep00_fold{i:02d}" for i in range(5))
    for world,ds in exp:
        for fold in all_fold_ids:
            for method in ("RLR","RLMT_v1","RLMT_DC_v0"):
                if (world,ds,fold,method) not in keys:
                    stop("missing frozen original fold "+str((world,ds,fold,method)))
    for row in inputrows:
        path=Path(row["source_path"]).resolve(strict=True)
        if sha256(path).lower()!=row["sha256"].lower():
            stop("source hash mismatch "+row["dataset"])
        if not (len(row["sha256"])==64 and
                row["reference_manifest_path"].endswith("/"+row["dataset"]+".arff")):
            stop("source manifest schema invalid "+row["dataset"])
    print("RLMT_CA_PS_RETRO_PREFLIGHT=PASS")
    print("ORIGINAL_REFERENCE_FITS=285")
    print("PUBLISHED_INPUT_SHA256_VERIFIED=19")
    print("OUTER_FOLDS_MATCHED=95")
    print("NEW_MODEL_FITS_PLANNED=95")
    print("NO_OLD_RESULT_WRITES=TRUE")
    print("NO_RESULT_DIR_CREATED=TRUE")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY")


if __name__=="__main__":
    main()
