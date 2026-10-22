#!/usr/bin/env python3
"""Independent NO-WRITE preflight for four-method RLR/RLMT 19x5 screen.

All 19 cohorts are EXPOSED from previous RLR/RLMT/CA-PS development.
Read archived original 285 models without touching their predictions;
verify original input SHA256, manifest 95 fold locations and frozen git SHA.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OLD_SHA="4873878f57f806bf5b52144756f8087b17fa540c"
PAPER="4079a5ac400f1d2cfded523c7000f655470925ce"
OLD_DIR="results/rlmt_dc_published_data_5fold_v1"
NEW_DIR="results/rlmt_v2_19x5_four_method_v1"
EXPECTED={
  "compcancer_9x10":("alizadeh-2000-v1","armstrong-2002-v1",
    "bittner-2000","chen-2002","golub-1999-v1",
    "gordon-2002","pomeroy-2002-v1","shipp-2002-v1","singh-2002"),
  "compcancer_2x2":("west-2001","yeoh-2002-v1"),
  "legacy8":("adenocarcinoma","breast2","cns","colon","dlbcl",
             "leukemia","ovarian","prostate")
}

def need(ok,message):
    if not ok:raise SystemExit("RLMT_V2_19X5_PREFLIGHT_STOP: "+message)

def git(repo,*args):
    return subprocess.run(["git",*args],cwd=repo,check=True,
                          capture_output=True,text=True).stdout.strip()

def rows(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def digest(path):
    hash=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            hash.update(block)
    return hash.hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-ref",required=True)
    ap.add_argument("--new-sha",required=True)
    a=ap.parse_args()
    original=Path(a.original_ref).resolve(strict=True)
    need(ROOT==Path.cwd().resolve(),"run from NEW v2 full-benchmark worktree")
    need(git(ROOT,"rev-parse","HEAD")==a.new_sha,
         "new experimental SHA mismatch")
    need(git(ROOT,"rev-parse","origin/feature/rlmt-v2-19x5-all-methods")
         ==a.new_sha,"experimental origin SHA mismatch")
    need(git(ROOT,"rev-parse","origin/main")==PAPER,
         "published main unexpectedly changed")
    need(git(original,"rev-parse","HEAD")==OLD_SHA,
         "frozen original source SHA changed")
    need(not git(ROOT,"-c","core.filemode=false","status",
                 "--porcelain","--untracked-files=no"),
         "tracked experimental source is dirty")
    need(not git(original,"-c","core.filemode=false","status",
                 "--porcelain","--untracked-files=no"),
         "original tracked source is dirty")
    need(not (ROOT/NEW_DIR).exists(),
         "new results directory exists; no automatic overwrite")
    old=original/OLD_DIR
    need((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "original 285-model run incomplete")
    for name in ("published_5fold_audit_v1/independent_fold_metrics.csv",
                 "published_5fold_accepted_replay_v1/accepted_models.csv"):
        need((old/name).is_file(),"historical audit missing: "+name)
    inputs=rows(old/"input_manifest.csv")
    models=rows(old/"run_manifest.csv")
    need(len(inputs)==19 and len(models)==285,
         "original dataset/model counts differ")
    expected=[(w,d) for w,ds in EXPECTED.items() for d in ds]
    need([(r["world"],r["dataset"]) for r in inputs]==expected,
         "cohort inventory/order mismatch")
    keys={(r["world"],r["dataset"],r["split_id"],r["method"])
          for r in models}
    need(len(keys)==285 and all(r["status"]=="OK" for r in models),
         "original model duplicates or failures")
    count=0
    for x in inputs:
        src=Path(x["source_path"]).resolve(strict=True)
        need(digest(src)==x["sha256"].lower(),
             "input SHA mismatch: "+x["dataset"])
        for i in range(5):
            s=f"rep00_fold{i:02d}"
            for m in ("RLR","RLMT_v1","RLMT_DC_v0"):
                need((x["world"],x["dataset"],s,m) in keys,
                     "missing original paired model "+str((x["dataset"],s,m)))
            count+=1
    need(count==95,"exact sample/fold grid mismatch")
    print("RLMT_V2_19X5_INPUT_SHA256=19/19_PASS")
    print("RLMT_V2_19X5_ORIGINAL_FROZEN_MODELS=285")
    print("RLMT_V2_19X5_ORIGINAL_MATCHED_FOLDS=95")
    print("RLMT_V2_19X5_PLANNED_NEW_FITS=380")
    print("RLMT_V2_19X5_PREFLIGHT=PASS")
    print("PREVIOUSLY_EXPOSED_COHORTS=TRUE")
    print("NO_OLD_RESULT_WRITES=TRUE")

if __name__=="__main__":main()
