#!/usr/bin/env python3
"""Read-only decomposition of verified 19x5 RLMT-DC logloss versus RLR.

No training, parameter tuning, saved-output modifications or report overwrites.
Reports all 19 dataset-level deltas, every accepted DC fold, the top holdout
logloss changes and the individual observed contributions to logloss harm.
"""
from __future__ import annotations
import argparse
import csv
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
REPORT=ROOT/"results/tspdt_recon_19x5_v1/four_method_comparison_v1"
REF_SUBDIR="results/rlmt_dc_published_data_5fold_v1"
DC="RLMT_DC_v0"


def fail(msg):
    raise SystemExit("RLMT_DC_5FOLD_DIAG_STOP: "+msg)


def read(path):
    if not path.is_file():
        fail("missing independently audited file "+str(path))
    with path.open(newline="",encoding="utf-8") as fd:
        return list(csv.DictReader(fd))


def floatval(d,k):
    v=float(d[k])
    if not math.isfinite(v):
        fail("nonfinite "+k)
    return v


def key(d):
    return d["world"],d["dataset"],d["split_id"]


def actual(boolstr):
    v=str(boolstr).strip().lower()
    if v not in ("true","false"):
        fail("not a recorded Boolean split indicator: "+repr(boolstr))
    return v=="true"


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-repo",required=True)
    ap.add_argument("--top-patients",type=int,default=10)
    args=ap.parse_args()
    if args.top_patients<0 or args.top_patients>50:
        fail("top-patients must be between 0 and 50")
    if Path.cwd().resolve()!=ROOT:
        fail("run inside separate TSPDT worktree")
    ref=Path(args.reference_repo).resolve(strict=True)
    refout=ref/REF_SUBDIR
    if not (refout/"status.txt").is_file() or \
       (refout/"status.txt").read_text().strip()!="RUNNER_EXIT=0":
        fail("285-fit RLR/RLMT benchmark not successful")
    dataset=read(REPORT/"matched_19_dataset_means.csv")
    folds=read(REPORT/"matched_95_folds.csv")
    manifest=read(refout/"run_manifest.csv")
    allmetrics=read(REPORT/"all_380_model_fold_metrics.csv")
    if not (len(dataset)==19 and len(folds)==95 and
            len(manifest)==285 and len(allmetrics)==380):
        fail("incomplete four-method audited data")
    mf={}
    for r in manifest:
        k=(*key(r),r["method"])
        if k in mf: fail("duplicate reference fit "+str(k))
        mf[k]=r
    foldmap={}
    for r in folds:
        k=key(r)
        if k in foldmap: fail("duplicate paired fold "+str(k))
        foldmap[k]=r
    drows=[]
    accepted=[]
    for d in dataset:
        w,ds=d["world"],d["dataset"]
        ds_folds=[r for r in folds if r["world"]==w and r["dataset"]==ds]
        if len(ds_folds)!=5 or len({r["split_id"] for r in ds_folds})!=5:
            fail("incomplete five-fold dataset "+ds)
        delta=floatval(d,f"{DC}__delta_mean_logloss_vs_RLR")
        fold_deltas=[floatval(r,f"{DC}__delta_logloss_vs_RLR")
                     for r in ds_folds]
        if abs(delta-statistics.mean(fold_deltas))>1e-9:
            fail("dataset aggregate differs from paired fold evidence: "+ds)
        nacc=0
        for row in ds_folds:
            k=key(row)
            match=mf.get((*k,DC))
            if match is None:
                fail("missing DC fit "+str(k))
            accepted_split=actual(match["accepted_split"])
            nacc+=int(accepted_split)
            fd=floatval(row,f"{DC}__delta_logloss_vs_RLR")
            if not accepted_split and abs(fd)>1e-9:
                fail("unexpected unsplit-tree deviation "+str(k))
            if accepted_split and int(match["realized_leaves"])!=2:
                fail("accepted DC tree does not have two leaves "+str(k))
            if accepted_split:
                accepted.append((
                    w,ds,row["split_id"],fd,
                    floatval(row,f"{DC}__delta_accuracy_vs_RLR"),
                    floatval(row,f"{DC}__delta_brier_vs_RLR"),
                    match
                ))
        drows.append((
            ds,w,delta,
            floatval(d,f"{DC}__delta_mean_accuracy_vs_RLR"),
            floatval(d,f"{DC}__delta_mean_balanced_accuracy_vs_RLR"),
            floatval(d,f"{DC}__delta_mean_brier_vs_RLR"),nacc,
            sum(z>1e-10 for z in fold_deltas),
            sum(z < -1e-10 for z in fold_deltas)
        ))

    mean=statistics.mean(d[2] for d in drows)
    print("=== ALL 19 DATASETS — RLMT-DC MINUS RLR ===")
    print(f"{'dataset':23} {'dLL':>11} {'dAcc':>10} "
          f"{'dBalAcc':>10} {'dBrier':>11} {'accepted':>8}")
    for d in sorted(drows,key=lambda r:-r[2]):
        print(f"{d[0]:23} {d[2]:+11.6f} {d[3]:+10.5f} "
              f"{d[4]:+10.5f} {d[5]:+11.6f} {d[6]:3d}/5")
    print(f"DATASET_MEAN_DC_MINUS_RLR_LOGLOSS={mean:+.9f}")
    print(f"DATASETS_LL_WORSE={sum(z[2]>1e-10 for z in drows)}")
    print(f"DATASETS_LL_BETTER={sum(z[2]< -1e-10 for z in drows)}")
    print(f"DATASETS_LL_TIED={sum(abs(z[2])<=1e-10 for z in drows)}")
    print("=== EVERY ACCEPTED TWO-LEAF DC FOLD ===")
    print(f"{'dataset':23} {'fold':>12} {'dLL':>11} "
          f"{'dAcc':>10} {'dBrier':>11}")
    for w,ds,f,ll,acc,bs,_ in sorted(accepted,key=lambda t:-t[3]):
        print(f"{ds:23} {f:>12} {ll:+11.6f} "
              f"{acc:+10.5f} {bs:+11.6f}")
    print(f"ACCEPTED_DC_FOLDS={len(accepted)}/95")
    print("=== WORST OBSERVED INDIVIDUAL DC LOGLOSS INCREASES ===")
    worst=[]
    for w,ds,f,delta,acc,bs,r_dc in accepted:
        k=(w,ds,f)
        r_baseline=mf.get((*k,"RLR"))
        if r_baseline is None:
            fail("missing base RLR model "+str(k))
        pair=[]
        for m,r in (("RLR",r_baseline),(DC,r_dc)):
            path=(ref/str(r["model_dir"])/"predictions_test.csv").resolve()
            if not path.is_relative_to(refout.resolve()):
                fail("model prediction outside immutable reference output")
            pair.append(read(path))
        baseline,candidate=pair
        if len(baseline)!=len(candidate):
            fail("unpaired heldout predictions "+str(k))
        order=foldmap[k]
        pos=order["class1_reference_ARFF_order"]
        sample_diffs=[]
        for ix,(a,b) in enumerate(zip(baseline,candidate),start=1):
            if (a["sample_id"]!=b["sample_id"] or
                a["y_true"]!=b["y_true"]):
                fail("unpaired sample ID/truth "+str(k))
            y=1 if a["y_true"]==pos else 0
            probs=[float(a["prob"]),float(b["prob"])]
            if not all(0<=p<=1 and math.isfinite(p) for p in probs):
                fail("invalid saved probability "+str(k))
            clipped=[min(1-1e-15,max(1e-15,p)) for p in probs]
            losses=[-(y*math.log(p)+(1-y)*math.log1p(-p))
                    for p in clipped]
            sample_diffs.append(losses[1]-losses[0])
            worst.append((losses[1]-losses[0],ds,f,ix,
                          a["y_true"],probs[0],probs[1]))
        if abs(statistics.mean(sample_diffs)-delta)>1e-9:
            fail("individual logloss deltas don't match fold "+str(k))
    print(f"{'dataset':22} {'fold':>12} {'index':>6} "
          f"{'true':>7} {'p_RLR':>10} {'p_DC':>10} {'dLoss':>11}")
    for loss,ds,f,ix,y,p0,p1 in sorted(worst,reverse=True)[:args.top_patients]:
        print(f"{ds:22} {f:>12} {ix:6d} {y:>7} "
              f"{p0:10.6f} {p1:10.6f} {loss:+11.6f}")
    print("RLMT_DC_READONLY_OUTCOME_DIAG=PASS")


if __name__=="__main__":
    main()
