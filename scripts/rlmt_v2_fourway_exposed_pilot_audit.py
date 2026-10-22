#!/usr/bin/env python3
"""Independent 4-way paper-paired audit (2 datasets × 3 folds × 4 methods).

SOURCE: 12 H1/H2 previous pilot fits + 12 new TSPDT_GINI/RLMT_LA
fits, ALL read-only. Original archived RLR 285-fold prediction files
are read-only too. Never invent missing models, never retrain here.
Scientific stage: previously exposed retrospective development.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
RESULT="results/rlmt_v2_fourway_pilot_2x3_v1"
EARLIER="results/rlmt_v2_exposed_pilot_2x3_v1"
ARCHIVE="results/rlmt_dc_published_data_5fold_v1"
COHORTS=("chen-2002","gordon-2002")
FOLDS=tuple(f"rep00_fold{i:02d}" for i in range(3))
METHODS=("H1","H2","TSPDT_GINI","RLMT_LA")

def need(ok,message):
    if not ok:
        raise SystemExit("FOURWAY_AUDIT_STOP: "+message)

def read(path):
    need(path.is_file(),"missing file "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        records=list(csv.DictReader(f))
    need(bool(records),"empty file "+str(path))
    return records

def numeric(value,context):
    try:
        z=float(value)
    except (ValueError,TypeError) as e:
        raise SystemExit("FOURWAY_AUDIT_STOP: non-numeric "+context) from e
    need(math.isfinite(z),"not finite "+context)
    return z

def metric(y,p,positive):
    need(len(y)==len(p)>0,"empty predictions")
    one=[int(v==positive) for v in y]
    need(set(one)=={0,1},"one-class outer fold not supported")
    need(all(0<=v<=1 for v in p),"probability out of [0,1]")
    q=[min(1-1e-15,max(1e-15,v)) for v in p]
    loss=-sum(a*math.log(b)+(1-a)*math.log1p(-b)
              for a,b in zip(one,q))/len(one)
    brier=sum((a-b)**2 for a,b in zip(one,p))/len(one)
    correct=sum((v>=.5)==bool(a) for v,a in zip(p,one))/len(one)
    return loss,brier,correct

def within(a,b,tol=1e-10):
    return len(a)==len(b) and all(abs(x-y)<=tol for x,y in zip(a,b))

def output_csv(path,records):
    with path.open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--reference-repo",required=True)
    parser.add_argument("--previous-pilot",required=True)
    a=parser.parse_args()
    need(ROOT==Path.cwd().resolve(),"must run from fourway worktree")
    old=Path(a.reference_repo).resolve(strict=True)/ARCHIVE
    previous=Path(a.previous_pilot).resolve(strict=True)/EARLIER
    here=(ROOT/RESULT).resolve(strict=True)
    need((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "archival RLR original incomplete")
    need((previous/"status.txt").read_text().strip()=="RLMT_V2_PILOT_EXIT=0",
         "previous H1/H2 original pilot not complete")
    newrows=read(here/"run_manifest.csv")
    prevrows=read(previous/"run_manifest.csv")
    originals=read(old/"run_manifest.csv")
    need(len(newrows)==12 and len(prevrows)==12 and len(originals)==285,
         "unexpected inventory lengths")
    archived={(r["world"],r["dataset"],r["split_id"],r["method"]):r
              for r in originals}
    need(len(archived)==285,"duplicate historical model IDs")
    all_rows={}
    for mode,rows,base in [
        ("oldH1H2",prevrows,previous.parent.parent),
        ("newTrees",newrows,ROOT)]:
        for record in rows:
            method=record["mode"]
            need(method in (("H1","H2") if mode=="oldH1H2" else
                            ("TSPDT_GINI","RLMT_LA")),
                 "unexpected method name or source")
            key=(record["dataset"],record["split_id"],method)
            need(key not in all_rows,"duplicate paired method")
            need(record["status"]=="OK" and
                 record["scientific_stage"]==
                 "EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
                 "missing result status")
            all_rows[key]=(record,base)
    need(set(all_rows)=={
        (dataset,fold,mode)
        for dataset in COHORTS for fold in FOLDS for mode in METHODS
    },"incomplete exact 2 x 3 x 4 design")

    combined=[]
    for dataset in COHORTS:
        for fold in FOLDS:
            models={}
            baselines={}
            yy=None
            sample_keys=None
            original_dir=None
            original_probs=None
            for mode in METHODS:
                record,base=all_rows[(dataset,fold,mode)]
                name=(base/record["model_dir"]).resolve(strict=True)
                need(name.is_relative_to(base.resolve()),
                     "experimental model path escapes model root")
                observed=read(name/"predictions_test.csv")
                keys=[(r["sample_id"],r["y_true"]) for r in observed]
                need(len(keys)==len(set(keys)),"duplicate outer-test samples")
                if sample_keys is None:
                    sample_keys=keys
                    yy=[r["y_true"] for r in observed]
                else:
                    need(keys==sample_keys,
                         "sample/true class mismatch across four modes")
                world=record["world"]
                oldrecord=archived.get((world,dataset,fold,"RLR"))
                need(oldrecord is not None,"archived original comparator missing")
                archive_model=(Path(a.reference_repo).resolve()/
                               oldrecord["model_dir"]).resolve(strict=True)
                need(archive_model.is_relative_to(old),
                     "archived original outside frozen archive")
                frozen=read(archive_model/"predictions_test.csv")
                need([(r["sample_id"],r["y_true"]) for r in frozen]==keys,
                     "archived RLR sample alignment mismatch")
                if original_dir is None:
                    original_dir=archive_model
                    original_probs=[numeric(r["prob"],"old RLR") for r in frozen]
                else:
                    need(archive_model==original_dir,
                         "different original RLR comparator among methods")
                models[mode]=[numeric(r["prob"],"fourway model")
                              for r in observed]
                baselines[mode]=[numeric(r["parent_prob"],"aligned parent")
                                 for r in observed]
                need(within(baselines["H1"],baselines[mode],tol=1e-12),
                     "aligned parent differs across H1/H2/Gini/LA")
                if mode in ("H1","H2"):
                    accepted=record["accepted"].lower()=="true"
                else:
                    accepted=record["accepted"].lower()=="true"
                need(record["accepted"].lower() in ("true","false"),
                     "unexpected accepted decision encoding")
                if not accepted:
                    need(models[mode]==baselines[mode],
                         "no-split fallback is not bit-exact")
            support=read(original_dir.parent/"class_support.csv")
            need(len(support)==2,"original class support missing")
            positive=support[1]["class"]
            original_m=metric(yy,original_probs,positive)
            root_m=metric(yy,baselines["H1"],positive)
            for mode in METHODS:
                record,_=all_rows[(dataset,fold,mode)]
                m=metric(yy,models[mode],positive)
                for key,val in [
                    ("new_logloss",m[0]),("new_brier",m[1]),
                    ("aligned_parent_logloss",root_m[0]),
                    ("aligned_brier",root_m[1]),
                    ("original_RLR_logloss",original_m[0]),
                    ("original_RLR_brier",original_m[1]),
                    ("delta_vs_aligned",m[0]-root_m[0]),
                    ("delta_vs_original",m[0]-original_m[0])]:
                    need(abs(numeric(record[key],mode+"/"+key)-val)<1e-9,
                         "stored metric/independent prediction mismatch "+mode+"/"+key)
                combined.append(dict(
                    dataset=dataset,world=record["world"],split_id=fold,
                    mode=mode,accepted=record["accepted"].lower()=="true",
                    n_test=len(yy),logloss=m[0],brier=m[1],
                    accuracy=m[2],aligned_logloss=root_m[0],
                    original_RLR_logloss=original_m[0],
                    delta_LL_vs_aligned=m[0]-root_m[0],
                    delta_LL_vs_original=m[0]-original_m[0],
                    scientific_stage="EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY"))
    need(len(combined)==24,"did not produce 24 paired fold/method rows")
    report=here/"independent_fourway_comparison"
    need(not report.exists(),"prior independent audit exists")
    report.mkdir()
    output_csv(report/"24_method_fold_metrics.csv",combined)
    summaries=[]
    for method in METHODS:
        r=[x for x in combined if x["mode"]==method]
        need(len(r)==6,"missing method fold")
        benefits=sum(x["delta_LL_vs_aligned"]< -1e-10 for x in r)
        harms=sum(x["delta_LL_vs_aligned"]> 1e-10 for x in r)
        per_cohort={dataset:statistics.mean(
                x["delta_LL_vs_aligned"] for x in r
                if x["dataset"]==dataset)
              for dataset in COHORTS}
        item=dict(method=method,folds=6,
                  accepted=sum(x["accepted"] for x in r),
                  improvement_folds=benefits,
                  harm_folds=harms,
                  tie_folds=6-benefits-harms,
                  mean_fold_delta_LL_vs_aligned=statistics.mean(
                      x["delta_LL_vs_aligned"] for x in r),
                  dataset_equal_mean_delta_LL_vs_aligned=statistics.mean(
                      per_cohort.values()),
                  mean_fold_delta_LL_vs_original=statistics.mean(
                      x["delta_LL_vs_original"] for x in r))
        summaries.append(item)
        print(f"{method}: accepted={item['accepted']}/6 "
              f"win/loss/tie={benefits}/{harms}/{6-benefits-harms} "
              f"mean_dLL_vs_aligned="
              f"{item['mean_fold_delta_LL_vs_aligned']:+.9f}")
    output_csv(report/"four_method_summary.csv",summaries)
    (report/"audit_meta.json").write_text(json.dumps({
        "stage":"EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
        "fold_count":6,"methods":list(METHODS),
        "n_method_fold_records":24,
        "reused_prior_H1H2_fits":12,
        "new_TSPDT_GINI_RLMT_LA_fits":12,
        "independent_saved_prediction_replay":True,
        "retrained_old_models":False,
        "source_of_published_RLR":"read-only frozen original 285 fits"
    },indent=2)+"\n",encoding="utf-8")
    print("RLMT_V2_FOURWAY_AUDITED_METHOD_FOLDS=24")
    print("RLMT_V2_FOURWAY_PAIRED_SPLITS=6")
    print("RLMT_V2_FOURWAY_INDEPENDENT_AUDIT=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY")
    print("REPORT="+str(report))

if __name__=="__main__":
    main()
