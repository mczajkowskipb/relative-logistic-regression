#!/usr/bin/env python3
"""Independent patient-level replay of 380 matched paper-fold RLR/RLMT fits.

All previously exposed 19 cohorts. All four real experimental methods
must be accounted for. Each of 95 folds must have 4 complete models.
Archived published original RLR remains read-only.
Does not retrain models, choose parameters, drop adverse folds or tune.
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
OUT="results/rlmt_v2_19x5_four_method_v1"
ORIG="results/rlmt_dc_published_data_5fold_v1"
MODES=("H1","H2","TSPDT_GINI","RLMT_LA")
EPS=1e-10

def need(ok,msg):
    if not ok:raise SystemExit("RLMT_V2_19X5_AUDIT_STOP: "+msg)

def records(file):
    need(file.is_file(),"missing "+str(file))
    with file.open(newline="",encoding="utf-8") as f:
        data=list(csv.DictReader(f))
    need(bool(data),"empty "+str(file))
    return data

def numeric(v,ctx):
    try:value=float(v)
    except (ValueError,TypeError) as e:
        raise SystemExit("RLMT_V2_19X5_AUDIT_STOP: numeric "+ctx) from e
    need(math.isfinite(value),"nonfinite "+ctx)
    return value

def metrics(y,p,positive):
    need(len(y)==len(p)>0,"empty predictions")
    yy=[int(t==positive) for t in y]
    need(set(yy)=={0,1},"one-class test labels")
    need(all(0<=v<=1 for v in p),"invalid probability")
    pp=[min(1-1e-15,max(1e-15,v)) for v in p]
    ll=-sum(t*math.log(v)+(1-t)*math.log1p(-v)
            for t,v in zip(yy,pp))/len(pp)
    brier=sum((v-t)**2 for t,v in zip(yy,p))/len(pp)
    acc=sum(int((v>=.5)==bool(t)) for t,v in zip(p,yy))/len(pp)
    n0=sum(t==0 for t in yy)
    n1=len(yy)-n0
    ba=(sum((v<.5) and t==0 for v,t in zip(p,yy))/n0+
        sum((v>=.5) and t==1 for v,t in zip(p,yy))/n1)/2
    return ll,brier,acc,ba

def emit(path,data):
    need(bool(data),"cannot write empty "+str(path))
    with path.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader()
        w.writerows(data)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-ref",required=True)
    a=ap.parse_args()
    need(Path.cwd().resolve()==ROOT,"run from v2 19x5 worktree")
    oldbase=Path(a.original_ref).resolve(strict=True)
    old=(oldbase/ORIG).resolve(strict=True)
    current=(ROOT/OUT).resolve(strict=True)
    need((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "historical original run incomplete")
    original=records(old/"run_manifest.csv")
    new=records(current/"run_manifest.csv")
    input_manifest=records(old/"input_manifest.csv")
    need(len(original)==285 and len(new)==380 and len(input_manifest)==19,
         "archived/incoming number of fit rows differs")
    need(all(r["status"]=="OK" for r in new),
         "at least one model failed; cannot claim 380 complete fits")
    legacy={(r["world"],r["dataset"],r["split_id"],r["method"]):r
            for r in original}
    four={(r["world"],r["dataset"],r["split_id"],r["mode"]):r
          for r in new}
    need(len(legacy)==285 and len(four)==380,
         "duplicate historical or new model-fit key")
    for x in input_manifest:
        world,dataset=x["world"],x["dataset"]
        for i in range(5):
            fold=f"rep00_fold{i:02d}"
            need((world,dataset,fold,"RLR") in legacy,
                 "missing original RLR fold")
            for method in MODES:
                need((world,dataset,fold,method) in four,
                     "missing 19x5 algorithm fold "+str((dataset,fold,method)))
    paired=[]
    for x in input_manifest:
        world,dataset=x["world"],x["dataset"]
        for i in range(5):
            sid=f"rep00_fold{i:02d}"
            original_item=legacy[(world,dataset,sid,"RLR")]
            rdir=(oldbase/original_item["model_dir"]).resolve(strict=True)
            need(rdir.is_relative_to(old),"original RLR model escapes archived root")
            frozen=records(rdir/"predictions_test.csv")
            support=records(rdir.parent/"class_support.csv")
            need(len(support)==2,"original class-support labels missing")
            positive=support[1]["class"]
            k=[(r["sample_id"],r["y_true"]) for r in frozen]
            need(len(k)==len(set(k)),"original duplicated patient/sample")
            y=[r["y_true"] for r in frozen]
            origp=[numeric(r["prob"],"original RLR probability") for r in frozen]
            origscore=metrics(y,origp,positive)
            aligned=None
            aligned_score=None
            for m in MODES:
                row=four[(world,dataset,sid,m)]
                path=(ROOT/row["model_dir"]).resolve(strict=True)
                need(path.is_relative_to(current),
                     "experimental model escapes root")
                observed=records(path/"predictions_test.csv")
                need([(r["sample_id"],r["y_true"]) for r in observed]==k,
                     "cross-model per-sample alignment failure")
                p=[numeric(r["prob"],"model prob") for r in observed]
                pp=[numeric(r["parent_prob"],"aligned parent") for r in observed]
                if aligned is None:
                    aligned=pp
                    aligned_score=metrics(y,pp,positive)
                else:
                    need(all(abs(a-b)<1e-12 for a,b in zip(aligned,pp)),
                         "different aligned parent within same fold")
                need(row["accepted"].lower() in ("true","false"),
                     "invalid accepted status")
                accepted=row["accepted"].lower()=="true"
                if not accepted:
                    need(p==pp,"rejected model not exact aligned parent")
                fresh=metrics(y,p,positive)
                for field,v in [
                  ("new_logloss",fresh[0]),("new_brier",fresh[1]),
                  ("aligned_parent_logloss",aligned_score[0]),
                  ("aligned_brier",aligned_score[1]),
                  ("original_RLR_logloss",origscore[0]),
                  ("original_RLR_brier",origscore[1]),
                  ("delta_vs_aligned",fresh[0]-aligned_score[0]),
                  ("delta_vs_original",fresh[0]-origscore[0])]:
                    need(abs(numeric(row[field],field)-v)<1e-9,
                         "metric mismatch for "+str((dataset,sid,m,field)))
                paired.append(dict(
                  world=world,dataset=dataset,split_id=sid,method=m,
                  n_test=len(y),accepted=accepted,
                  logloss=fresh[0],brier=fresh[1],
                  accuracy=fresh[2],balanced_accuracy=fresh[3],
                  aligned_logloss=aligned_score[0],
                  original_RLR_logloss=origscore[0],
                  delta_LL_vs_aligned=fresh[0]-aligned_score[0],
                  delta_LL_vs_original=fresh[0]-origscore[0],
                  delta_Brier_vs_aligned=fresh[1]-aligned_score[1],
                  warning_total=int(row["warning_total"]),
                  low_class_warning_total=int(row["warning_low_class"])))
    need(len(paired)==380,"unexpected number of verified predictions")
    report=current/"independent_19x5_audit"
    need(not report.exists(),"prior independent report root exists")
    report.mkdir()
    emit(report/"380_fold_method_metrics.csv",paired)
    summary=[]
    per_dataset=[]
    for method in MODES:
        entries=[r for r in paired if r["method"]==method]
        need(len(entries)==95,"incomplete model 95-fold inventory "+method)
        for x in input_manifest:
            world,dataset=x["world"],x["dataset"]
            z=[r for r in entries if r["world"]==world and
                r["dataset"]==dataset]
            need(len(z)==5,"not exactly 5 paired folds per dataset")
            mean_ll=statistics.mean(r["delta_LL_vs_aligned"] for r in z)
            mean_orig=statistics.mean(r["delta_LL_vs_original"] for r in z)
            per_dataset.append(dict(method=method,
                world=world,dataset=dataset,
                accepted_folds=sum(r["accepted"] for r in z),
                mean_delta_ll_vs_aligned=mean_ll,
                mean_delta_ll_vs_original=mean_orig,
                mean_delta_brier_vs_aligned=statistics.mean(
                  r["delta_Brier_vs_aligned"] for r in z),
                warning_count=sum(r["warning_total"] for r in z),
                low_class_warnings=sum(r["low_class_warning_total"] for r in z)))
        ds=[r for r in per_dataset if r["method"]==method]
        wins=sum(r["mean_delta_ll_vs_aligned"]< -EPS for r in ds)
        harms=sum(r["mean_delta_ll_vs_aligned"]>EPS for r in ds)
        summary.append(dict(method=method,folds=95,datasets=19,
            accepted_folds=sum(r["accepted"] for r in entries),
            dataset_improved=wins,dataset_worsened=harms,
            dataset_tied=19-wins-harms,
            mean_equal_dataset_delta_ll_vs_aligned=statistics.mean(
                r["mean_delta_ll_vs_aligned"] for r in ds),
            mean_equal_dataset_delta_ll_vs_original=statistics.mean(
                r["mean_delta_ll_vs_original"] for r in ds),
            total_warnings=sum(r["warning_total"] for r in entries),
            total_low_class_warnings=sum(
                r["low_class_warning_total"] for r in entries)))
    emit(report/"method_19_dataset_means.csv",per_dataset)
    emit(report/"four_method_summary.csv",summary)
    (report/"audit_metadata.json").write_text(json.dumps({
      "study":"19 previously exposed biomedical datasets x5 old outer folds",
      "methods":list(MODES),"new_fits":380,"original_models_retrained":0,
      "patient_test_prediction_alignment":True,
      "frozen_original_RLR_predictions_reused_readonly":True,
      "scientific_stage":"EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
      "warning_reports_captured":True},indent=2)+"\n",encoding="utf-8")
    for r in summary:
        print(f"{r['method']}: accepted={r['accepted_folds']}/95 "
              f"datasetLL win/loss/tie="
              f"{r['dataset_improved']}/{r['dataset_worsened']}/{r['dataset_tied']} "
              f"mean_dLL_vs_ALIGNED={r['mean_equal_dataset_delta_ll_vs_aligned']:+.9f} "
              f"mean_dLL_vs_ORIGINAL_RLR="
              f"{r['mean_equal_dataset_delta_ll_vs_original']:+.9f} "
              f"warnings={r['total_warnings']} "
              f"low_class_warnings={r['total_low_class_warnings']}")
    print("RLMT_V2_19X5_PATIENT_LEVEL_PREDICTION_AUDIT=PASS")
    print("RLMT_V2_19X5_ALL_380_NEW_MODELS_AUDITED=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY")
    print("REPORT="+str(report))

if __name__=="__main__":main()
