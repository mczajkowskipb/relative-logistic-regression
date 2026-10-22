#!/usr/bin/env python3
"""Independent read-only patient-level audit, fixed 19×5×5 RLMT-v4.

Audits saved historical-original-RLR parity AND computes original v4
model logloss/Brier on precisely the same 95 heldout fold patient IDs.
Equal-dataset mean, leave-one-dataset-out and conditional rank repair
effect vs global calibration; never omit losses.
"""
from __future__ import annotations
import argparse
import csv
import math
import statistics as stats
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/rlmt_v4_rank_repair_19x5_v1"
OLD="results/rlmt_dc_published_data_5fold_v1"
METHODS=("ORIGINAL_RLR","GLOBAL_CAL",
         "CONF_ROUTE","ACTIVITY_ROUTE","SOFT_GAP_ROUTE")
TOL=1e-10

def require(condition,message):
    if not condition:
        raise SystemExit("RLMT4_19X5_AUDIT_STOP: "+str(message))

def csvread(path):
    require(path.is_file(),"missing "+str(path))
    with path.open("r",newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    require(bool(rows),"empty CSV "+str(path))
    return rows

def number(x):
    try:
        v=float(x)
    except (ValueError,TypeError) as e:
        raise SystemExit("RLMT4_19X5_AUDIT_STOP: nonnumeric") from e
    require(math.isfinite(v),"nonfinite value")
    return v

def proper(y,prob,positive):
    require(len(y)==len(prob)>0,"inconsistent predictions")
    require(all(0<=p<=1 for p in prob),"probability outside [0,1]")
    binary=[int(t==positive) for t in y]
    require(set(binary)=={0,1},"heldout fold missing class")
    q=[max(1e-15,min(1-1e-15,p)) for p in prob]
    ll=-sum(t*math.log(p)+(1-t)*math.log1p(-p)
             for t,p in zip(binary,q))/len(y)
    brier=sum((t-p)**2 for t,p in zip(binary,prob))/len(y)
    return ll,brier

def classification(y,prob,positive):
    """Binary balanced accuracy and macro-F1 at fixed 0.5 threshold."""
    t=[int(z==positive) for z in y]
    h=[int(p>=0.5) for p in prob]
    recalls=[]
    f1=[]
    for label in (0,1):
        tp=sum(x==label and z==label for x,z in zip(h,t))
        fp=sum(x==label and z!=label for x,z in zip(h,t))
        fn=sum(x!=label and z==label for x,z in zip(h,t))
        require(tp+fn>0,"missing class at fixed threshold evaluation")
        recalls.append(tp/(tp+fn))
        f1.append(2*tp/(2*tp+fp+fn) if (2*tp+fp+fn)>0 else 0)
    return stats.mean(recalls),stats.mean(f1)

def create_csv(path,rows):
    require(len(rows)>0,"cannot save empty audit")
    with path.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-ref",required=True)
    args=ap.parse_args()
    ref=Path(args.original_ref).resolve(strict=True)
    old=(ref/OLD).resolve(strict=True)
    out=OUT.resolve(strict=True)
    require(Path.cwd().resolve()==ROOT,
            "workdir must equal v4 19x5 source root")
    require((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
            "original archived data incomplete")
    manifest=csvread(out/"run_manifest.csv")
    old_rows=csvread(old/"run_manifest.csv")
    inputs=csvread(old/"input_manifest.csv")
    require(len(manifest)==475 and len(old_rows)==285 and len(inputs)==19,
            "expected 475 v4, 285 archived models, 19 cohorts")
    lookup={(r["world"],r["dataset"],r["split_id"],r["method"]):r
            for r in manifest}
    legacy={(r["world"],r["dataset"],r["split_id"],r["method"]):r
            for r in old_rows}
    need={(r["world"],r["dataset"],f"rep00_fold{i:02d}",method)
          for r in inputs for i in range(5) for method in METHODS}
    require(len(lookup)==475 and set(lookup)==need,
            "wrong, duplicate or missing method-fold keys")
    require(len(legacy)==285,"duplicate archived RLR model keys")

    fold_audit=[]
    dataset_audit=[]
    audited=0
    for item in inputs:
        world=item["world"]
        dataset=item["dataset"]
        per_method={name:[] for name in METHODS}
        for fold_index in range(5):
            fold=f"rep00_fold{fold_index:02d}"
            archived=legacy.get((world,dataset,fold,"RLR"))
            require(archived is not None,"missing original RLR fold")
            frozen=(ref/archived["model_dir"]).resolve(strict=True)
            require(frozen.is_relative_to(old),
                    "archived model path outside frozen reference")
            frozenp=csvread(frozen/"predictions_test.csv")
            classes=csvread(frozen.parent/"class_support.csv")
            require(len(classes)==2,"original positive class unavailable")
            positive=classes[1]["class"]
            frozen_ids=[(r["sample_id"],r["y_true"]) for r in frozenp]
            require(len(frozen_ids)==len(set(frozen_ids)),
                    "duplicate original patient IDs")
            y=[r["y_true"] for r in frozenp]
            pold=[number(r["prob"]) for r in frozenp]
            original_ll,original_brier=proper(y,pold,positive)
            original_ba,original_f1=classification(y,pold,positive)
            probs={}
            rowfolds=[]
            for method in METHODS:
                key=(world,dataset,fold,method)
                record=lookup[key]
                require(record["status"]=="OK" and
                        record["scientific_stage"]==
                        "EXPOSED_RETROSPECTIVE_DEV_ONLY",
                        "result not OK or stage mislabeled")
                path=(ROOT/record["model_dir"]).resolve(strict=True)
                require(path.is_relative_to(out),
                        "new predictions escape v4 result root")
                predictions=csvread(path/"predictions_test.csv")
                require([(r["sample_id"],r["y_true"]) for r in predictions]
                        ==frozen_ids,"heldout patient grid/label changed")
                p=[number(r["prob"]) for r in predictions]
                parent=[number(r["original_RLR_prob"]) for r in predictions]
                require(max(abs(a-b) for a,b in zip(parent,pold))<1e-9,
                        "not exactly original archived RLR parent")
                if method=="ORIGINAL_RLR":
                    require(max(abs(a-b) for a,b in zip(p,pold))<1e-9,
                            "ORIGINAL_RLR v4 differs from archived paper")
                ll,brier=proper(y,p,positive)
                ba,f1=classification(y,p,positive)
                for name,val in [
                    ("logloss",ll),
                    ("original_logloss",original_ll),
                    ("delta_vs_original",ll-original_ll),
                    ("brier",brier)]:
                    require(abs(number(record[name])-val)<1e-9,
                            "saved "+name+" not independently reproduced")
                require(int(record["n_test"])==len(p),
                        "heldout size mismatch")
                split=int(record["route_n"])
                require(split in (0,1),"not a stump model")
                if method in ("ORIGINAL_RLR","GLOBAL_CAL"):
                    require(split==0,"non-DT control has route")
                if split:
                    require(record["gene_a"] and record["gene_b"]
                            and record["gene_a"]!=record["gene_b"],
                            "split genes missing/invalid")
                require(abs(number(record["parity_max_abs"]))<=1e-9,
                        "paper parent parity error unacceptable")
                score=dict(
                    world=world,dataset=dataset,split_id=fold,
                    method=method,n_test=len(p),
                    route_n=split,logloss=ll,brier=brier,
                    balanced_accuracy=ba,macro_f1=f1,
                    original_logloss=original_ll,
                    original_brier=original_brier,
                    original_balanced_accuracy=original_ba,
                    original_macro_f1=original_f1,
                    delta_balanced_accuracy_vs_original=ba-original_ba,
                    delta_macro_f1_vs_original=f1-original_f1,
                    delta_ll_vs_original=ll-original_ll,
                    delta_ll_vs_global=0.0,
                    delta_brier_vs_original=brier-original_brier,
                    delta_brier_vs_global=0.0,
                    gamma=number(record["gamma"]),
                    warning_count=int(record["warning_count"]))
                rowfolds.append(score)
                probs[method]=(ll,brier)
                audited+=1
            global_ll,global_brier=probs["GLOBAL_CAL"]
            for row in rowfolds:
                row["delta_ll_vs_global"]=row["logloss"]-global_ll
                row["delta_brier_vs_global"]=row["brier"]-global_brier
                fold_audit.append(row)
                per_method[row["method"]].append(row)
        for method in METHODS:
            z=per_method[method]
            require(len(z)==5,"not five comparable folds per method")
            dataset_audit.append(dict(
                world=world,dataset=dataset,method=method,
                actual_route_folds=sum(r["route_n"] for r in z),
                mean_delta_ll_vs_original=stats.mean(
                    r["delta_ll_vs_original"] for r in z),
                mean_delta_ll_vs_global=stats.mean(
                    r["delta_ll_vs_global"] for r in z),
                mean_delta_brier_vs_original=stats.mean(
                    r["delta_brier_vs_original"] for r in z),
                mean_delta_brier_vs_global=stats.mean(
                    r["delta_brier_vs_global"] for r in z),
                mean_delta_ba_vs_original=stats.mean(
                    r["delta_balanced_accuracy_vs_original"] for r in z),
                mean_delta_macro_f1_vs_original=stats.mean(
                    r["delta_macro_f1_vs_original"] for r in z),
                mean_abs_gamma=stats.mean(abs(r["gamma"]) for r in z),
                shared_fold_warnings=sum(r["warning_count"] for r in z)))
    require(audited==475 and len(fold_audit)==475 and
            len(dataset_audit)==95,
            "not all 475 method-folds / 95 dataset methods audited")
    summary=[]
    loo=[]
    for method in METHODS:
        z=[r for r in dataset_audit if r["method"]==method]
        require(len(z)==19,"expected exactly 19 source cohorts")
        delta_original=[r["mean_delta_ll_vs_original"] for r in z]
        delta_cal=[r["mean_delta_ll_vs_global"] for r in z]
        win_orig=sum(v < -TOL for v in delta_original)
        loss_orig=sum(v > TOL for v in delta_original)
        win_cal=sum(v < -TOL for v in delta_cal)
        loss_cal=sum(v > TOL for v in delta_cal)
        avg_original=stats.mean(delta_original)
        avg_cal=stats.mean(delta_cal)
        summary.append(dict(
            method=method,n_datasets=19,n_folds=95,
            actual_split_folds=sum(r["actual_route_folds"] for r in z),
            datasets_win_vs_original=win_orig,
            datasets_loss_vs_original=loss_orig,
            datasets_tie_vs_original=19-win_orig-loss_orig,
            datasets_win_vs_global=win_cal,
            datasets_loss_vs_global=loss_cal,
            datasets_tie_vs_global=19-win_cal-loss_cal,
            mean_deltaLL_vs_original=avg_original,
            mean_deltaLL_vs_global=avg_cal,
            mean_deltaBrier_vs_original=stats.mean(
                r["mean_delta_brier_vs_original"] for r in z),
            mean_deltaBrier_vs_global=stats.mean(
                r["mean_delta_brier_vs_global"] for r in z),
            mean_deltaBA_vs_original=stats.mean(
                r["mean_delta_ba_vs_original"] for r in z),
            mean_deltaMacroF1_vs_original=stats.mean(
                r["mean_delta_macro_f1_vs_original"] for r in z)))
        for i,r in enumerate(z):
            loo.append(dict(method=method,omitted_dataset=r["dataset"],
                mean_deltaLL_vs_original_excluded=
                stats.mean(v for j,v in enumerate(delta_original) if j!=i),
                mean_deltaLL_vs_global_excluded=
                stats.mean(v for j,v in enumerate(delta_cal) if j!=i)))
    audit=out/"independent_v4_19x5_audit"
    require(not audit.exists(),"audit already exists, refuse overwrite")
    audit.mkdir()
    create_csv(audit/"475_method_fold_metrics.csv",fold_audit)
    create_csv(audit/"19_dataset_method_means.csv",dataset_audit)
    create_csv(audit/"5_method_summary.csv",summary)
    create_csv(audit/"leave_one_dataset_out.csv",loo)
    for r in summary:
        print(
          f"{r['method']}: routes={r['actual_split_folds']}/95 "
          f"W/L/T_vs_ORIGINAL={r['datasets_win_vs_original']}/"
          f"{r['datasets_loss_vs_original']}/"
          f"{r['datasets_tie_vs_original']} "
          f"W/L/T_vs_GLOBAL_CAL={r['datasets_win_vs_global']}/"
          f"{r['datasets_loss_vs_global']}/"
          f"{r['datasets_tie_vs_global']} "
          f"mean_dLL_vs_ORIGINAL={r['mean_deltaLL_vs_original']:+.9f} "
          f"mean_dLL_vs_GLOBAL={r['mean_deltaLL_vs_global']:+.9f}")
    print("RLMT4_19X5_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS")
    print("RLMT4_19X5_475_NEW_PREDICTION_VECTORS_AUDITED=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
    print("REPORT="+str(audit))
if __name__=="__main__":
    main()
