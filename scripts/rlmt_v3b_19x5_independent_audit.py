#!/usr/bin/env python3
"""Independent 19×5 patient prediction replay for 4 calibration controls.

Verifies immutable archived original RLR per-patient predictions,
4 exact samples/method/fold, identical aligned root, LL/Brier scores,
explicit selected relational model counts, and SHARED warning inventory.
Writes dataset-equal averages, paired route-vs-global-control differences.
No models refit; no method/fold exclusion.
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
RESULT="results/rlmt_v3b_calibration_19x5_v1"
ORIG="results/rlmt_dc_published_data_5fold_v1"
METHODS=("GLOBAL_RIDGE_CAL","STUMP_TSP_RIDGE",
         "STUMP_GINI_RIDGE","STUMP_RESID_RIDGE")
EPS=1e-10

def need(x,msg):
    if not x:raise SystemExit("V3B_AUDIT_STOP: "+msg)

def records(p):
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    need(bool(rows),"empty "+str(p))
    return rows

def value(x):
    try:v=float(x)
    except (TypeError,ValueError) as e:
        raise SystemExit("V3B_AUDIT_STOP: nonnumerical") from e
    need(math.isfinite(v),"nonfinite model output")
    return v

def metrics(y,p,positive):
    need(len(y)==len(p)>0 and all(0<=x<=1 for x in p),
         "probability range error")
    y01=[int(z==positive) for z in y]
    need(set(y01)=={0,1},"missing heldout class")
    q=[max(1e-15,min(1-1e-15,x)) for x in p]
    ll=-sum(t*math.log(z)+(1-t)*math.log1p(-z)
        for t,z in zip(y01,q))/len(y01)
    brier=sum((x-t)**2 for x,t in zip(p,y01))/len(y01)
    return ll,brier

def write(p,rows):
    need(bool(rows),"cannot write empty audit")
    with p.open("x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-ref",required=True)
    a=ap.parse_args()
    ref=Path(a.original_ref).resolve(strict=True)
    source=(ref/ORIG).resolve(strict=True)
    out=(ROOT/RESULT).resolve(strict=True)
    need(Path.cwd().resolve()==ROOT,"run from new v3B worktree")
    need((source/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "original 285-model reference not complete")
    fit=records(out/"run_manifest.csv")
    archived=records(source/"run_manifest.csv")
    dataset_manifest=records(source/"input_manifest.csv")
    need(len(fit)==380 and len(archived)==285 and
         len(dataset_manifest)==19,
         "incomplete source or 380 new results")
    keys={(r["world"],r["dataset"],r["split_id"],r["method"]):r
          for r in fit}
    originals={(r["world"],r["dataset"],r["split_id"],r["method"]):r
               for r in archived}
    need(len(keys)==380 and len(originals)==285,
         "duplicate experiment identifiers")
    needs={
      (r["world"],r["dataset"],f"rep00_fold{j:02d}",method)
      for r in dataset_manifest for j in range(5)
      for method in METHODS
    }
    need(set(keys)==needs,"not exact 19x5x4 predeclared design")
    warning_path=out/"shared_training_warning_inventory.csv"
    counts=defaultdict(int)
    low=defaultdict(int)
    if warning_path.exists():
        for x in records(warning_path):
            k=(x["world"],x["dataset"],x["split_id"])
            n=int(x["frequency"])
            need(n>0,"invalid warning frequency")
            counts[k]+=n
            if "fewer than 8" in x["message"]:
                low[k]+=n
    verified=[]
    fold_comparisons=[]
    for input_row in dataset_manifest:
      world=input_row["world"]
      dataset=input_row["dataset"]
      for fold_idx in range(5):
        fold=f"rep00_fold{fold_idx:02d}"
        old=originals.get((world,dataset,fold,"RLR"))
        need(old is not None,"missing old RLR comparator")
        folder=(ref/old["model_dir"]).resolve(strict=True)
        need(folder.is_relative_to(source),
             "original model outside frozen archive")
        rp=records(folder/"predictions_test.csv")
        support=records(folder.parent/"class_support.csv")
        need(len(support)==2,"historical positive label missing")
        positive=support[1]["class"]
        old_keys=[(x["sample_id"],x["y_true"]) for x in rp]
        need(len(old_keys)==len(set(old_keys)),
             "original duplicated test sample")
        truth=[x["y_true"] for x in rp]
        legacy_ll,legacy_brier=metrics(
            truth,[value(x["prob"]) for x in rp],positive)
        parent=None
        parentll=None
        results={}
        for method in METHODS:
            key=(world,dataset,fold,method)
            m=keys[key]
            need(m["status"]=="OK" and
                 m["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEV_ONLY",
                 "method fit status or stage differs")
            modeldir=(ROOT/m["model_dir"]).resolve(strict=True)
            need(modeldir.is_relative_to(out),
                 "new model output escapes experimental root")
            pred=records(modeldir/"predictions_test.csv")
            need([(x["sample_id"],x["y_true"]) for x in pred]==old_keys,
                 "outer test patient/label mismatch")
            probability=[value(x["prob"]) for x in pred]
            root=[value(x["parent_prob"]) for x in pred]
            if parent is None:
                parent=root
                parentll,_=metrics(truth,root,positive)
            else:
                need(len(parent)==len(root) and all(
                     abs(x-y)<1e-12 for x,y in zip(parent,root)),
                     "global aligned parent changed between method fits")
            ll,brier=metrics(truth,probability,positive)
            for field,calculated in [
                ("scored_logloss",ll),
                ("parent_logloss",parentll),
                ("original_RLR_logloss",legacy_ll),
                ("delta_vs_parent",ll-parentll),
                ("delta_vs_original",ll-legacy_ll),
                ("brier",brier)]:
                need(abs(value(m[field])-calculated)<1e-9,
                     "saved metric inconsistent "+str((key,field)))
            n=int(m["internal_nodes"])
            need(n==0 if method=="GLOBAL_RIDGE_CAL" else n in (0,1),
                 "control incorrectly has split or tree node count invalid")
            if n:
                need(m["gene_a"] and m["gene_b"] and
                     m["gene_a"]!=m["gene_b"],
                     "selected split missing genes")
                need(all(x["leaf_id"] in ("2","3") for x in pred),
                     "route leaf_id invalid")
            else:
                need(all(x["leaf_id"]=="1" for x in pred),
                     "non-route produced multiple leaves")
            need(int(m["n_test"])==len(pred),
                 "test count inconsistent")
            need(int(m["shared_parent_warnings"])==
                 counts[(world,dataset,fold)] and
                 int(m["shared_parent_low_class_warnings"])==
                 low[(world,dataset,fold)],
                 "shared parent warning count inconsistent")
            results[method]=(ll,brier)
            verified.append(dict(
                world=world,dataset=dataset,split_id=fold,
                method=method,n_test=len(pred),
                logloss=ll,brier=brier,
                parent_logloss=parentll,
                original_RLR_logloss=legacy_ll,
                internal_nodes=n,
                route_effect=value(m["route_effect"]),
                delta_ll_vs_parent=ll-parentll,
                delta_ll_vs_original=ll-legacy_ll,
                delta_ll_vs_global_cal=0.0,
                shared_warnings=counts[(world,dataset,fold)]))
        global_ll=results["GLOBAL_RIDGE_CAL"][0]
        for row in verified[-4:]:
            row["delta_ll_vs_global_cal"]=row["logloss"]-global_ll
        fold_comparisons.append(
            (world,dataset,fold,results,legacy_ll,parentll))
    need(len(verified)==380 and len(fold_comparisons)==95,
         "unexpected audited fit count")
    dataset_report=[]
    summary=[]
    for method in METHODS:
        for entry in dataset_manifest:
            world,dataset=entry["world"],entry["dataset"]
            rr=[x for x in verified if x["world"]==world and
                x["dataset"]==dataset and x["method"]==method]
            need(len(rr)==5,"expected five matched folds per dataset")
            dataset_report.append(dict(
               method=method,world=world,dataset=dataset,
               actual_split_folds=sum(x["internal_nodes"]>0 for x in rr),
               mean_delta_vs_parent=statistics.mean(
                   x["delta_ll_vs_parent"] for x in rr),
               mean_delta_vs_original=statistics.mean(
                   x["delta_ll_vs_original"] for x in rr),
               mean_delta_vs_global_cal=statistics.mean(
                   x["delta_ll_vs_global_cal"] for x in rr),
               mean_abs_route_effect=statistics.mean(
                   abs(x["route_effect"]) for x in rr)))
        dd=[x for x in dataset_report if x["method"]==method]
        need(len(dd)==19,"dataset count must equal 19")
        improved=sum(x["mean_delta_vs_global_cal"]< -EPS for x in dd)
        worsened=sum(x["mean_delta_vs_global_cal"]>EPS for x in dd)
        summary.append(dict(
            method=method,datasets=19,folds=95,
            actual_split_folds=sum(x["actual_split_folds"] for x in dd),
            datasets_win_vs_global=improved,
            datasets_loss_vs_global=worsened,
            datasets_tie_vs_global=19-improved-worsened,
            mean_delta_vs_aligned=statistics.mean(
                x["mean_delta_vs_parent"] for x in dd),
            mean_delta_vs_original_RLR=statistics.mean(
                x["mean_delta_vs_original"] for x in dd),
            mean_delta_vs_global_cal=statistics.mean(
                x["mean_delta_vs_global_cal"] for x in dd)))
    report=out/"independent_v3b_audit"
    need(not report.exists(),"old audit root must not be overwritten")
    report.mkdir()
    write(report/"380_fold_method_metrics.csv",verified)
    write(report/"19_dataset_method_metrics.csv",dataset_report)
    write(report/"four_method_summary.csv",summary)
    (report/"design.json").write_text(json.dumps({
      "study":"95 matched folds, 19 previously exposed datasets",
      "methods":METHODS,
      "original_RLR_results":"immutable read only",
      "one_shared_parent_OOF_per_fold":True,
      "global_calibration_control":True,
      "same_outer_test_samples_all_4_models":True,
      "data_independent_validation":False,
      "warning_counts_shared_per_fold_not_per_method":True},
      indent=2)+"\n")
    for r in summary:
        print(
          f"{r['method']}: split={r['actual_split_folds']}/95 "
          f"dataset vs calibrated control W/L/T="
          f"{r['datasets_win_vs_global']}/"
          f"{r['datasets_loss_vs_global']}/"
          f"{r['datasets_tie_vs_global']} "
          f"deltaLL_vs_ORIGINAL_RLR="
          f"{r['mean_delta_vs_original_RLR']:+.9f} "
          f"deltaLL_vs_GLOBAL_CAL="
          f"{r['mean_delta_vs_global_cal']:+.9f}")
    print("RLMT_V3B_19X5_INDEPENDENT_PATIENT_AUDIT=PASS")
    print("RLMT_V3B_19X5_AUDITED_FOLD_METHODS=380")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
    print("REPORT="+str(report))

if __name__=="__main__":main()
