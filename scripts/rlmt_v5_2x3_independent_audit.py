#!/usr/bin/env python3
"""Strict patient-level V5 gate and original RLR parity replay, 2×3."""
from __future__ import annotations
import argparse
import csv
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/rlmt_v5_honest_router_2x3_v1"
OLD="results/rlmt_dc_published_data_5fold_v1"
DS=("chen-2002","gordon-2002")
FOLDS=("rep00_fold00","rep00_fold01","rep00_fold02")
MODES=("ORIGINAL_RLR","GLOBAL_CAL","CONF_ROUTE","ACTIVITY_ROUTE")
def ensure(cond,msg):
    if not cond: raise SystemExit("RLMT5_AUDIT_STOP: "+msg)
def read(path):
    ensure(path.is_file(),"not found "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        records=list(csv.DictReader(f))
    ensure(bool(records),"CSV empty "+str(path))
    return records
def number(s):
    try:r=float(s)
    except (ValueError,TypeError) as e:
        raise SystemExit("RLMT5_AUDIT_STOP: numeric parsing") from e
    ensure(math.isfinite(r),"nonfinite value")
    return r
def scores(y,prob,positive):
    ensure(len(y)==len(prob)>0 and all(0<=p<=1 for p in prob),
           "invalid probabilities or patient alignment")
    yy=[int(t==positive) for t in y]
    ensure(set(yy)=={0,1},"both heldout classes required")
    pp=[max(1e-15,min(1-1e-15,p)) for p in prob]
    ll=-sum(t*math.log(p)+(1-t)*math.log1p(-p)
            for t,p in zip(yy,pp))/len(y)
    brier=sum((p-t)**2 for p,t in zip(prob,yy))/len(y)
    return ll,brier

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-ref",required=True)
    args=ap.parse_args()
    ref=Path(args.original_ref).resolve(strict=True)
    archive=(ref/OLD).resolve(strict=True)
    out=OUT.resolve(strict=True)
    ensure(Path.cwd().resolve()==ROOT,"wrong v5 worktree")
    ensure((archive/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
           "original RLR archive not finished")
    results=read(out/"run_manifest.csv")
    original=read(archive/"run_manifest.csv")
    ensure(len(results)==6 and len(original)==285,"wrong inventory")
    bykey={(r["dataset"],r["split_id"]):r for r in results}
    byorig={(r["dataset"],r["split_id"],r["method"]):r for r in original}
    ensure(len(bykey)==6 and
           set(bykey)=={(d,s) for d in DS for s in FOLDS},
           "pilot fold key missing/duplicated")
    records=[]
    for d in DS:
      for f in FOLDS:
        row=bykey[(d,f)]
        refrow=byorig[(d,f,"RLR")]
        folder=(ref/refrow["model_dir"]).resolve(strict=True)
        ensure(folder.is_relative_to(archive),"old model dir escaped")
        old=read(folder/"predictions_test.csv")
        cls=read(folder.parent/"class_support.csv")
        ensure(len(cls)==2,"old positive class list changed")
        positive=cls[1]["class"]
        target=[(p["sample_id"],p["y_true"]) for p in old]
        ensure(len(target)==len(set(target)),"duplicate original patients")
        originalprob=[number(p["prob"]) for p in old]
        root_ll,root_brier=scores(
            [p["y_true"] for p in old],originalprob,positive)
        loc=(ROOT/row["model_dir"]).resolve(strict=True)
        ensure(loc.is_relative_to(out),"new model dir escaped")
        p=read(loc/"predictions_test.csv")
        ensure([(x["sample_id"],x["y_true"]) for x in p]==target,
               "outer heldout sample IDs or true labels changed")
        newprob=[number(x["prob"]) for x in p]
        parent=[number(x["original_RLR_prob"]) for x in p]
        ensure(max(abs(a-b) for a,b in zip(parent,originalprob))<1e-9,
               "original RLR probabilities not identical")
        ensure(number(row["max_parent_parity"])<1e-9,
               "reported parent parity exceeds 1e-9")
        ll,brier=scores([x["y_true"] for x in p],newprob,positive)
        ensure(abs(ll-number(row["logloss"]))<1e-9 and
               abs(root_ll-number(row["original_logloss"]))<1e-9 and
               abs((ll-root_ll)-number(row["delta_vs_original"]))<1e-9 and
               abs(brier-number(row["brier"]))<1e-9 and
               abs(root_brier-number(row["original_brier"]))<1e-9,
               "saved loss inconsistent with independent heldout replay")
        method=row["selected_method"]
        w=number(row["selected_weight"])
        ensure(method in MODES and 0<=w<=1 and
               (method=="ORIGINAL_RLR" if w==0 else
                method!="ORIGINAL_RLR"),
               "invalid train-only chosen method or weight")
        ensure(all(p0["selected_method"]==method and
                   abs(number(p0["selected_weight"])-w)<1e-14
                   for p0 in p),
               "saved patient-level router inconsistent")
        ensure((row["accepted"]=="TRUE")== (w>0) and
               int(row["route_n"])==int(w>0 and
                   method in ("CONF_ROUTE","ACTIVITY_ROUTE")),
               "selected method vs route count/acceptance inconsistent")
        if row["route_n"]=="1":
            ensure(row["gene_a"] and row["gene_b"] and
                   row["gene_a"]!=row["gene_b"],
                   "actual accepted relation lacks valid genes")
        ensure(row["status"]=="OK" and
               row["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEV_ONLY",
               "non-passing or mislabeled study")
        records.append((d,f,method,w,ll-root_ll,brier-root_brier))
    ensure(len(records)==6,"pilot incomplete")
    for d in DS:
        vals=[r for r in records if r[0]==d]
        print(f"{d}: mean_delta_LL={statistics.mean(r[4] for r in vals):+.9f} "+
              "decisions="+",".join(
                 f"{r[1]}:{r[2]}@{r[3]:.2f}" for r in vals))
    print(f"RLMT5_MEAN_DELTA_LL_VS_ORIGINAL={statistics.mean(r[4] for r in records):+.9f}")
    print(f"RLMT5_ACTIVATED_INNER_VALIDATED_CORRECTIONS={sum(r[3]>0 for r in records)}/6")
    print("RLMT5_PILOT_INDEPENDENT_PATIENT_AUDIT=PASS")
    print("RLMT5_PILOT_ORIGINAL_RLR_PARENT_PARITY=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")

if __name__=="__main__":main()
