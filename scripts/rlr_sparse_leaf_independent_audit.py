#!/usr/bin/env python3
"""Independent patient-level audit of genuine CompCancer original-RLR sparse pilot."""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path("results/rlr_sparse_leaf_real_pilot_v1")
SETS=("chowdary-2006","singh-2002")
METHODS=("ORIGINAL_RLR","GLOBAL_2","GLOBAL_5",
         "SHARED_2","SHARED_5","INDEPENDENT_2","INDEPENDENT_5")

def need(ok,msg):
    if not ok:raise SystemExit("SPARSE_AUDIT_STOP: "+msg)
def read(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))
def metric(y,p):
    need(0<=p<=1 and math.isfinite(p),"invalid predicted probability")
    p=min(1-1e-15,max(1e-15,p))
    return -(y*math.log(p)+(1-y)*math.log1p(-p))
def main():
    runs=read(ROOT/"run_manifest.csv")
    sm=read(ROOT/"method_metrics.csv")
    pr=read(ROOT/"patient_predictions.csv")
    methods={(r["dataset"],int(r["outer_fold"]),r["method"]):r for r in sm}
    expected={(ds,f,m) for ds in SETS for f in range(2) for m in METHODS}
    need(set(methods)==expected and len(sm)==len(expected)==28,
         "28 method-study rows missing or duplicates")
    expected_r={(ds,f) for ds in SETS for f in range(2)}
    need({(r["dataset"],int(r["outer_fold"])) for r in runs}==expected_r and
         len(runs)==4,"dataset/split manifest not complete")
    group=defaultdict(list)
    for x in pr:
        k=(x["dataset"],int(x["outer_fold"]),x["method"])
        need(k in expected,"unknown patient/method group")
        group[k].append(x)
    need(set(group)==expected,"missing patient/method block")
    recomputed=0
    for ds,fold in expected_r:
        source_n=104 if ds=="chowdary-2006" else 102
        run=next(x for x in runs
                 if x["dataset"]==ds and int(x["outer_fold"])==fold)
        need(int(run["n_train"])+int(run["n_E"])==source_n and
             int(run["D_n"])+int(run["V_n"])==int(run["n_train"]) and
             int(run["n_E"])>=30 and
             run["original_model_identity"]=="FROZEN_ORIGINAL_RLR_K10_M10_D2",
             "original-RLR identity/split invalid")
        base=group[ds,fold,"ORIGINAL_RLR"]
        ids={r["sample_id"] for r in base}
        need(len(ids)==len(base)==int(run["n_E"]),"baseline E patient IDs inconsistent")
        lookup={x["sample_id"]:x for x in base}
        for m in METHODS:
            key=(ds,fold,m)
            r=methods[key]
            a=group[key]
            need(len(a)==len(base) and
                 {x["sample_id"] for x in a}==ids and
                 len(set(x["sample_id"] for x in a))==len(a),
                 "heldout patient identity misaligned")
            losses=[]
            briers=[]
            gains=[]
            isaccepted=r["accepted"]=="1"
            for x in a:
                true=int(x["truth"])
                ref=lookup[x["sample_id"]]
                need(true in (0,1) and true==int(ref["truth"]),
                     "actual E phenotype inconsistent")
                old=float(x["original_prob"])
                parent=float(ref["method_prob"])
                p=float(x["method_prob"])
                need(old==parent,"original RLR patient-prediction parity violation")
                if m=="ORIGINAL_RLR" or not isaccepted:
                    need(p==parent,"unaccepted fallback does not equal RLR exactly")
                losses.append(metric(true,p))
                briers.append((p-true)**2)
                gains.append(metric(true,p)-metric(true,parent))
                recomputed+=1
            E_ll=sum(losses)/len(losses)
            E_brier=sum(briers)/len(briers)
            E_delta=sum(gains)/len(gains)
            for field,val in (("E_logloss",E_ll),
                              ("E_brier",E_brier),
                              ("delta_E_logloss",E_delta)):
                need(abs(float(r[field])-val)<1e-9,
                     "patient metric independently recomputed differs: "+str(key)+field)
            need(r["scientific_stage"]=="EXPOSED_BIOMEDICAL_DEV_ONLY",
                 "unwarranted scientific validation label")
    print("=== SPARSE RLR REAL PUBLIC BIOLOGICAL PILOT ===")
    for m in METHODS:
        rr=[methods[ds,fold,m] for ds,fold in sorted(expected_r)]
        accepted=sum(r["accepted"]=="1" for r in rr)
        d=sum(float(r["delta_E_logloss"]) for r in rr)/len(rr)
        print(f"{m:16s} V_accepted={accepted}/4 mean_E_deltaLL={d:+.6f}")
    print("SPARSE_REAL_28_METHOD_STUDIES_INDEPENDENTLY_AUDITED=PASS")
    print(f"SPARSE_REAL_{recomputed}_PATIENT_PREDICTION_METRICS_REPLAYED=PASS")
    print("SPARSE_REAL_EXACT_PATIENT_LEVEL_RLR_FALLBACK=PASS")
    print("SPARSE_REAL_NOT_INDEPENDENT_BIOMEDICAL_CONFIRMATION")
if __name__=="__main__":main()
