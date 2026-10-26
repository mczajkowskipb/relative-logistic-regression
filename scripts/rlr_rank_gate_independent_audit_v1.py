#!/usr/bin/env python3
"""Independent heldout RLR rank-stump metric/parity audit; no training code."""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path

OUT=Path("results/rlr_rank_gate_v1")
N={"chowdary-2006":104,"singh-2002":102}
ARMS=("ORIGINAL_RLR","GLOBAL_MAIN_CONTROL",
      "RANK_STUMP_UNGATED","RANK_STUMP_V_GATED")
def fail(cond,msg):
    if not cond:
        raise SystemExit("RANK_GATE_AUDIT_STOP: "+msg)
def read(name):
    with (OUT/name).open(encoding="utf-8",newline="") as fh:
        return list(csv.DictReader(fh))
def ll(y,p):
    q=min(1-1e-15,max(1e-15,p))
    return -y*math.log(q)-(1-y)*math.log1p(-q)
def ba(rows):
    return (sum(int(int(r["prob"] and float(r["prob"])>=.5)==int(r["truth"]))
                for r in rows if int(r["truth"])==0) /
            sum(int(r["truth"])==0 for r in rows) +
            sum(int(float(r["prob"])>=.5) for r in rows if int(r["truth"])==1) /
            sum(int(r["truth"])==1 for r in rows))/2

def main():
    m=read("method_metrics.csv")
    p=read("patient_predictions.csv")
    v=read("split_manifest.csv")
    s=read("selection_audit.csv")
    expected={(ds,f,a) for ds in N for f in range(3) for a in ARMS}
    metric={(r["dataset"],int(r["fold"]),r["method"]):r for r in m}
    fail(len(metric)==len(m)==len(expected)==24 and set(metric)==expected,
         "unexpected method/fold matrix")
    grp=defaultdict(list)
    for row in p:
        k=row["dataset"],int(row["fold"]),row["method"]
        fail(k in expected,"unplanned method")
        grp[k].append(row)
    fail(set(grp)==expected and len(p)==824,"prediction inventory")
    selection={(r["dataset"],int(r["fold"])):r for r in s}
    manifest={(r["dataset"],int(r["fold"])):r for r in v}
    fail(len(s)==len(v)==len(selection)==len(manifest)==6,
         "selection and split manifest incomplete")
    for ds,total in N.items():
        coverage=set()
        for f in range(3):
            key=(ds,f)
            row=manifest[key]
            tr=set(row["train_ids"].split("|"))
            te=set(row["E_ids"].split("|"))
            fail(len(tr)==int(row["n_train"]) and
                 len(te)==int(row["n_E"]) and not tr&te,
                 "train/test contamination")
            fail(not coverage&te,"reused E patient")
            coverage.update(te)
            fail(len(tr|te)==total,"incomplete cohort split")
            base={r["sample_id"]:r for r in grp[ds,f,"ORIGINAL_RLR"]}
            fail(set(base)==te,"baseline E participants mismatch")
            accepted=int(selection[key]["accepted"])
            fail(accepted in (0,1),"invalid V gate")
            for a in ARMS:
                arr=grp[ds,f,a]
                seen={r["sample_id"]:r for r in arr}
                fail(len(arr)==len(seen)==len(te) and set(seen)==te,
                     "missing/duplicated E patient")
                ys=[];probs=[];losses=[];bs=[]
                for sid,r in seen.items():
                    y=int(r["truth"])
                    q=float(r["prob"])
                    q0=float(base[sid]["prob"])
                    fail(y in (0,1) and y==int(base[sid]["truth"]) and
                         0<=q<=1 and math.isfinite(q) and
                         abs(float(r["baseline_prob"])-q0)<1e-13,
                         "probability or label parity")
                    if a=="RANK_STUMP_V_GATED":
                        p_tree=float(grp[ds,f,"RANK_STUMP_UNGATED"][
                             next(i for i,z in enumerate(grp[ds,f,"RANK_STUMP_UNGATED"])
                                  if z["sample_id"]==sid)]["prob"])
                        target=p_tree if accepted else q0
                        fail(abs(q-target)<1e-13,
                             "V-gated model does not respect fallback")
                    ys.append(y);probs.append(q)
                    losses.append(ll(y,q));bs.append((y-q)**2)
                fail(min(ys)==0 and max(ys)==1,"fold lacks a class")
                mm=metric[ds,f,a]
                delta=sum(losses)/len(losses) - sum(
                  ll(int(base[r["sample_id"]]["truth"]),
                     float(base[r["sample_id"]]["prob"])) for r in arr)/len(arr)
                actual_ba=.5*(sum(int(q<.5) for q,y in zip(probs,ys) if y==0)/ys.count(0) +
                                  sum(int(q>=.5) for q,y in zip(probs,ys) if y==1)/ys.count(1))
                for col,val in (("logloss",sum(losses)/len(losses)),
                                ("brier",sum(bs)/len(bs)),
                                ("ba",actual_ba),("delta_logloss",delta)):
                    fail(abs(float(mm[col])-val)<1e-9,
                         "metric recomputation mismatch "+str((ds,f,a,col)))
                fail(int(mm["accepted"])==accepted,"selection recording")
                fail(mm["stage"]=="REUSED_ORIGINAL_PAPER_DATA_DEV_ONLY",
                     "invalid evidence stage")
        fail(len(coverage)==total,"entire cohort E coverage failed")
    print("RANK_GATE_824_E_PATIENT_ROWS_INDEPENDENTLY_REPLAYED=PASS")
    print("RANK_GATE_24_METHOD_FOLD_METRICS_AUDITED=PASS")
    print("RANK_GATE_E_COVERAGE_AND_FALLBACK=PASS")
    for a in ARMS:
        vals=[float(metric[k]["delta_logloss"]) for k in expected if k[2]==a]
        print(f"RANK_GATE_SUMMARY {a} mean_fold_delta_logloss={sum(vals)/len(vals):+.6f} improved={sum(x < -1e-9 for x in vals)}/6")
    print("RANK_GATE_EXPOSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY")
if __name__=="__main__":
    main()
