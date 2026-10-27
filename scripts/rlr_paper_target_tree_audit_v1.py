#!/usr/bin/env python3
"""Independent matched kTSP / original RLR / greedy RLR tree diagnostic audit."""
from __future__ import annotations
import csv,math
from pathlib import Path
from collections import defaultdict
OUT=Path("results/rlr_paper_target_tree_v1")
SETS={"chowdary-2006":104,"singh-2002":102}
ARMS=("TSP1","KTSP10","KTSP20","KTSP50","RLR10","RLR20",
      "RLR50","GREEDY_RLR20_TREE","GREEDY_RLR50_TREE","VALUE_LR")
def ensure(c,m):
    if not c:raise SystemExit("PAPER_TARGET_AUDIT_STOP: "+m)
def load(p):
    with (OUT/p).open(newline="",encoding="utf-8") as h:
        return list(csv.DictReader(h))
def logloss(y,p):
    p=max(1e-15,min(1-1e-15,p))
    return -(y*math.log(p)+(1-y)*math.log1p(-p))
def auc(y,p):
    order=sorted(range(len(p)),key=lambda i:p[i])
    ranks=[0.0]*len(p)
    for start in range(len(p)):
        pass
    i=0
    while i<len(p):
        j=i+1
        while j<len(p) and p[order[j]]==p[order[i]]:j+=1
        r=(i+1+j)/2
        for k in range(i,j):ranks[order[k]]=r
        i=j
    pos=[i for i,v in enumerate(y) if v==1]
    neg=[i for i,v in enumerate(y) if v==0]
    ensure(pos and neg,"missing class")
    return (sum(ranks[i] for i in pos)-
            len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg))
def f1(y,h,c):
    tp=sum(a==c and b==c for a,b in zip(y,h))
    fp=sum(a!=c and b==c for a,b in zip(y,h))
    fn=sum(a==c and b!=c for a,b in zip(y,h))
    return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0
def compute(y,p,h):
    return dict(auc=auc(y,p),
      macro_f1=(f1(y,h,0)+f1(y,h,1))/2,
      positive_f1=f1(y,h,1),
      ba=(sum(a==0 and b==0 for a,b in zip(y,h))/y.count(0)+
          sum(a==1 and b==1 for a,b in zip(y,h))/y.count(1))/2,
      brier=sum((a-b)**2 for a,b in zip(y,p))/len(y),
      logloss=sum(logloss(a,b) for a,b in zip(y,p))/len(y))
def main():
    m=load("method_metrics.csv")
    p=load("patient_predictions.csv")
    f=load("split_manifest.csv")
    keys={(d,i,a) for d in SETS for i in range(3) for a in ARMS}
    metrics={(r["dataset"],int(r["fold"]),r["method"]):r for r in m}
    ensure(len(m)==len(metrics)==60 and set(metrics)==keys,
           "method-fold inventory")
    rows=defaultdict(list)
    for r in p:
        k=(r["dataset"],int(r["fold"]),r["method"])
        ensure(k in keys,"unplanned method")
        rows[k].append(r)
    ensure(len(p)==2060 and set(rows)==keys,
           "prediction row inventory")
    manifests={(r["dataset"],int(r["fold"])):r for r in f}
    ensure(len(f)==len(manifests)==6,"manifest cases")
    for ds,n in SETS.items():
        visited=set()
        for fold in range(3):
            entry=manifests[ds,fold]
            tr=set(entry["train_ids"].split("|"))
            te=set(entry["E_ids"].split("|"))
            ensure(len(tr)==int(entry["n_train"]) and
                   len(te)==int(entry["n_E"]) and
                   not tr&te and len(tr|te)==n and
                   not te&visited,"subject coverage/leakage")
            visited|=te
            ref={r["sample_id"]:r for r in rows[ds,fold,"RLR20"]}
            ensure(set(ref)==te,"RLR baseline alignment")
            for method in ARMS:
                pred=rows[ds,fold,method]
                byid={r["sample_id"]:r for r in pred}
                ensure(len(pred)==len(byid)==len(te) and set(byid)==te,
                       "duplicate or missing E patients")
                yy=[];pp=[];hh=[]
                for pid in sorted(te):
                    r=byid[pid]
                    y=int(r["truth"]);pr=float(r["prob"])
                    h=int(r["pred_label"])
                    ensure(y in (0,1) and h in (0,1) and
                           math.isfinite(pr) and 0<=pr<=1 and
                           y==int(ref[pid]["truth"]),
                           "patient class/probability mismatch")
                    yy.append(y);pp.append(pr);hh.append(h)
                val=compute(yy,pp,hh)
                row=metrics[ds,fold,method]
                for k,z in val.items():
                    ensure(abs(float(row[k])-z)<1e-9,
                           f"metric mismatch {ds}/{fold}/{method}/{k}")
                ensure(row["stage"]=="EXPOSED_COHORT_DEVELOPMENT_ONLY",
                       "evidence stage mislabeled")
        ensure(len(visited)==n,"incomplete E cohort coverage")
    print("PAPER_TARGET_2060_E_PATIENT_ROWS_REPLAYED=PASS")
    print("PAPER_TARGET_60_METHOD_FOLD_AUC_F1_BA_LOSS_BRIER_AUDITED=PASS")
    print("PAPER_TARGET_ALL_PATIENTS_E_COVERAGE_DISJOINT=PASS")
    for arm in ARMS:
        vals=[metrics[k] for k in keys if k[2]==arm]
        print(f"PAPER_TARGET_AUDIT_SUMMARY {arm:20s} AUC={sum(float(v['auc']) for v in vals)/6:.5f} MacroF1={sum(float(v['macro_f1']) for v in vals)/6:.5f} logloss={sum(float(v['logloss']) for v in vals)/6:.5f}")
    print("PAPER_TARGET_EXPOSED_TWO_COHORT_DIAGNOSTIC_ONLY")
if __name__=="__main__":main()
