#!/usr/bin/env python3
"""Independent patient-level replay of OAS-DT V1 original-paper E output."""
from __future__ import annotations
import csv,math
from collections import defaultdict
from pathlib import Path

OUT=Path("results/rlr_ordinal_anchor_score_tree_v1")
SETS={"chowdary-2006":104,"singh-2002":102}
ARMS=("ORIGINAL_RLR","ANCHOR_SCORE_STUMP","ANCHOR_SCORE_TREE_D2",
      "SAME_ROOT_BEST_SINGLE_PAIR","SAME_COUNT_RIDGE",
      "V_SELECTED_ANCHOR_SCORE_TREE","RLR_V_GATED_ANCHOR_SCORE_TREE")
def check(cond,msg):
    if not cond:raise SystemExit("OASDT_AUDIT_STOP: "+msg)
def rd(name):
    with (OUT/name).open(encoding="utf-8",newline="") as h:
        return list(csv.DictReader(h))
def logloss(y,p):
    q=max(1e-15,min(1-1e-15,p))
    return -y*math.log(q)-(1-y)*math.log1p(-q)
def foldmetrics(y,p):
    pred=[int(x>=.5) for x in p]
    N0=y.count(0);N1=y.count(1)
    check(N0>0 and N1>0,"one E class absent")
    ba=.5*(sum(int(z==0 and v==0) for z,v in zip(y,pred))/N0+
           sum(int(z==1 and v==1) for z,v in zip(y,pred))/N1)
    f1=[]
    for c in (0,1):
        tp=sum(int(a==c and b==c) for a,b in zip(y,pred))
        fp=sum(int(a!=c and b==c) for a,b in zip(y,pred))
        fn=sum(int(a==c and b!=c) for a,b in zip(y,pred))
        f1.append(2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0)
    return dict(logloss=sum(logloss(a,b) for a,b in zip(y,p))/len(y),
                brier=sum((a-b)**2 for a,b in zip(y,p))/len(y),
                balanced_accuracy=ba,macro_f1=sum(f1)/2)
def main():
    met=rd("method_metrics.csv")
    pred=rd("patient_predictions.csv")
    choices=rd("selection_audit.csv")
    manifests=rd("split_manifest.csv")
    expected={(d,f,a) for d in SETS for f in range(3) for a in ARMS}
    m={(r["dataset"],int(r["fold"]),r["method"]):r for r in met}
    check(len(m)==len(met)==len(expected)==42 and set(m)==expected,
          "incomplete experiment arm grid")
    g=defaultdict(list)
    for r in pred:
        key=(r["dataset"],int(r["fold"]),r["method"])
        check(key in expected,"unexpected dataset or model")
        g[key].append(r)
    check(len(pred)==1442 and set(g)==expected,
          "missing method patient probability")
    c={(r["dataset"],int(r["fold"])):r for r in choices}
    v={(r["dataset"],int(r["fold"])):r for r in manifests}
    check(len(c)==len(choices)==len(v)==len(manifests)==6,
          "duplicate / missing V and outer manifests")
    for ds,total in SETS.items():
        coverage=set()
        for f in range(3):
            vrow=v[ds,f]
            T=set(vrow["training_ids"].split("|"))
            E=set(vrow["E_ids"].split("|"))
            check(not T&E and len(T)==int(vrow["n_train"]) and
                  len(E)==int(vrow["n_E"]) and len(T|E)==total,
                  "training/E contamination or source drift")
            check(not coverage&E,"outer E participant reuse")
            coverage.update(E)
            base={r["sample_id"]:r for r in g[ds,f,"ORIGINAL_RLR"]}
            check(len(base)==len(E) and set(base)==E,"original model E alignment")
            chosen=c[ds,f]["chosen"]
            accepted=int(c[ds,f]["accepted"])
            check(chosen in ("NONE","STUMP","DEPTH2") and accepted in (0,1),
                  "bad selected structure")
            for arm in ARMS:
                rr=g[ds,f,arm]
                check(len(rr)==len(E),"method E coverage count")
                arr={r["sample_id"]:r for r in rr}
                check(len(arr)==len(E) and set(arr)==E,
                      "duplicate/missing patient identifiers")
                yy=[]; pp=[]; dd=[]
                for pid in sorted(E):
                    r=arr[pid];rb=base[pid]
                    y=int(r["truth"]);p=float(r["prob"])
                    p0=float(rb["prob"])
                    check(y in (0,1) and y==int(rb["truth"]) and
                          math.isfinite(p) and 0<=p<=1 and
                          abs(float(r["original_prob"])-p0)<1e-13,
                          "patient class or parent probability mismatch")
                    if arm=="V_SELECTED_ANCHOR_SCORE_TREE":
                        requested={"STUMP":"ANCHOR_SCORE_STUMP",
                                   "DEPTH2":"ANCHOR_SCORE_TREE_D2",
                                   "NONE":"ORIGINAL_RLR"}[chosen]
                        other=next(t for t in g[ds,f,requested]
                                   if t["sample_id"]==pid)
                        check(abs(p-float(other["prob"]))<1e-13,
                              "selected tree has wrong depth")
                    if arm=="RLR_V_GATED_ANCHOR_SCORE_TREE":
                        mode="V_SELECTED_ANCHOR_SCORE_TREE" if accepted else "ORIGINAL_RLR"
                        other=next(t for t in g[ds,f,mode]
                                   if t["sample_id"]==pid)
                        check(abs(p-float(other["prob"]))<1e-13,
                              "rejected tree did not exactly fallback")
                    yy.append(y);pp.append(p)
                    dd.append(logloss(y,p)-logloss(y,p0))
                values=foldmetrics(yy,pp)
                values["delta_logloss"]=sum(dd)/len(dd)
                row=m[ds,f,arm]
                for name,value in values.items():
                    check(abs(float(row[name])-value)<1e-9,
                          f"metric mismatch {ds}/{f}/{arm}/{name}")
                check(int(row["accepted"])==accepted and
                      row["chosen"]==chosen and
                      row["stage"]=="EXPOSED_ORIGINAL_PAPER_DEVELOPMENT_ONLY",
                      "selection/evidence metadata drift")
        check(len(coverage)==total,"not all cohort patients tested")
    print("OASDT_1442_HELDOUT_PATIENT_SCORES_REPLAYED=PASS")
    print("OASDT_42_METHOD_FOLD_MACROF1_BA_LOGLOSS_BRIER_AUDITED=PASS")
    print("OASDT_ALL_PATIENTS_COVERED_AND_FALLBACK_EXACT=PASS")
    for arm in ARMS:
        cases=[float(r["delta_logloss"]) for k,r in m.items() if k[2]==arm]
        print(f"OASDT_SUMMARY {arm} mean_E_delta_logloss={sum(cases)/len(cases):+.7f} improved={sum(x<-1e-9 for x in cases)}/6")
    print("OASDT_EXPOSED_DEVELOPMENT_NOT_FRESH_VALIDATION")
if __name__=="__main__":main()
