#!/usr/bin/env python3
"""Independent E-patient validation of anchor lookahead RLR model stump V1."""
from __future__ import annotations
import csv,math
from collections import defaultdict
from pathlib import Path

OUT=Path("results/rlr_anchor_lookahead_child_rlr_v1")
COHORT={"chowdary-2006":104,"singh-2002":102}
ARMS=("ORIGINAL_RLR","GREEDY_ROOT_CHILD_RLR",
      "GREEDY_SAME_PAIRS_GLOBAL","GREEDY_GATE_ONLY",
      "LOOKAHEAD_CHILD_RLR","LOOKAHEAD_SAME_PAIRS_GLOBAL",
      "LOOKAHEAD_GATE_ONLY","V_GATED_LOOKAHEAD_CHILD_RLR")
def check(x,msg):
    if not x:raise SystemExit("LOOKAHEAD_AUDIT_STOP: "+msg)
def load(name):
    with (OUT/name).open(encoding="utf-8",newline="") as h:
        return list(csv.DictReader(h))
def logloss(y,p):
    q=max(1e-15,min(1-1e-15,p))
    return -y*math.log(q)-(1-y)*math.log1p(-q)
def values(y,p):
    pred=[int(z>=.5) for z in p]
    n0=y.count(0);n1=y.count(1)
    check(n0>0 and n1>0,"one E class missing")
    ba=.5*(sum(a==0 and b==0 for a,b in zip(y,pred))/n0+
           sum(a==1 and b==1 for a,b in zip(y,pred))/n1)
    fs=[]
    for c in (0,1):
        tp=sum(a==c and b==c for a,b in zip(y,pred))
        fp=sum(a!=c and b==c for a,b in zip(y,pred))
        fn=sum(a==c and b!=c for a,b in zip(y,pred))
        fs.append(2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0)
    return dict(logloss=sum(logloss(a,b) for a,b in zip(y,p))/len(y),
                brier=sum((a-b)**2 for a,b in zip(y,p))/len(y),
                ba=ba,macro_f1=sum(fs)/2)
def main():
    metrics=load("method_metrics.csv")
    pred=load("patient_predictions.csv")
    selections=load("selection_audit.csv")
    split=load("split_manifest.csv")
    keys={(d,f,a) for d in COHORT for f in range(3) for a in ARMS}
    md={(r["dataset"],int(r["fold"]),r["method"]):r for r in metrics}
    check(len(md)==len(metrics)==len(keys)==48 and set(md)==keys,
          "method-fold grid mismatch")
    rows=defaultdict(list)
    for r in pred:
        key=(r["dataset"],int(r["fold"]),r["method"])
        check(key in keys,"unexpected model/fold")
        rows[key].append(r)
    check(len(pred)==1648 and set(rows)==keys,
          "incomplete E patient predictions")
    sd={(r["dataset"],int(r["fold"])):r for r in selections}
    mf={(r["dataset"],int(r["fold"])):r for r in split}
    check(len(sd)==len(selections)==len(mf)==len(split)==6,
          "selection/split manifest duplicate")
    for ds,n in COHORT.items():
        seen_E=set()
        for f in range(3):
            manifest=mf[ds,f]
            T=set(manifest["T_ids"].split("|"))
            E=set(manifest["E_ids"].split("|"))
            check(len(T)==int(manifest["n_T"]) and
                  len(E)==int(manifest["n_E"]) and not T&E and
                  len(T|E)==n and not seen_E&E,
                  "E coverage / overlap / T contamination")
            seen_E.update(E)
            parent={r["sample_id"]:r for r in rows[ds,f,"ORIGINAL_RLR"]}
            check(len(parent)==len(E) and set(parent)==E,
                  "original E IDs mismatch")
            sel=sd[ds,f]
            accepted=int(sel["accepted"])
            check(accepted in (0,1),"invalid V acceptance")
            for arm in ARMS:
                rr=rows[ds,f,arm]
                rmap={r["sample_id"]:r for r in rr}
                check(len(rr)==len(rmap)==len(E) and set(rmap)==E,
                      "missing or repeated E patients")
                ys=[];ps=[];delta=[]
                for pid in sorted(E):
                    r=rmap[pid];b=parent[pid]
                    y=int(r["truth"]);p=float(r["prob"]);p0=float(b["prob"])
                    check(y in (0,1) and y==int(b["truth"]) and
                          math.isfinite(p) and 0<=p<=1 and
                          abs(float(r["original_prob"])-p0)<1e-13,
                          "E patient label, probability or parent parity")
                    if arm=="V_GATED_LOOKAHEAD_CHILD_RLR":
                        target=("LOOKAHEAD_CHILD_RLR" if accepted else
                                 "ORIGINAL_RLR")
                        comp=next(x for x in rows[ds,f,target]
                                  if x["sample_id"]==pid)
                        check(abs(p-float(comp["prob"]))<1e-13,
                              "V gating failed exact parent fallback")
                    ys.append(y);ps.append(p)
                    delta.append(logloss(y,p)-logloss(y,p0))
                stats=values(ys,ps)
                stats["delta_logloss"]=sum(delta)/len(delta)
                stored=md[ds,f,arm]
                for k,val in stats.items():
                    check(abs(float(stored[k])-val)<1e-9,
                          f"E metric replay mismatch {ds}/{f}/{arm}/{k}")
                check(int(stored["accepted"])==accepted and
                     stored["stage"]=="EXPOSED_RLR_PAPER_DEVELOPMENT_ONLY",
                     "method/fold stage or V gate discrepancy")
        check(len(seen_E)==n,"cohort E patient union incomplete")
    print("LOOKAHEAD_1648_E_PROBS_FULL_INDEPENDENT_REPLAY=PASS")
    print("LOOKAHEAD_48_METHOD_FOLD_MACROF1_BA_LOGLOSS_BRIER=PASS")
    print("LOOKAHEAD_PATIENT_SPLITS_AND_EXACT_RLR_FALLBACK=PASS")
    for arm in ARMS:
        a=[float(r["delta_logloss"]) for k,r in md.items()
           if k[2]==arm]
        print(f"LOOKAHEAD_AUDIT_SUMMARY {arm} mean_E_deltaLL={sum(a)/len(a):+.7f} E_wins={sum(x<-1e-9 for x in a)}/6")
    print("LOOKAHEAD_REUSED_COHORTS_NOT_FRESH_CONFIRMATION")
if __name__=="__main__": main()
