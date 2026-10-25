#!/usr/bin/env python3
"""Stage E: independent patient-level 24-simulation arm comparison audit."""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v9_d_crossfit_2x3_v1"
METHODS=("ORIGINAL_RLR","GLOBAL_CAL","GATE_ONLY",
         "ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")
SCENARIOS=("null","global_pair","gate_only","xor_interaction")

def need(ok,why):
    if not ok: raise SystemExit("CRIT9_SYNTH_AUDIT_STOP: "+why)

def read(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        a=list(csv.DictReader(f))
    need(bool(a),"empty "+str(path))
    return a

def num(v):
    v=float(v)
    need(math.isfinite(v),"nonfinite")
    return v

def calc(y,p):
    need(len(y)==len(p)>0,"missing patient data")
    clipped=[max(1e-15,min(1-1e-15,x)) for x in p]
    return (
       mean(-(a*math.log(b)+(1-a)*math.log1p(-b))
            for a,b in zip(y,clipped)),
       mean((a-b)**2 for a,b in zip(y,p))
    )

def main():
    summary=read(OUT/"synthetic_panel.csv")
    patients=read(OUT/"synthetic_panel_patients.csv")
    need(len(summary)==144 and len(patients)==17280,
         "expected 24 independent studies x 6 arms x 120 patients")
    k={(r["scenario"],int(r["replicate"]),r["method"]):r for r in summary}
    expected={(s,n,m) for s in SCENARIOS for n in range(1,7) for m in METHODS}
    need(len(k)==144 and set(k)==expected,
         "synthetic panel unexpected method inventory")
    groups=defaultdict(list)
    for r in patients:
        key=(r["scenario"],int(r["replicate"]),r["method"])
        need(key in expected,"extraneous patient")
        groups[key].append(r)
    need(set(groups)==expected,"missing model patient predictions")

    for scenario in SCENARIOS:
      for replicate in range(1,7):
        key=(scenario,replicate,"ORIGINAL_RLR")
        root=groups[key]
        ids=[(r["sample_id"],r["y_true"]) for r in root]
        yy=[int(r["y_positive"]) for r in root]
        need(len(root)==120 and len(set(ids))==120 and set(yy)=={0,1},
             "invalid synthetic outer-test patients")
        need(all((r["y_true"]=="P")==bool(int(r["y_positive"])) for r in root),
             "synthetic positive-class orientation mismatch")
        rootprob=[num(r["prob"]) for r in root]
        rootll,_=calc(yy,rootprob)
        selections=set()
        for method in METHODS:
            rsum=k[scenario,replicate,method]
            rpat=groups[scenario,replicate,method]
            need([(r["sample_id"],r["y_true"]) for r in rpat]==ids,
                 "outer synthetic patient alignment mismatch")
            probs=[num(r["prob"]) for r in rpat]
            orig=[num(r["original_prob"]) for r in rpat]
            need(all(0<=p<=1 for p in probs) and
                 max(abs(a-b) for a,b in zip(orig,rootprob))<1e-10,
                 "invalid probabilities or parent parity")
            gotll,gotbrier=calc(yy,probs)
            need(abs(gotll-num(rsum["logloss"]))<1e-9 and
                 abs(gotbrier-num(rsum["brier"]))<1e-9 and
                 abs(rootll-num(rsum["original_logloss"]))<1e-9 and
                 abs(gotll-rootll-num(rsum["delta_vs_original"]))<1e-9,
                 "patient metrics do not replay")
            accept=int(rsum["accepted"])
            need(accept in (0,1) and
                 all(int(r["accepted"])==accept for r in rpat),
                 "synthetic acceptance inconsistencies")
            if not accept:
                need(max(abs(a-b) for a,b in zip(probs,rootprob))<1e-12,
                     "rejected synthetic model not original RLR")
            nominee=int(rsum["n_D_candidates"])
            sel=rsum["selected_D_rank"]
            need(0<=nominee<=6 and
                 (sel in ("","NA") if not nominee else
                  sel not in ("","NA") and 1<=int(sel)<=nominee),
                 "D-only nomination inventory invalid")
            selections.add((nominee,sel))
            need(rsum["scientific_stage"]=="NEW_INDEPENDENT_SYNTHETIC_DEV_ONLY",
                 "development label absent")
        need(len(selections)==1,"D nominee differs across model arms")
    print("=== INDEPENDENT SYNTHETIC PATIENT REPLAY ===")
    for scenario in SCENARIOS:
        for method in METHODS[1:]:
            subset=[k[scenario,rep,method] for rep in range(1,7)]
            ds=[num(r["delta_vs_original"]) for r in subset]
            count=sum(int(r["accepted"]) for r in subset)
            print(f"{scenario:18s} {method:22s} "
                  f"accepted={count}/6 mean_deltaLL={mean(ds):+.6f} "
                  f"W/L/T={sum(d<-1e-10 for d in ds)}/"
                  f"{sum(d>1e-10 for d in ds)}/"
                  f"{sum(abs(d)<=1e-10 for d in ds)}")
    print("CRIT9_17280_SYNTHETIC_PATIENT_PREDICTIONS_AUDITED=PASS")
    print("CRIT9_24_FRESH_SYNTHETIC_STUDIES_AUDITED=PASS")
    print("CRIT9_SYNTHETIC_SCIENTIFIC_STAGE=DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
