#!/usr/bin/env python3
"""Independent 34,560 patient audit of 32-gene complete-pool ablation."""
from __future__ import annotations
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v9_stageG_exhaustive32_v1"
SCENARIOS=("null","global_pair","gate_only","xor_interaction")
GROUPS=("REUSED_DEV","FRESH_DEV")
METHODS=("ORIGINAL_RLR","GLOBAL_CAL","GATE_ONLY",
         "ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")

def check(ok,why):
    if not ok:raise SystemExit("STAGEG_AUDIT_STOP: "+why)

def load(path):
    check(path.is_file(),"missing file "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    check(bool(rows),"empty file "+str(path))
    return rows

def f(z):
    z=float(z)
    check(math.isfinite(z),"nonfinite numerical")
    return z

def metrics(y,p):
    check(len(y)==len(p)==120 and all(0<=x<=1 for x in p),
          "bad study patient probabilities")
    pc=[max(1e-15,min(1-1e-15,x)) for x in p]
    return (mean(-(z*math.log(q)+(1-z)*math.log1p(-q))
                 for z,q in zip(y,pc)),
            mean((z-q)**2 for z,q in zip(y,p)))

def main():
    check(len(sys.argv)==2,"usage: auditor.py CRIT9_FROZEN_WORKTREE")
    v9=Path(sys.argv[1]).resolve(strict=True)
    oldroot=v9/"results/crit_v9_d_crossfit_2x3_v1"
    check((oldroot/"status.txt").read_text().strip()=="CRIT9_BATCH_EXIT=0",
          "frozen v9 not completed")
    old=load(oldroot/"synthetic_panel_patients.csv")
    summaries=load(OUT/"study_summary.csv")
    patients=load(OUT/"test_patients.csv")
    candidates=load(OUT/"candidate_provenance.csv")
    check(len(summaries)==288 and len(patients)==34560 and len(old)==17280,
          "48 study by six methods inventory incomplete")
    expected={(s,n,m) for s in SCENARIOS for n in range(1,13) for m in METHODS}
    m={(x["scenario"],int(x["replicate"]),x["method"]):x for x in summaries}
    check(len(m)==288 and set(m)==expected,"missing or duplicate results")
    patientblocks=defaultdict(list)
    for x in patients:
        k=(x["scenario"],int(x["replicate"]),x["method"])
        check(k in expected,"unexpected patient group")
        patientblocks[k].append(x)
    check(set(patientblocks)==expected,"missing patient result")
    oldblocks=defaultdict(list)
    for x in old:
        oldblocks[(x["scenario"],int(x["replicate"]),x["method"])].append(x)
    nominated=defaultdict(list)
    for x in candidates:
        k=(x["scenario"],int(x["replicate"]))
        check(k[:1][0] in SCENARIOS and k[1] in range(1,13),
              "unknown candidate")
        nominated[k].append(x)
    check(len(nominated)==48,"missing D candidate provenance")
    for s in SCENARIOS:
      for rep in range(1,13):
        k=(s,rep,"ORIGINAL_RLR")
        root=patientblocks[k]
        check(len(root)==120,"study not 120 patients")
        ids=[(x["sample_id"],x["y_true"]) for x in root]
        check(len(set(ids))==120,"duplicated patient")
        y=[int(x["positive"]) for x in root]
        check(set(y)=={0,1} and
              all(z==int(x["y_true"]=="P") for x,z in zip(root,y)),
              "positive label orientation incorrect")
        orig=[f(x["prob"]) for x in root]
        ll0,_=metrics(y,orig)
        g="REUSED_DEV" if rep<=6 else "FRESH_DEV"
        if rep<=6:
            archived=oldblocks[k]
            check([(x["sample_id"],x["y_true"]) for x in archived]==ids and
                  max(abs(a-f(b["prob"])) for a,b in zip(orig,archived))<1e-10,
                  "exact original RLR synthetic patient parity failed")
        cs=sorted(nominated[s,rep],key=lambda x:int(x["D_rank"]))
        check(1<=len(cs)<=6 and
              [int(x["D_rank"]) for x in cs]==list(range(1,len(cs)+1)),
              "shortlist rank sequence incorrect")
        chosen=[x for x in cs if int(x["chosen_by_D"])==1]
        check(len(chosen)==1 and
              int(chosen[0]["D_rank"]) ==
              1+max(range(len(cs)),key=lambda i:f(cs[i]["D_cv_lower"])),
              "D-only rank maximum inconsistency")
        for method in METHODS:
            r=m[s,rep,method]
            block=patientblocks[s,rep,method]
            check(len(block)==120 and r["group"]==g and
                  all(x["group"]==g for x in block),
                  "synthetic age/group inventory wrong")
            check([(x["sample_id"],x["y_true"]) for x in block]==ids,
                  "patient labels/IDs don't match")
            prob=[f(x["prob"]) for x in block]
            base=[f(x["original_prob"]) for x in block]
            check(max(abs(a-b) for a,b in zip(base,orig))<1e-12,
                  "baseline patients do not align")
            ll,br=metrics(y,prob)
            check(abs(ll-f(r["logloss"]))<1e-9 and
                  abs(br-f(r["brier"]))<1e-9 and
                  abs(ll0-f(r["original_logloss"]))<1e-9 and
                  abs(ll-ll0-f(r["delta_vs_original"]))<1e-9,
                  "test metrics not independently reproducible")
            acc=int(r["accepted"])
            check(acc in (0,1) and all(int(x["accepted"])==acc for x in block),
                  "acceptance inconsistency")
            if not acc:
                check(all(a==b for a,b in zip(prob,orig)),
                      "rejected model does not equal original RLR exactly")
            check(int(r["exhaustive_pair_count"])==496 and
                  int(r["exhaustive_true_gate_present"])==1 and
                  int(r["exhaustive_true_relation_present"])==1 and
                  int(r["original_pool_size"])<496 and
                  int(r["shortlisted_n"])==len(cs) and
                  int(r["D_selected_rank"])==int(chosen[0]["D_rank"]),
                  "complete pool/nomination invariant failed")
            check(int(r["shortlisted_planted_gate"]) ==
                  max(int(x["planted_gate"]) for x in cs) and
                  int(r["shortlisted_planted_joint"]) ==
                  max(int(x["planted_joint"]) for x in cs) and
                  int(r["shortlisted_planted_both_leaves"]) ==
                  max(int(x["planted_both_leaves"]) for x in cs) and
                  int(r["D_selected_planted_gate"]) ==
                  int(chosen[0]["planted_gate"]) and
                  int(r["D_selected_planted_joint"]) ==
                  int(chosen[0]["planted_joint"]),
                  "shortlist and planted truth provenance mismatch")
            check(r["scientific_stage"]=="SYNTHETIC_DEVELOPMENT_ONLY",
                  "incorrect scientific stage")
    print("=== EXHAUSTIVE32 SELECTION SUMMARY ===")
    for group in GROUPS:
      for scenario in SCENARIOS:
        rr=[m[scenario,n,"CONDITIONAL_TREE"]
            for n in range(1,13)
            if ("REUSED_DEV" if n<=6 else "FRESH_DEV")==group]
        check(len(rr)==6,"missing panel group")
        print(f"{group:10s} {scenario:17s} "
              f"orig_pool_gate={sum(int(x['original_gate_present']) for x in rr)}/6 "
              f"proposed_gate={sum(int(x['shortlisted_planted_gate']) for x in rr)}/6 "
              f"proposed_joint={sum(int(x['shortlisted_planted_joint']) for x in rr)}/6 "
              f"selected_joint={sum(int(x['D_selected_planted_joint']) for x in rr)}/6 "
              f"tree_accept={sum(int(x['accepted']) for x in rr)}/6 "
              f"mean_deltaLL={mean(f(x['delta_vs_original']) for x in rr):+.7f}")
    print("STAGEG_34560_PATIENT_METRICS_INDEPENDENTLY_AUDITED=PASS")
    print("STAGEG_48_STUDY_NOMINATION_AND_FALLBACK_AUDITED=PASS")
    print("STAGEG_24_REUSED_STUDIES_EXACT_ORIGINAL_RLR_PARITY=PASS")
    print("STAGEG_SCIENTIFIC_STAGE=SYNTHETIC_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
