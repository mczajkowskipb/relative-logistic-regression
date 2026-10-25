#!/usr/bin/env python3
"""Stage F: independent 8,640 patient metric replay and pinned v9 oracle provenance.

This audits saved outputs. It cannot audit the original fitting math
independently of the R implementation; oracle source is fixed and synthetic.
"""
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path
from statistics import mean
from collections import defaultdict

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v9_stageF_xor_oracle_v1"
SCENARIOS=("null","xor_interaction")
METHODS=("ORIGINAL_RLR","DISCOVERED_CRIT9",
         "oracle_tree_unconditional","oracle_tree_V_accepted",
         "oracle_global_unconditional","oracle_global_V_accepted")

def need(ok,why):
    if not ok:raise SystemExit("STAGEF_AUDIT_STOP: "+why)

def load(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        rr=list(csv.DictReader(f))
    need(bool(rr),"empty "+str(path))
    return rr

def val(z):
    f=float(z)
    need(math.isfinite(f),"nonfinite statistic")
    return f

def metric(y,p):
    need(len(y)==len(p)>0 and all(0<=x<=1 for x in p),
         "bad patient probability")
    c=[max(1e-15,min(1-1e-15,v)) for v in p]
    ll=mean(-(z*math.log(q)+(1-z)*math.log1p(-q))
             for z,q in zip(y,c))
    br=mean((z-q)**2 for z,q in zip(y,p))
    return ll,br

def main():
    need(len(sys.argv)==2,"usage independent_audit.py PINNED_V9_WORKTREE")
    v9=Path(sys.argv[1]).resolve(strict=True)
    need((v9/"results/crit_v9_d_crossfit_2x3_v1/status.txt").read_text().strip()
          =="CRIT9_BATCH_EXIT=0","v9 archive incomplete")
    old=load(v9/"results/crit_v9_d_crossfit_2x3_v1/synthetic_panel_patients.csv")
    old_summary=load(v9/"results/crit_v9_d_crossfit_2x3_v1/synthetic_panel.csv")
    rr=load(OUT/"mechanism_test_patients.csv")
    summ=load(OUT/"mechanism_summary.csv")
    candidates=load(OUT/"D_candidate_provenance.csv")
    need(len(rr)==8640 and len(summ)==72 and len(old)==17280 and
         len(old_summary)==144,"unexpected output counts")
    expect={(s,n,m) for s in SCENARIOS for n in range(1,7) for m in METHODS}
    summary={(r["scenario"],int(r["replicate"]),r["method"]):r for r in summ}
    need(len(summary)==72 and set(summary)==expect,
         "duplicate/missing study method summaries")
    groups=defaultdict(list)
    for x in rr:
        key=(x["scenario"],int(x["replicate"]),x["method"])
        need(key in expect,"invalid study/method patient record")
        groups[key].append(x)
    need(set(groups)==expect,"missing patient blocks")
    old_group=defaultdict(list)
    for x in old:
        old_group[(x["scenario"],int(x["replicate"]),x["method"])].append(x)
    old_summary_map={(x["scenario"],int(x["replicate"]),x["method"]):x
                     for x in old_summary}
    cgroup=defaultdict(list)
    for r in candidates:
        key=(r["scenario"],int(r["replicate"]))
        need(key in {(s,n) for s in SCENARIOS for n in range(1,7)},
             "unknown D candidate study")
        cgroup[key].append(r)
    need(len(cgroup)==12,"missing D candidate provenance")
    for s in SCENARIOS:
      for rep in range(1,7):
        key=(s,rep,"ORIGINAL_RLR")
        ref=groups[key]
        need(len(ref)==120,"not 120 E patients")
        ids=[(x["sample_id"],x["y_true"]) for x in ref]
        y=[int(x["positive"]) for x in ref]
        need(len(set(ids))==120 and set(y)=={0,1},"test labels invalid")
        need(all(int(row["positive"])==int(row["y_true"]=="P")
                 for row in ref),"positive-class orientation invalid")
        p0=[val(x["prob"]) for x in ref]
        old_base=old_group[key]
        need([(x["sample_id"],x["y_true"]) for x in old_base]==ids and
             max(abs(a-val(x["prob"])) for a,x in zip(p0,old_base))<1e-10,
             "frozen v9 exact parent patient parity failed")
        ll0,_=metric(y,p0)
        sm=summary[key]
        source_old=old_summary_map[key]
        need(abs(ll0-val(source_old["logloss"]))<1e-10,
             "frozen baseline summary differs")
        c=cgroup[s,rep]
        ranks=sorted(int(x["D_rank"]) for x in c)
        need(1<=len(c)<=6 and ranks==list(range(1,len(c)+1)),
             "D candidate source incomplete")
        selected=[x for x in c if int(x["selected_by_D"])==1]
        need(len(selected)==1,"not precisely one D-selected candidate")
        selected_rank=int(selected[0]["D_rank"])
        samekeys=("pool_size","true_gate_in_pool","true_relation_in_pool",
                  "shortlisted_candidates","shortlisted_true_gate",
                  "shortlisted_true_gate_relation",
                  "shortlisted_true_both_leaves","selected_D_rank",
                  "selected_true_gate","selected_true_gate_relation",
                  "oracle_tree_V_accept","oracle_global_V_accept")
        for method in METHODS:
            key=(s,rep,method)
            r=summary[key]
            block=groups[key]
            need(len(block)==120 and
                 [(x["sample_id"],x["y_true"]) for x in block]==ids,
                 "E patients misaligned")
            need(all(r[k]==sm[k] for k in samekeys),
                 "method-level provenance discrepancies")
            need(int(r["selected_D_rank"])==selected_rank and
                 int(r["shortlisted_candidates"])==len(c) and
                 int(r["shortlisted_true_gate"])==
                    max(int(x["gate_is_planted"]) for x in c) and
                 int(r["shortlisted_true_gate_relation"])==
                    max(int(x["gate_interaction_is_planted"]) for x in c) and
                 int(r["shortlisted_true_both_leaves"])==
                    max(int(x["both_leaf_relations_planted"]) for x in c) and
                 int(r["selected_true_gate"])==
                    int(selected[0]["gate_is_planted"]) and
                 int(r["selected_true_gate_relation"])==
                    int(selected[0]["gate_interaction_is_planted"]),
                 "candidate provenance summary does not replay")
            prob=[val(x["prob"]) for x in block]
            base=[val(x["original_prob"]) for x in block]
            need(max(abs(a-b) for a,b in zip(base,p0))<1e-12,
                 "baseline patient probs misaligned")
            ll,br=metric(y,prob)
            need(abs(ll-val(r["logloss"]))<1e-9 and
                 abs(br-val(r["brier"]))<1e-9 and
                 abs(ll0-val(r["original_logloss"]))<1e-9 and
                 abs(ll-ll0-val(r["delta_vs_original"]))<1e-9,
                 "patient metrics failed independent replay")
            acc=int(r["accepted"])
            need(acc in (0,1) and
                 all(int(x["accepted"])==acc for x in block),
                 "unmatched acceptance flag")
            if method=="ORIGINAL_RLR":
                need(acc==0,"baseline cannot be accepted")
            if method=="DISCOVERED_CRIT9":
                archive=old_group[s,rep,"CONDITIONAL_TREE"]
                need(max(abs(a-val(b["prob"])) for a,b in
                         zip(prob,archive))<1e-10,
                     "CRIT9 v9 tree was not reproduced exactly")
            if method=="oracle_tree_V_accepted":
                need(acc==int(r["oracle_tree_V_accept"]) and
                     acc==int(val(r["oracle_tree_V_lower"])>0.002),
                     "oracle tree V test not one shot")
            if method=="oracle_global_V_accepted":
                need(acc==int(r["oracle_global_V_accept"]) and
                     acc==int(val(r["oracle_global_V_lower"])>0.002),
                     "oracle GLM V test not one shot")
            if method.endswith("_V_accepted") and not acc:
                need(max(abs(a-b) for a,b in zip(prob,p0))<1e-12,
                     "rejected oracle fails RLR fallback")
            if method.endswith("_unconditional"):
                need(acc==1,"unconditional oracle should always predict")
            need(r["scientific_stage"]==
                 "REUSED_SYNTHETIC_DEV_MECHANISM_ONLY",
                 "honest stage label missing")
        print(f"STAGEF_AUDIT {s}:{rep} pool_true="
              f"{sm['true_gate_in_pool']}/{sm['true_relation_in_pool']} "
              f"D_top6_gate={sm['shortlisted_true_gate']} "
              f"D_selected_true={sm['selected_true_gate_relation']} "
              f"oracle_tree_V_accept={sm['oracle_tree_V_accept']} "
              f"oracle_tree_E_delta={val(summary[s,rep,'oracle_tree_unconditional']['delta_vs_original']):+.6f}")
    for s in SCENARIOS:
        print("=== "+s+" ===")
        for method in METHODS:
            r=[summary[s,n,method] for n in range(1,7)]
            deltas=[val(x["delta_vs_original"]) for x in r]
            count=sum(int(x["accepted"]) for x in r)
            print(f"{method:34s} accept={count}/6 mean_deltaLL={mean(deltas):+.7f}")
    print("STAGEF_8640_INDEPENDENT_PATIENT_REPLAY=PASS")
    print("STAGEF_FROZEN_V9_PATIENT_PARITY=PASS")
    print("STAGEF_12_STUDY_DISCOVERY_RANK_PROVENANCE=PASS")
    print("STAGEF_SCIENTIFIC_STAGE=REUSED_SYNTHETIC_DEV_ONLY")

if __name__=="__main__":
    main()
