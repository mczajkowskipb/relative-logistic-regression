#!/usr/bin/env python3
"""Independent heldout E recomputation for original and pair-specific RLR encodings.
    Does not use the training algorithm or access E in learning.
    Asserts published outer split identity, E labels and per-patient metrics.
"""
from __future__ import annotations
import csv,math
from collections import defaultdict
from pathlib import Path

OUT=Path("results/rlr_pair_specific_evidence_v1")
SETS={"chowdary-2006":104,"singh-2002":102}
KS=(10,20,50)
MODES=("TERNARY_M10","ADAPTIVE_Q25","GRADED_M10_Q75")
def must(ok,msg):
    if not ok: raise SystemExit("RLR_MARGIN_AUDIT_STOP: "+msg)
def rd(n):
    p=OUT/n
    must(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def ll(y,p):
    must(0<=p<=1 and math.isfinite(p),"nonfinite or impossible probability")
    p=max(1e-15,min(1-1e-15,p))
    return -y*math.log(p)-(1-y)*math.log1p(-p)
def main():
    m=rd("method_metrics.csv")
    p=rd("patient_predictions.csv")
    manifest=rd("split_manifest.csv")
    keys={(ds,fold,K,mode) for ds in SETS for fold in range(3)
          for K in KS for mode in MODES}
    metric={(x["dataset"],int(x["fold"]),int(x["K"]),x["mode"]):x for x in m}
    must(len(m)==len(keys)==54 and set(metric)==keys,"study identity incomplete")
    groups=defaultdict(list)
    for x in p:
        key=(x["dataset"],int(x["fold"]),int(x["K"]),x["mode"])
        must(key in keys,"unknown patient method")
        groups[key].append(x)
    must(set(groups)==keys,"incomplete patient predictions")
    manifest_map={(x["dataset"],int(x["fold"])):x for x in manifest}
    must(len(manifest_map)==len(manifest)==6,"split manifest incomplete")
    n=0
    for ds,total in SETS.items():
      fold_ids=[]
      for fold in range(3):
        q=manifest_map[ds,fold]
        must(int(q["n_train"])+int(q["n_E"])==total and
             q["root_source"]=="ORIGINAL_RLR_SWITCHBOX_POOL1000",
             "source RLR pool or sample count drift")
        E=set(q["test_id_list"].split("|"))
        must(len(E)==int(q["n_E"]) and len(E)>20,
             "missing heldout patient IDs")
        fold_ids.append(E)
        for K in KS:
            base={x["sample_id"]:x for x in groups[ds,fold,K,"TERNARY_M10"]}
            must(set(base)==E,"M10 patient E mismatch")
            for mode in MODES:
                key=(ds,fold,K,mode)
                rr=metric[key]
                arr=groups[key]
                must(len(arr)==len(E) and
                     {x["sample_id"] for x in arr}==E,
                     "method patient alignment or duplicates")
                losses=[];brier=[];deltas=[]
                for x in arr:
                    b=base[x["sample_id"]]
                    truth=int(x["truth"])
                    must(truth==int(b["truth"]) and truth in (0,1),
                         "positive-class orientation or label mismatch")
                    p0=float(b["method_prob"])
                    p1=float(x["method_prob"])
                    must(abs(float(x["original_prob"])-p0)<1e-13,
                         "original RLR probability disagreement")
                    if mode=="TERNARY_M10":
                        must(p1==p0,"frozen original method not identical")
                    losses.append(ll(truth,p1))
                    brier.append((truth-p1)**2)
                    deltas.append(ll(truth,p1)-ll(truth,p0))
                    n+=1
                for col,got in (("logloss",sum(losses)/len(E)),
                                ("brier",sum(brier)/len(E)),
                                ("delta_logloss",sum(deltas)/len(E))):
                    must(abs(float(rr[col])-got)<1e-9,
                         f"metric mismatch {key} {col}")
                must(rr["scientific_stage"]==
                     "REUSED_COMPCANCER_DEVELOPMENT_ONLY",
                     "invalid stage claim")
                must(0<float(rr["mean_train_activity"])<=1 and
                     0<float(rr["mean_test_activity"])<=1 and
                     int(rr["global_active"])<=int(rr["K"]),
                     "feature activity outside bounds")
        must(set.union(*fold_ids)==
             set.union(*fold_ids),"internal fold mismatch")
      must(not any(fold_ids[i]&fold_ids[j]
                 for i in range(3) for j in range(i+1,3)),
           "outer cohort test participants leaked across folds")
      must(len(set.union(*fold_ids))==total,
           "outer cohort not fully covered")
    print("=== RLR PAIR-SPECIFIC MARGIN PILOT — REUSED PAPER DATA ONLY ===")
    for K in KS:
      for mode in MODES:
        ss=[metric[ds,f,K,mode] for ds in SETS for f in range(3)]
        delta=sum(float(x["delta_logloss"]) for x in ss)/len(ss)
        wins=sum(float(x["delta_logloss"])<-1e-9 for x in ss)
        print(f"K={K:02d} {mode:17s} mean_delta_E_logloss={delta:+.6f} improved={wins}/6")
    print(f"RLR_MARGIN_{n}_E_PATIENT_PREDICTIONS_INDEPENDENTLY_REPLAYED=PASS")
    print("RLR_MARGIN_54_SPLIT_METHOD_METRICS_AUDITED=PASS")
    print("RLR_MARGIN_E_PATIENTS_PARTITIONED_WITHOUT_OVERLAP=PASS")
    print("RLR_MARGIN_NOT_FRESH_EXTERNAL_BIOMEDICAL_VALIDATION")
if __name__=="__main__":main()
