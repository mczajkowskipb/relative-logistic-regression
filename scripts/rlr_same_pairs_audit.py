#!/usr/bin/env python3
"""Independent heldout E-patient score audit for same-pair equal vote vs RLR."""
import csv,math
from collections import defaultdict
from pathlib import Path
OUT=Path("results/rlr_same_pairs_weighting_v1")
DATASETS={"chowdary-2006":104,"singh-2002":102}
KS=(10,20,50)
METHODS=("ORIGINAL_RLR","CALIBRATED_EQUAL_VOTE")
def need(ok,msg):
    if not ok: raise SystemExit("SAME_PAIR_AUDIT_STOP: "+msg)
def load(filename):
    p=OUT/filename
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def logloss(y,p):
    need(math.isfinite(p) and 0<=p<=1,"invalid E probability")
    p=min(1-1e-15,max(1e-15,p))
    return -y*math.log(p)-(1-y)*math.log1p(-p)
def main():
    metrics=load("method_metrics.csv")
    preds=load("patient_predictions.csv")
    split=load("split_manifest.csv")
    keys={(ds,fold,k,m) for ds in DATASETS for fold in range(3)
          for k in KS for m in METHODS}
    idx={(x["dataset"],int(x["outer_fold"]),int(x["K"]),x["method"]):x
         for x in metrics}
    need(len(metrics)==len(keys)==36 and set(idx)==keys,
         "incorrect K/method/cohort split inventory")
    group=defaultdict(list)
    for x in preds:
        key=(x["dataset"],int(x["outer_fold"]),int(x["K"]),x["method"])
        need(key in keys,"unexpected E method block")
        group[key].append(x)
    need(set(group)==keys and len(preds)==1236,
         "E patient inventory incomplete")
    manifest={(x["dataset"],int(x["outer_fold"])):x for x in split}
    need(len(split)==6 and len(manifest)==6,
         "split manifest incorrect")
    audited=0
    for dataset,n in DATASETS.items():
      all_E=[]
      for fold in range(3):
        sp=manifest[dataset,fold]
        expected=set(sp["test_sample_ids"].split("|"))
        need(int(sp["n_train"])+int(sp["n_E"])==n and
             len(expected)==int(sp["n_E"]),"sample split mismatch")
        all_E.append(expected)
        for k in KS:
          original={x["sample_id"]:x for x in
                    group[dataset,fold,k,"ORIGINAL_RLR"]}
          need(set(original)==expected and len(original)==len(expected),
               "original patient identity mismatch")
          for method in METHODS:
            key=(dataset,fold,k,method)
            r=idx[key]
            ps=group[key]
            need(len(ps)==len(expected) and
                 {x["sample_id"] for x in ps}==expected,
                 "calibrated vote patient mismatch")
            losses=[];briers=[];delta=[]
            for x in ps:
              sample=x["sample_id"]
              truth=int(x["truth"])
              need(truth==int(original[sample]["truth"]) and truth in (0,1),
                   "E patient class inconsistency")
              p0=float(original[sample]["method_prob"])
              p1=float(x["method_prob"])
              need(abs(float(x["baseline_prob"])-p0)<1e-14,
                   "original score not identical across arms")
              if method=="ORIGINAL_RLR":
                  need(p0==p1,"baseline probability parity broken")
              losses.append(logloss(truth,p1))
              briers.append((p1-truth)**2)
              delta.append(logloss(truth,p1)-logloss(truth,p0))
              audited+=1
            for name,calculated in (("logloss",sum(losses)/len(ps)),
                                    ("brier",sum(briers)/len(ps)),
                                    ("delta_logloss",sum(delta)/len(ps))):
              need(abs(float(r[name])-calculated)<1e-9,
                   f"computed E metric discrepancy {key} {name}")
            need(r["stage"]=="EXPOSED_ORIGINAL_PAPER_COMPCANCER_DEV_ONLY",
                 "scientific stage not preserved")
      need(all(not (all_E[i]&all_E[j])
               for i in range(3) for j in range(i+1,3)),
           "overlapping E fold patient IDs")
      need(len(set.union(*all_E))==n,
           "patient E partitions do not cover cohort")
    print("=== EXACTLY MATCHED RLR PAIRS / CODING: CALIBRATED EQUAL WEIGHT VS RLR ===")
    for k in KS:
        rr=[idx[ds,f,k,"CALIBRATED_EQUAL_VOTE"]
            for ds in DATASETS for f in range(3)]
        d=sum(float(x["delta_logloss"]) for x in rr)/6
        rw=sum(float(x["delta_logloss"])>1e-9 for x in rr)
        print(f"K={k} EQUAL_minus_RLR_logloss={d:+.7f} RLR_better={rw}/6")
    print("SAME_PAIRS_36_METHOD_STUDIES_REPLAYED=PASS")
    print(f"SAME_PAIRS_{audited}_PATIENT_E_SCORES_REPLAYED=PASS")
    print("SAME_PAIRS_ALL_E_SPLITS_DISJOINT_AND_COMPLETE=PASS")
    print("SAME_PAIRS_EXPOSED_PUBLICATION_COHORTS_NOT_FRESH_VALIDATION")
if __name__=="__main__":main()
