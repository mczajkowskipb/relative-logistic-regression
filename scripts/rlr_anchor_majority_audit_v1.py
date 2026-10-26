#!/usr/bin/env python3
"""Independent original-paper E-patient probability / fold / fallback audit."""
from __future__ import annotations
import csv, math
from collections import defaultdict
from pathlib import Path

OUT=Path("results/rlr_anchor_majority_tree_v1")
SETS={"chowdary-2006":104,"singh-2002":102}
ARMS=("ORIGINAL_RLR","ANCHOR_STUMP","ANCHOR_TREE_D2",
      "SAME_COMPARISONS_RIDGE","V_SELECTED_ANCHOR_TREE",
      "RLR_GATED_ANCHOR_TREE")
def check(flag,msg):
    if not flag:raise SystemExit("ANCHOR_AUDIT_STOP: "+msg)
def read(name):
    with (OUT/name).open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))
def ll(y,p):
    q=max(1e-15,min(1-1e-15,p))
    return -y*math.log(q)-(1-y)*math.log1p(-q)
def main():
    metrics=read("method_metrics.csv")
    preds=read("patient_predictions.csv")
    splits=read("split_manifest.csv")
    choices=read("selection_audit.csv")
    expected={(ds,f,a) for ds in SETS for f in range(3) for a in ARMS}
    met={(x["dataset"],int(x["fold"]),x["method"]):x for x in metrics}
    check(len(met)==len(metrics)==len(expected)==36 and set(met)==expected,
          "metrics cases missing or duplicated")
    rec=defaultdict(list)
    for x in preds:
        key=(x["dataset"],int(x["fold"]),x["method"])
        check(key in expected,"unknown arm")
        rec[key].append(x)
    check(len(preds)==1236 and set(rec)==expected,
          "wrong archived E probability inventory")
    smap={(x["dataset"],int(x["fold"])):x for x in splits}
    cmap={(x["dataset"],int(x["fold"])):x for x in choices}
    check(len(splits)==len(choices)==len(smap)==len(cmap)==6,
          "missing or repeated fold manifest")
    for ds,total in SETS.items():
        E_seen=set()
        for f in range(3):
            split=smap[ds,f]
            train_ids=set(split["training_ids"].split("|"))
            E_ids=set(split["E_ids"].split("|"))
            check(len(train_ids)==int(split["n_train"]) and
                  len(E_ids)==int(split["n_E"]) and
                  len(train_ids|E_ids)==total and not train_ids&E_ids,
                  "training/test ID leakage")
            check(not E_seen&E_ids,"duplicate E patients")
            E_seen.update(E_ids)
            base={r["sample_id"]:r for r in rec[ds,f,"ORIGINAL_RLR"]}
            check(set(base)==E_ids,"baseline ID mismatch")
            c=cmap[ds,f]
            accepted=int(c["V_accepted"])
            check(accepted in (0,1),"unknown gate")
            mode=c["inner_mode"]
            check(mode in ("NONE","STUMP","DEPTH2"),"invalid tree mode")
            for a in ARMS:
                group=rec[ds,f,a]
                rows={r["sample_id"]:r for r in group}
                check(len(group)==len(rows)==len(E_ids) and set(rows)==E_ids,
                      "missing/duplicated patient prediction")
                losses=[];br=[];truths=[];probs=[];delta=[]
                for sid in E_ids:
                    r=rows[sid]
                    y=int(r["truth"])
                    pp=float(r["prob"])
                    p0=float(base[sid]["prob"])
                    check(y in (0,1) and y==int(base[sid]["truth"]) and
                          math.isfinite(pp) and 0<=pp<=1 and
                          abs(float(r["original_prob"])-p0)<1e-13,
                          "truth/baseline class/parity issue")
                    p_sel=float(next(q["prob"] for q in rec[ds,f,
                       "V_SELECTED_ANCHOR_TREE"] if q["sample_id"]==sid))
                    p_stump=float(next(q["prob"] for q in rec[ds,f,
                       "ANCHOR_STUMP"] if q["sample_id"]==sid))
                    p_d2=float(next(q["prob"] for q in rec[ds,f,
                       "ANCHOR_TREE_D2"] if q["sample_id"]==sid))
                    if a=="V_SELECTED_ANCHOR_TREE":
                        target=p_d2 if mode=="DEPTH2" else p_stump if mode=="STUMP" else p0
                        check(abs(pp-target)<1e-13,"wrong V selected tree")
                    if a=="RLR_GATED_ANCHOR_TREE":
                        target=p_sel if accepted else p0
                        check(abs(pp-target)<1e-13,"V rejection fallback differs")
                    losses.append(ll(y,pp));br.append((y-pp)**2)
                    delta.append(ll(y,pp)-ll(y,p0))
                    truths.append(y);probs.append(pp)
                check(min(truths)==0 and max(truths)==1,"one class missing E")
                BA=.5*(sum(1 for y,p in zip(truths,probs)
                              if y==0 and p<.5)/truths.count(0)+
                       sum(1 for y,p in zip(truths,probs)
                              if y==1 and p>=.5)/truths.count(1))
                mm=met[ds,f,a]
                check(int(mm["accepted"])==accepted,
                      "V selection manifest mismatch")
                for col,v in (("logloss",sum(losses)/len(losses)),
                              ("brier",sum(br)/len(br)),
                              ("ba",BA),
                              ("delta_logloss",sum(delta)/len(delta))):
                    check(abs(float(mm[col])-v)<1e-9,
                          "metric mismatch "+str((ds,f,a,col)))
                check(mm["evidence_stage"]==
                      "EXPOSED_ORIGINAL_RLR_PAPER_DATA_DEVELOPMENT_ONLY",
                      "scientific stage mismatch")
        check(len(E_seen)==total,"incomplete cohort E coverage")
    print("ANCHOR_TREE_1236_PATIENT_E_PROBS_INDEPENDENTLY_REPLAYED=PASS")
    print("ANCHOR_TREE_36_METHOD_FOLD_METRICS_AUDITED=PASS")
    print("ANCHOR_TREE_EXACT_GATED_FALLBACK_AND_CLASS_PARITY=PASS")
    print("ANCHOR_TREE_ALL_COHORT_E_PATIENTS_PARTITIONED=PASS")
    for arm in ARMS:
        delta=[float(v["delta_logloss"]) for k,v in met.items()
               if k[2]==arm]
        print(f"ANCHOR_TREE_SUMMARY {arm} mean_deltaLL={sum(delta)/len(delta):+.7f} "
              f"E_wins={sum(x < -1e-9 for x in delta)}/6")
    print("ANCHOR_TREE_EXPOSED_DEVELOPMENT_NOT_FRESH_VALIDATION")
if __name__=="__main__":
    main()
