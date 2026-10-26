#!/usr/bin/env python3
"""Independent heldout patient prediction, complexity and 3-fold integrity audit."""
import csv,math,argparse
from pathlib import Path
from collections import defaultdict
n_map={"chen-2002":179,"singh-2002":102,"yeoh-2002-v1":248}
arms=("ORIGINAL_RLR","TSP_1","kTSP_9","TSPDT_GREEDY_PROXY",
      "ARLET_DETERMINISTIC","ARLET_EA")
def need(x,m):
 if not x:raise SystemExit("ARLET_INDEPENDENT_STOP: "+m)
def table(p):
 need(p.is_file(),"missing "+str(p))
 with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))
def score(y,p):
 pred=[int(v>=.5) for v in p]
 r=[sum(pred[j]==c for j in range(len(y)) if y[j]==c)/
    max(1,sum(v==c for v in y)) for c in (0,1)]
 f=[2*sum(pred[j]==c and y[j]==c for j in range(len(y)))/
    max(1,sum(v==c for v in pred)+sum(v==c for v in y))
    for c in (0,1)]
 prob=[max(1e-7,min(1-1e-7,v)) for v in p]
 return dict(BA=sum(r)/2,macro_F1=sum(f)/2,
   brier=sum((v-t)**2 for v,t in zip(p,y))/len(y),
   logloss=sum(-(t*math.log(v)+(1-t)*math.log1p(-v))
               for t,v in zip(y,prob))/len(y))
def main():
 q=argparse.ArgumentParser();q.add_argument("--dataset",choices=tuple(n_map),required=True)
 ds=q.parse_args().dataset
 root=Path("results/rlr_arlet_v1")/ds/"pilot"
 mets=table(root/"method_metrics.csv")
 preds=table(root/"patient_predictions.csv")
 splits=table(root/"split_manifest.csv")
 expect={(f"rep00_fold{i:02d}",m) for i in range(3) for m in arms}
 get={(r["split_id"],r["method"]):r for r in mets}
 need(len(mets)==18 and set(get)==expect,"method by outerfold missing/duplicates")
 grouped=defaultdict(list)
 for p in preds:grouped[p["split_id"],p["method"]].append(p)
 need(set(grouped)==expect,"patient groups incomplete")
 have=set();replayed=0
 for sp in splits:
  fold=sp["split_id"];ids=sp["E_patient_ids"].split("|")
  need(len(set(ids))==len(ids)==int(sp["n_E"]),"E patient IDs duplicated")
  need(not have.intersection(ids),"same patient in multiple E folds")
  have.update(ids)
  originals={x["sample_id"]:x for x in grouped[fold,"ORIGINAL_RLR"]}
  for arm in arms:
   rr=grouped[fold,arm]
   need(len(rr)==len(ids) and {x["sample_id"] for x in rr}==set(ids),
        "method E patient sets differ")
   truth=[];prob=[]
   for x in rr:
    y=int(x["truth"]);p=float(x["prob"])
    need(y in (0,1) and 0<=p<=1 and math.isfinite(p),
         "score or class invalid")
    need(y==int(originals[x["sample_id"]]["truth"]),
         "paired positive-class orientation mismatch")
    if arm.startswith("ARLET"):
     active=int(x["active_relations"])
     need(0<=active<=12,"too many active relationships per path")
     need(x["route"] in ("CONFIDENT_CLASS0","CONFIDENT_CLASS1","UNRESOLVED"),
          "missing interpretable 3-branch route")
    truth.append(y);prob.append(p)
   independent=score(truth,prob)
   ref=get[fold,arm]
   for metric,v in independent.items():
    need(abs(v-float(ref[metric]))<1e-8,
        f"metric disagreement {ds}/{fold}/{arm}/{metric}: {v}")
   need(ref["evidence_stage"]==
    "EXPOSED_PUBLICATION_SOURCE_RECONSTRUCTION_NOT_ARFF_SHA",
    "incorrect scientific evidence class")
   replayed+=len(rr)
 need(len(splits)==3 and len(have)==n_map[ds],
      "outer 3fold partitions incomplete")
 need(replayed==6*n_map[ds],"total independent patient predictions")
 print(f"ARLET_INDEPENDENT_DATASET={ds}")
 print(f"ARLET_INDEPENDENT_REPLAYED={replayed}")
 for arm in arms:
  x=[get[f"rep00_fold{i:02d}",arm] for i in range(3)]
  print("ARLET_INDEPENDENT_ARM",arm,
        f"BA={sum(float(v['BA']) for v in x)/3:.5f}",
        f"macroF1={sum(float(v['macro_F1']) for v in x)/3:.5f}",
        f"logloss={sum(float(v['logloss']) for v in x)/3:.5f}")
 print("ARLET_PATIENT_LEVEL_3_FOLDS_6_METHODS_REPLAYED=PASS")
 print("ARLET_E_PATIENTS_COVER_COHORT_WITHOUT_OVERLAP=PASS")
 print("ARLET_SOURCE_ARFF_SHA256_PARITY=NOT_VERIFIED")
if __name__=="__main__":main()
