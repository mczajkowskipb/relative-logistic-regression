#!/usr/bin/env python3
"""Per-dataset independent audit of original-paper full11 reconstructed-source CV.

Uses only written heldout patient scores, checks no patient overlap within
each 5-fold repetition; does NOT claim original processed ARFF hash parity.
"""
import argparse,csv,math
from pathlib import Path
from collections import defaultdict
ALL={"alizadeh-2000-v1":(42,10),"armstrong-2002-v1":(72,10),
 "bittner-2000":(38,10),"chen-2002":(179,10),
 "golub-1999-v1":(72,10),"gordon-2002":(181,10),
 "pomeroy-2002-v1":(34,10),"shipp-2002-v1":(77,10),
 "singh-2002":(102,10),"west-2001":(49,2),
 "yeoh-2002-v1":(248,2)}
MODES=("ORIGINAL_RLR","REGULARIZED_EQUAL_VOTE")
KS=(10,20,50)
def need(x,msg):
  if not x:raise SystemExit("FULL11_AUDIT_STOP: "+msg)
def csvread(p):
  need(p.is_file(),"missing "+str(p))
  with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def logloss(y,p):
  need(math.isfinite(p) and 0<=p<=1,"probability invalid")
  p=min(1-1e-15,max(1e-15,p))
  return -(y*math.log(p)+(1-y)*math.log1p(-p))
def main():
  parser=argparse.ArgumentParser()
  parser.add_argument("--dataset",required=True,choices=tuple(ALL))
  parser.add_argument("--mode",choices=("smoke","full"),required=True)
  a=parser.parse_args()
  n,reps=ALL[a.dataset]
  expected=5*reps if a.mode=="full" else 1
  root=Path("results/rlr_full11_same_pairs_v1")/a.dataset/a.mode
  mets=csvread(root/"method_metrics.csv")
  predictions=csvread(root/"patient_predictions.csv")
  splits=csvread(root/"split_manifest.csv")
  keys={(f"rep{r:02d}_fold{f:02d}",k,method)
        for r in range(reps) for f in range(5)
        for k in KS for method in MODES}
  if a.mode=="smoke":
    keys={("rep00_fold00",k,method) for k in KS for method in MODES}
  lookup={(x["split_id"],int(x["K"]),x["method"]):x for x in mets}
  need(len(mets)==len(lookup)==6*expected and
       set(lookup)==keys,"method-fold inventory is not frozen")
  splitid={x["split_id"]:x for x in splits}
  need(len(splits)==len(splitid)==expected,"missing split manifests")
  pred=defaultdict(list)
  for x in predictions:
    k=(x["split_id"],int(x["K"]),x["method"])
    need(k in keys,"unknown method-patient score")
    pred[k].append(x)
  need(set(pred)==keys,"missing prediction group")
  total_scored=0
  foldtest=defaultdict(list)
  for sid,sp in splitid.items():
    rep=int(sp["rep"]);fold=int(sp["fold"])
    need(sid==f"rep{rep:02d}_fold{fold:02d}" and
         int(sp["n_train"])+int(sp["n_E"])==n,
         "partition name, train/test sizes invalid")
    cohort_E=set(sp["test_ids"].split("|"))
    need(len(cohort_E)==int(sp["n_E"]),"duplicate heldout ID")
    foldtest[rep].append(cohort_E)
    for k in KS:
      src=pred[sid,k,"ORIGINAL_RLR"]
      srcids={x["sample_id"]:x for x in src}
      need(len(srcids)==len(src)==len(cohort_E) and
           set(srcids)==cohort_E,"original score/E identity mismatch")
      for method in MODES:
        entries=pred[sid,k,method]
        need(len(entries)==len(cohort_E) and
             {x["sample_id"] for x in entries}==cohort_E,
             "E scores not patient-aligned")
        m=lookup[sid,k,method]
        ls=[];br=[];delta=[]
        for x in entries:
          y=int(x["truth"])
          need(y in (0,1),"nonbinary class")
          parent=srcids[x["sample_id"]]
          need(y==int(parent["truth"]),"class labels disagree")
          p0=float(parent["method_prob"])
          p=float(x["method_prob"])
          need(abs(float(x["original_prob"])-p0)<1e-13,
               "score arrays not paired to identical source")
          if method=="ORIGINAL_RLR":
             need(p==p0,"literal original RLR score changed")
          ls.append(logloss(y,p))
          br.append((p-y)**2)
          delta.append(logloss(y,p)-logloss(y,p0))
          total_scored+=1
        for key,val in (("logloss",sum(ls)/len(ls)),
                        ("brier",sum(br)/len(br)),
                        ("delta_logloss",sum(delta)/len(delta))):
          need(abs(float(m[key])-val)<1e-9,
               f"patient level replay mismatch {sid}/{k}/{method}/{key}")
        need(m["scientific_stage"]=="SOURCE_RECONSTRUCTION_NOT_ARFF_PARITY",
             "misleading exact publication parity claim")
  for rep,folds in foldtest.items():
    need(all(not (a&b) for i,a in enumerate(folds) for b in folds[i+1:]),
         "patient present in more than one outer E fold within repetition")
    if a.mode=="full":
      need(len(folds)==5 and len(set.union(*folds))==n,
           "full 5fold CV incomplete")
  for k in KS:
    r=[lookup[sid,k,"REGULARIZED_EQUAL_VOTE"] for sid in splitid]
    d=sum(float(x["delta_logloss"]) for x in r)/len(r)
    print(f"FULL11_AUDIT dataset={a.dataset} K={k} mean_equal_minus_RLR_logloss={d:+.7f} RLR_better={sum(float(x['delta_logloss'])>0 for x in r)}/{len(r)}")
  print(f"FULL11_AUDIT_{a.dataset}_SPLITS={expected}")
  print(f"FULL11_AUDIT_E_PREDICTIONS_REPLAYED={total_scored}")
  print("FULL11_AUDIT_SPLITS_PATIENT_LEVEL_CONSISTENT=PASS")
  print("FULL11_AUDIT_ORIGINAL_ARFF_SHA_PARITY=NOT_VERIFIED")
if __name__=="__main__":main()
