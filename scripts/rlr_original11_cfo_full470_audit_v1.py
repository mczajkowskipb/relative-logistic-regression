#!/usr/bin/env python3
"""Independent original 11-cohort same-E patient/method audit from CI artifacts."""
import argparse,csv,math
from collections import defaultdict
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True)
ap.add_argument("--out",required=True)
ap.add_argument("--expected-splits",type=int,default=1)
ap.add_argument("--full470",action="store_true")
a=ap.parse_args()
cohorts=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000","chen-2002",
         "golub-1999-v1","gordon-2002","pomeroy-2002-v1","shipp-2002-v1",
         "singh-2002","west-2001","yeoh-2002-v1")
arms=("TSP1","KTSP10","KTSP20","KTSP50",
      "RLR10","RLR20","RLR50","VALUE_LR",
      "CFO_TREE_V_GATED","CFO_SAMEPAIRS_GLOBAL_V_GATED",
      "CFO_TREE_UNGATED_DIAGNOSTIC")
def check(ok,msg):
 if not ok:raise SystemExit("CFO_OUTER_AUDIT_STOP "+str(msg))
def read(path):
 with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def metric(yy,pp,hh):
 def f(cl):
  tp=sum(y==cl and h==cl for y,h in zip(yy,hh))
  fp=sum(y!=cl and h==cl for y,h in zip(yy,hh))
  fn=sum(y==cl and h!=cl for y,h in zip(yy,hh))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0
 n1=yy.count(1);n0=yy.count(0)
 check(n0>0 and n1>0,"no both E classes")
 ix=sorted(range(len(pp)),key=lambda i:pp[i])
 ranks=[0.]*len(pp)
 i=0
 while i<len(pp):
  j=i+1
  while j<len(pp) and pp[ix[i]]==pp[ix[j]]:j+=1
  q=(i+1+j)/2
  for k in range(i,j):ranks[ix[k]]=q
  i=j
 auc=(sum(ranks[k] for k in range(len(yy)) if yy[k]==1)-
      n1*(n1+1)/2)/(n1*n0)
 ll=sum(-(y*math.log(max(1e-15,min(1-1e-15,p)))+
           (1-y)*math.log1p(-max(1e-15,min(1-1e-15,p))))
       for y,p in zip(yy,pp))/len(yy)
 return {"positive_f1":f(1),"macro_f1":(f(0)+f(1))/2,
   "ba":(sum(y==0 and h==0 for y,h in zip(yy,hh))/n0+
         sum(y==1 and h==1 for y,h in zip(yy,hh))/n1)/2,
   "auc":auc,"logloss":ll,
   "brier":sum((y-p)**2 for y,p in zip(yy,pp))/len(yy)}
metrics=[];meta=[];allP=[]
E_by_rep=defaultdict(dict)
N_by_ds={}
for ds in cohorts:
 files=list(Path(a.input).rglob(f"{ds}/rep*/method_metrics.csv"))
 expected= (10 if ds in ("west-2001","yeoh-2002-v1") else 50) if a.full470 else a.expected_splits
 check(len(files)==expected,
       f"{ds} missing repeated-fold archives expected {expected} got {len(files)}")
 for file in files:
  root=file.parent
  mm=read(file)
  m=read(root/"split_manifest.csv")
  pp=read(root/"patient_predictions.csv")
  check(len(m)==1 and len(mm)==11, f"manifest/model count {ds}")
  ma=m[0];sid=ma["split_id"]
  stage="ORIGINAL_PAPER_NAMED_PUBLIC_SOURCE_RECONSTRUCTION_DEVELOPMENT" if a.full470 else "ALREADY_EXPLORED_DATA_FIRST_OUTER_E_DEVELOPMENT"
  check(ma["dataset"]==ds and ma["evidence"]==stage,
        f"bad evidence stage {ds}")
  T=ma["T_ids"].split("|");E=ma["E_ids"].split("|")
  N_by_ds[ds]=int(ma["n_total"])
  rep=int(sid.split("_")[0].replace("rep",""))
  fold=int(sid.split("_")[1].replace("fold",""))
  check(fold in range(5) and rep in range(2 if ds in ("west-2001","yeoh-2002-v1") else 10),
        "split outside original 470 plan")
  check(fold not in E_by_rep[ds,rep],"duplicate original outer fold")
  E_by_rep[ds,rep][fold]=set(E)
  check(len(T)==int(ma["n_outer_T"]) and len(set(T))==len(T) and
        len(E)==int(ma["n_outer_E"]) and len(set(E))==len(E) and
        not (set(T)&set(E)) and len(T)+len(E)==int(ma["n_total"]),
        f"outer T/E leak or missing subject {ds}")
  bm={r["method"]:r for r in mm}
  check(set(bm)==set(arms) and len(pp)==len(E)*len(arms),
        f"wrong original 11 baselines or row count {ds}")
  preds=defaultdict(list)
  for p in pp:preds[p["method"]].append(p)
  check(set(preds)==set(arms),f"missing score arm {ds}")
  for name in arms:
   ar=preds[name];pr={r["sample_id"]:r for r in ar}
   check(len(pr)==len(E) and set(pr)==set(E),f"prediction IDs {ds}/{name}")
   p=[];y=[];h=[]
   for patient in sorted(E):
    row=pr[patient]
    check(row["dataset"]==ds and row["split_id"]==sid,
         f"prediction wrong split {ds}/{name}")
    yi=int(row["truth"]);pi=float(row["prob"]);hi=int(row["label"])
    check(yi in (0,1) and hi in (0,1) and
          math.isfinite(pi) and 0<=pi<=1,f"invalid patient score {ds}/{name}")
    y.append(yi);p.append(pi);h.append(hi)
   val=metric(y,p,h)
   for key,z in val.items():
    check(abs(float(bm[name][key])-z)<1e-9,
          f"metric mismatch {ds}/{sid}/{name}/{key}")
   allP.append((ds,sid,name,len(p)))
  ref=preds["RLR50"]
  for arm in ("CFO_TREE_V_GATED","CFO_SAMEPAIRS_GLOBAL_V_GATED"):
   refp={r["sample_id"]:r for r in ref}
   cmp={r["sample_id"]:r for r in preds[arm]}
   if ma["accepted"]=="0":
    check(all(float(cmp[k]["prob"])==float(refp[k]["prob"]) and
              cmp[k]["label"]==refp[k]["label"] for k in E),
          f"no exact original RLR fallback {ds}/{arm}")
  if ma["accepted"]=="1":
   check(float(ma["V_parent_logloss_gain"])>=.005 and
         float(ma["V_global_logloss_gain"])>=.002 and
         float(ma["V_ba_change"])>=-.02 and
         float(ma["V_macro_f1_change"])>=-.01,
         f"gate accepted without train evidence {ds}")
  meta.append(ma);metrics.extend(mm)
if a.full470:
 for ds in cohorts:
  repetitions=2 if ds in ("west-2001","yeoh-2002-v1") else 10
  for rep in range(repetitions):
   parts=E_by_rep[ds,rep]
   check(set(parts)==set(range(5)),"missing original 5fold CV repeat "+ds)
   taken=set()
   for i in range(5):
    check(not(taken & parts[i]),"same patient tested twice per repeat "+ds)
    taken |= parts[i]
   check(len(taken)==N_by_ds[ds],"not full patient coverage per repeat "+ds)
 check(len(meta)==470 and len(metrics)==5170 and
       len(allP)==5170,"wrong full original 470 method count")
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all_cohort_method_fold_metrics.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=metrics[0].keys());w.writeheader();w.writerows(metrics)
with (out/"all_cohort_selection_manifest.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=meta[0].keys());w.writeheader();w.writerows(meta)
print("CFO_OUTER_11_COHORTS_E_MANIFEST_AND_NATIVE_KTSP_LABELS=PASS")
print("CFO_OUTER_"+str(len(metrics))+"_MODEL_FOLD_METRICS_INDEPENDENT_REPLAY=PASS")
print("CFO_OUTER_"+str(len(allP))+"_METHOD_FOLD_PREDICTION_BLOCKS_VERIFIED=PASS")
print("CFO_OUTER_PROBABILITIES_AND_V_GATED_EXACT_FALLBACK=PASS")
print("CFO_OUTER_TRAIN_ONLY_STRUCTURE_NO_ORIGINAL_ARFF_PARITY_CLAIM")
for arm in arms:
 rows=[x for x in metrics if x["method"]==arm]
 print(f"CFO_OUTER_ARM {arm} F1_pos={sum(float(z['positive_f1']) for z in rows)/len(rows):.6f} MacroF1={sum(float(z['macro_f1']) for z in rows)/len(rows):.6f} AUC={sum(float(z['auc']) for z in rows)/len(rows):.6f} Logloss={sum(float(z['logloss']) for z in rows)/len(rows):.6f}")
print("CFO_OUTER_INNER_V_ACCEPTED="+str(sum(z["accepted"]=="1" for z in meta))+"/"+str(len(meta)))
if a.full470:
 perds=[]
 for ds in cohorts:
  for arm in arms:
   these=[r for r in metrics if r["dataset"]==ds and r["method"]==arm]
   check(len(these)==(10 if ds in ("west-2001","yeoh-2002-v1") else 50),
         "dataset-method fold inventory mismatch "+ds+"/"+arm)
   item={"dataset":ds,"method":arm,"n_outer_folds":len(these)}
   for k in ("positive_f1","macro_f1","ba","auc","logloss","brier"):
    item[k]=sum(float(r[k]) for r in these)/len(these)
   perds.append(item)
 with (out/"all11_dataset_level_means.csv").open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=perds[0].keys());w.writeheader();w.writerows(perds)
 print("CFO_OUTER_ORIGINAL_470_VALIDATION_FOLDS=PASS")
 print("CFO_OUTER_EACH_DATASET_REP_PATIENT_COVERAGE_DISJOINT=PASS")
 for arm in arms:
  v=[r for r in perds if r["method"]==arm]
  print("CFO_FULL470_DATASET_MEAN "+arm+
        " posF1="+format(sum(r["positive_f1"] for r in v)/11,".6f")+
        " macroF1="+format(sum(r["macro_f1"] for r in v)/11,".6f")+
        " auc="+format(sum(r["auc"] for r in v)/11,".6f")+
        " logloss="+format(sum(r["logloss"] for r in v)/11,".6f"))

