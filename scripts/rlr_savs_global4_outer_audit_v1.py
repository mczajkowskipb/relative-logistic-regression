#!/usr/bin/env python3
"""Independent patient-replay for original CompCancer firstfold or 470 E folds."""
import argparse,csv,math,collections
from pathlib import Path
a=argparse.ArgumentParser()
a.add_argument("--input",required=True)
a.add_argument("--out",required=True)
a.add_argument("--full470",action="store_true")
q=a.parse_args()
sets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
arms=("TSP1","KTSP4","KTSP10","KTSP20","KTSP50",
      "RLR10","RLR20","RLR50","VALUE_LR",
      "SELECTIVE4_M10","SELECTIVE4_M10PCT")
def need(flag,msg):
 if not flag:raise SystemExit("SAVS470_AUDIT_STOP "+msg)
def read(path):
 with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def score(y,p,h):
 n1=y.count(1);n0=y.count(0)
 need(min(n1,n0)>0,"empty E class")
 def f(cl):
  tp=sum(x==cl and z==cl for x,z in zip(y,h))
  fp=sum(x!=cl and z==cl for x,z in zip(y,h))
  fn=sum(x==cl and z!=cl for x,z in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 order=sorted(range(len(p)),key=lambda i:p[i])
 rank=[0.]*len(p);i=0
 while i<len(p):
  j=i+1
  while j<len(p) and p[order[j]]==p[order[i]]:j+=1
  r=(i+1+j)/2
  for k in range(i,j):rank[order[k]]=r
  i=j
 auc=(sum(rank[i] for i in range(len(y)) if y[i]==1)-
      n1*(n1+1)/2)/(n1*n0)
 cp=[max(1e-15,min(1-1e-15,z)) for z in p]
 return dict(positive_f1=f(1),macro_f1=(f(0)+f(1))/2,
    ba=(sum(x==0 and z==0 for x,z in zip(y,h))/n0+
        sum(x==1 and z==1 for x,z in zip(y,h))/n1)/2,
    auc=auc,
    logloss=sum(-(x*math.log(z)+(1-x)*math.log1p(-z))
                for x,z in zip(y,cp))/len(y),
    brier=sum((x-z)**2 for x,z in zip(y,p))/len(y))
full=[];manifests=[];reps=collections.defaultdict(dict)
for ds in sets:
 folders=list(Path(q.input).rglob(ds+"/rep*/method_metrics.csv"))
 expected=(10 if ds in ("west-2001","yeoh-2002-v1") else 50) if q.full470 else 1
 need(len(folders)==expected, f"expected {expected} folds for {ds}, got {len(folders)}")
 for src in folders:
  folder=src.parent
  rows=read(src)
  pat=read(folder/"patient_predictions.csv")
  meta=read(folder/"manifest.csv")
  need(len(rows)==11 and len(meta)==1 and
        {r["method"] for r in rows}==set(arms),
        "missing methods or manifest "+ds)
  m=meta[0];sid=m["split_id"]
  need(m["dataset"]==ds and
      m["evidence"]=="ALL_FEATURE_SELECTION_OUTER_T_ONLY_ORIGINAL_E_TEST_ONLY",
      "manifest evidence mismatch "+ds)
  T=m["T_ids"].split("|");E=m["E_ids"].split("|")
  need(len(T)==int(m["n_T"]) and len(E)==int(m["n_E"]) and
      len(set(T))==len(T) and len(set(E))==len(E) and
      not(set(T)&set(E)) and len(T)+len(E)==int(m["n_total"]),
      "patient leakage or incomplete T/E "+ds)
  rep=int(sid[3:5]);fold=int(sid[-2:])
  need(0<=fold<=4 and 0<=rep<(2 if ds in ("west-2001","yeoh-2002-v1") else 10),
      "fold id out of original scope")
  need(fold not in reps[ds,rep],"duplicate cohort fold")
  reps[ds,rep][fold]=set(E)
  groups=collections.defaultdict(list)
  for row in pat:
   need(row["dataset"]==ds and row["split_id"]==sid,"mixed dataset/split")
   groups[row["method"]].append(row)
  need(set(groups)==set(arms),"missing per-patient method "+ds)
  for method in arms:
   z=next(t for t in rows if t["method"]==method)
   pp=groups[method]
   need(len(pp)==len(E) and {r["sample_id"] for r in pp}==set(E),
        "test patient mismatch "+ds+"/"+method)
   per={r["sample_id"]:r for r in pp}
   ordered=[per[id] for id in sorted(E)]
   y=[int(r["truth"]) for r in ordered]
   p=[float(r["prob"]) for r in ordered]
   h=[int(r["label"]) for r in ordered]
   need(all(yi in (0,1) and hi in (0,1) and
        math.isfinite(pi) and 0<=pi<=1
        for yi,hi,pi in zip(y,h,p)),"invalid probability")
   evals=score(y,p,h)
   for key,val in evals.items():
    need(abs(float(z[key])-val)<1e-9,
         f"metric mismatch {ds}/{sid}/{method}/{key}")
   if method.startswith("SELECTIVE4"):
    aa=[int(r["n_active"]) for r in ordered]
    bb=[int(r["n_tested"]) for r in ordered]
    need(all(0<=a<=b<=4 for a,b in zip(aa,bb)),
         "abstention activity inconsistent "+ds)
    need(abs(sum(aa)/len(aa)-float(z["mean_active"]))<1e-9 and
         abs(sum(bb)/len(bb)-float(z["mean_tested"]))<1e-9,
         "activity mean mismatch "+ds)
   if method.startswith("KTSP") or method=="TSP1":
    nominal=int(method[4:]) if method.startswith("KTSP") else 1
    need(all(int(r["n_tested"])==nominal for r in ordered),
         "kTSP native voting test count missing "+ds)
  full.extend(rows);manifests.append(m)
if q.full470:
 need(len(manifests)==470 and len(full)==5170,"not full original 470 model folds")
 for ds in sets:
  rmax=2 if ds in ("west-2001","yeoh-2002-v1") else 10
  for rep in range(rmax):
   folds=reps[ds,rep]
   need(set(folds)==set(range(5)),f"incomplete repetition {ds}/{rep}")
   taken=set()
   for fold in range(5):
    need(not(taken&folds[fold]),"duplicate E patient within repetition")
    taken |= folds[fold]
   N=next(int(m["n_total"]) for m in manifests if m["dataset"]==ds)
   need(len(taken)==N,"not each patient tested once per repetition")
out=Path(q.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all_fold_method_scores.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=full[0].keys());w.writeheader();w.writerows(full)
with (out/"all_fold_manifests.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=manifests[0].keys());w.writeheader();w.writerows(manifests)
dataset_means=[]
for ds in sets:
 for method in arms:
  z=[r for r in full if r["dataset"]==ds and r["method"]==method]
  d=dict(dataset=ds,method=method,n_folds=len(z))
  for metric in ("positive_f1","macro_f1","ba","auc","logloss","brier"):
   d[metric]=sum(float(x[metric]) for x in z)/len(z)
  for metric in ("mean_active","mean_tested"):
   good=[float(x[metric]) for x in z if x[metric] not in ("NA","")]
   d[metric]=sum(good)/len(good) if good else ""
  dataset_means.append(d)
with (out/"eleven_dataset_method_means.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=dataset_means[0].keys())
 w.writeheader();w.writerows(dataset_means)
print("SAVS470_ALL_11_ORIGINAL_COHORTS_SOURCE_RECONSTRUCTION=PASS")
print("SAVS470_"+str(len(full))+"_METHOD_FOLDS_PATIENT_METRICS_REPLAYED=PASS")
print("SAVS470_PERSONWISE_ABSTENTION_MATCHES_RULE_COUNTS=PASS")
if q.full470: print("SAVS470_ORIGINAL_470_REPETITION_PATIENT_COVERAGE=PASS")
for method in arms:
 z=[r for r in dataset_means if r["method"]==method]
 print("SAVS470_EQUAL_11_DATASET_MEAN",method,
  " ".join(k+"="+format(sum(r[k] for r in z)/11,".6f")
           for k in ("positive_f1","macro_f1","auc","logloss"))+
  (" mean_active="+format(sum(r["mean_active"] for r in z)/11,".6f")
    if all(isinstance(r["mean_active"],float) for r in z) else ""))
print("SAVS470_SOURCE_ORIGINAL_PROCESSED_ARFF_PARITY_NOT_VERIFIED")
