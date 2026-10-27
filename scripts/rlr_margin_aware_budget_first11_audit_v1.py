#!/usr/bin/env python3
"""Independent patient-probability replay of 11 original paper E folds.

Do not treat these previously explored public Schliep source reconstructions
as SHA-verified original ARFF or untouched external confirmation.
"""
import argparse,csv,math,collections
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True)
p.add_argument("--out",required=True)
a=p.parse_args()
datasets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
"chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
"shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
arms=("RLR20","MARGIN_RLR20","RLR50","MARGIN_RLR50",
"KTSP10","KTSP15","KTSP20","KTSP25","KTSP_M20_GLOBAL",
"KTSP_M20_LOCAL","KTSP_M50_GLOBAL","KTSP_M50_LOCAL")
def check(cond,msg):
 if not cond: raise SystemExit("MB_RLR_AUDIT_STOP "+msg)
def rd(path):
 with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def stats(y,p,h):
 def f1(cl):
  tp=sum(x==cl and z==cl for x,z in zip(y,h))
  fp=sum(x!=cl and z==cl for x,z in zip(y,h))
  fn=sum(x==cl and z!=cl for x,z in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 n=len(y);n1=y.count(1);n0=y.count(0)
 check(min(n0,n1)>0,"invalid E class distribution")
 ix=sorted(range(n),key=lambda i:p[i])
 ranks=[0.]*n;i=0
 while i<n:
  j=i+1
  while j<n and p[ix[j]]==p[ix[i]]:j+=1
  mid=(i+1+j)/2.
  for k in range(i,j):ranks[ix[k]]=mid
  i=j
 cp=[min(1-1e-15,max(1e-15,z)) for z in p]
 return {"positive_f1":f1(1),"macro_f1":(f1(1)+f1(0))/2,
     "ba":(sum(yi==0 and hi==0 for yi,hi in zip(y,h))/n0+
          sum(yi==1 and hi==1 for yi,hi in zip(y,h))/n1)/2,
     "auc":(sum(ranks[i] for i,yi in enumerate(y) if yi==1)-
            n1*(n1+1)/2)/(n1*n0),
     "logloss":sum(-(yi*math.log(pi)+(1-yi)*math.log1p(-pi))
                   for yi,pi in zip(y,cp))/n,
     "brier":sum((yi-pi)**2 for yi,pi in zip(y,p))/n}
def budget(v):
 return max(2,min(50,int(math.floor(v+0.5))))
allrows=[]
for ds in datasets:
 paths=list(Path(a.input).rglob(ds+"/method_metrics.csv"))
 check(len(paths)==1,ds+" needs exactly one cohort archive")
 root=paths[0].parent
 m=rd(paths[0]);probs=rd(root/"patient_E.csv")
 manifest=rd(root/"train_budget_manifest.csv")
 check(len(m)==len(arms) and {x["method"] for x in m}==set(arms),
      ds+" arm count wrong")
 check(len(manifest)==1,ds+" one train manifest")
 z=manifest[0]
 check(z["dataset"]==ds and z["split_id"]=="rep00_fold00" and
       z["provenance"]=="K_MATCH_ALWAYS_CHOSEN_FROM_ORIGINAL_RLR_OUTER_T_ONLY",
       ds+" wrong matching protocol")
 T=z["T_ids"].split("|");E=z["E_ids"].split("|")
 check(len(T)==int(z["n_T"]) and len(E)==int(z["n_E"]) and
       len(set(T))==len(T) and len(set(E))==len(E) and
       not(set(T)&set(E)),ds+" training/test overlap")
 ranges={"KTSP10":10,"KTSP15":15,"KTSP20":20,"KTSP25":25,
  "KTSP_M20_GLOBAL":budget(float(z["K20_global"])),
  "KTSP_M20_LOCAL":budget(float(z["K20_train_patient"])),
  "KTSP_M50_GLOBAL":budget(float(z["K50_global"])),
  "KTSP_M50_LOCAL":budget(float(z["K50_train_patient"]))}
 group=collections.defaultdict(list)
 for row in probs:group[row["method"]].append(row)
 check(set(group)==set(arms),ds+" missing patient model")
 for model in arms:
  meta=next(x for x in m if x["method"]==model)
  rr=group[model]
  check(len(rr)==len(E) and {x["sample_id"] for x in rr}==set(E),
        ds+"/"+model+" missing E patient")
  rr=sorted(rr,key=lambda x:x["sample_id"])
  y=[int(x["truth"]) for x in rr]
  pp=[float(x["prob"]) for x in rr]
  h=[int(x["label"]) for x in rr]
  activity=[int(x["active"]) for x in rr]
  cost=[int(x["tested"]) for x in rr]
  check(all(x in (0,1) for x in y+h) and
        all(math.isfinite(x) and 0<=x<=1 for x in pp),
        ds+"/"+model+" invalid class probability")
  check(all(0<=v<=u<=50 for v,u in zip(activity,cost)),
        ds+"/"+model+" wrong activity")
  if model.startswith("RLR") or model.startswith("MARGIN_RLR"):
   check(all(hi==int(pi>=.5) for hi,pi in zip(h,pp)),
         ds+" RLR threshold not standard")
  else:
   check(1<=int(meta["nominal_pairs"])<=ranges[model] and
         len(set(cost))==1 and cost[0]==int(meta["nominal_pairs"]) and
         all(x==cost[0] for x in activity),
         ds+" real kTSP realized K mismatch "+model)
  actual=stats(y,pp,h)
  for key,val in actual.items():
   check(abs(float(meta[key])-val)<1e-9,
         f"{ds}/{model}/{key} mismatch")
  for key,vals in (("test_mean_active_pairs",activity),
                   ("test_mean_tested_pairs",cost)):
   check(abs(float(meta[key])-sum(vals)/len(vals))<1e-9,
         f"{ds}/{model}/{key} mismatch")
  allrows.append(meta)
dest=Path(a.out);dest.mkdir(parents=True,exist_ok=True)
with (dest/"all11_firstfold_matched_budget_results.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=allrows[0].keys());w.writeheader();w.writerows(allrows)
print("MB_RLR_ALL_11_DATASETS_144_METHOD_ROWS_REPLAYED=PASS" if len(allrows)==144 else "MB_RLR_AUDIT_ERROR_COUNT")
check(len(allrows)==11*len(arms),"not 11 x12 arms")
print("MB_RLR_REAL_KTSP_MATCH_K_SELECTED_FROM_ORIGINAL_RLR_T_ONLY=PASS")
print("MB_RLR_ORIGINAL_E_PATIENTS_ALIGNED_AND_SCORES_INDEPENDENT=PASS")
for model in arms:
 rows=[x for x in allrows if x["method"]==model]
 def mean(key):return sum(float(x[key]) for x in rows)/len(rows)
 print(f"MB_RLR_MEAN {model} F1={mean('positive_f1'):.6f} "+
       f"MacroF1={mean('macro_f1'):.6f} AUC={mean('auc'):.6f} "+
       f"LL={mean('logloss'):.6f} n_global={mean('global_nonzero_pairs'):.3f} "+
       f"n_patient={mean('test_mean_active_pairs'):.3f}")
print("MB_RLR_ORIGINAL_PROCESSED_ARFF_PARITY_NOT_VERIFIED")
