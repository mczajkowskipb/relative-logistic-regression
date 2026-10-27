#!/usr/bin/env python3
"""SAVS3 independent patient replay and proper honest-S/H gate audit."""
import argparse,csv,math,collections
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True)
p.add_argument("--out",required=True)
arg=p.parse_args()
sets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
modes=("RLR10","RELATIVE10PCT")
arms=("HONEST_LOOKAHEAD_V_GATED","V1_ROOT_FIRST_TREE",
      "V1_GLOBAL_SAME_RULES","GLOBAL_FOUR_RULES")
def req(ok,msg):
 if not ok:raise SystemExit("SAVS3_REPLAY_FAIL: "+msg)
def read(path):
 with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def metrics(rows):
 truth=[int(x["truth"]) for x in rows]
 p=[float(x["prob"]) for x in rows]
 labels=[int(x["label"]) for x in rows]
 active=[int(x["n_active"]) for x in rows]
 considered=[int(x["n_considered"]) for x in rows]
 roots=[int(x["n_root"]) for x in rows]
 n=len(rows)
 req(n>0 and min(truth.count(0),truth.count(1))>0,"one class")
 req(all(0<=x<=k<=8 and r in (0,3,5) for x,k,r in
     zip(active,considered,roots)),"invalid activity")
 req(all(math.isfinite(x) and 0<x<1 for x in p),"invalid p")
 req(all(l==int(x>=.5) for l,x in zip(labels,p)),"label threshold")
 def f1(cl):
  tp=sum(a==cl and b==cl for a,b in zip(truth,labels))
  fp=sum(a!=cl and b==cl for a,b in zip(truth,labels))
  fn=sum(a==cl and b!=cl for a,b in zip(truth,labels))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 pp=[min(1-1e-15,max(1e-15,x)) for x in p]
 z={
   "positive_f1":f1(1),
   "macro_f1":(f1(0)+f1(1))/2,
   "ba":(sum(a==0 and b==0 for a,b in zip(truth,labels))/truth.count(0)+
         sum(a==1 and b==1 for a,b in zip(truth,labels))/truth.count(1))/2,
   "logloss":sum(-(a*math.log(b)+(1-a)*math.log1p(-b))
                  for a,b in zip(truth,pp))/n,
   "mean_active":sum(active)/n,
   "median_active":sorted(active)[n//2] if n%2 else
     (sorted(active)[n//2-1]+sorted(active)[n//2])/2,
   "mean_comparisons":sum(a+b for a,b in zip(considered,roots))/n,
   "pct_zero_active":sum(a==0 for a in active)/n,
   "pct_le2_active":sum(a<=2 for a in active)/n
 }
 return z
summary=[]
for ds in sets:
 paths=list(Path(arg.input).rglob(ds+"/metrics_V.csv"))
 req(len(paths)==1,f"one metrics file expected {ds}")
 root=paths[0].parent
 mm=read(root/"metrics_V.csv")
 pp=read(root/"patient_V.csv")
 ss=read(root/"honest_selection_D.csv")
 req(len(mm)==8 and len(ss)==2,"wrong method counts "+ds)
 req(set((x["margin_mode"],x["method"]) for x in mm)==
      {(z,a) for z in modes for a in arms},"missing expected variants "+ds)
 for sel in ss:
  req(int(sel["n_S"])+int(sel["n_H"])==int(sel["n_D"]) and
      int(sel["n_D"])+int(sel["n_V"])==int(sel["n_T"]),
      "overlap inventory "+ds)
  req(sel["evidence"]=="S_ONLY_ROOT_POOL_H_SELECTION_D_FINALFIT_V_TEST_NO_ORIGINAL_E",
      "not honest original outer E condition "+ds)
  if int(sel["honest_tree_accepted"]):
   req(float(sel["H_tree_macroF1"])>=float(sel["H_global_macroF1"])+.01-1e-10 and
       float(sel["H_tree_BA"])>=float(sel["H_global_BA"])-.02-1e-10,
       "incorrect acceptance "+ds)
 for z in mm:
  matching=[x for x in pp if
    x["margin_mode"]==z["margin_mode"] and x["method"]==z["method"]]
  req(len(matching)==int(z["n_V"]) and
      len({x["sample_id"] for x in matching})==len(matching),
      "patient record mismatch "+ds)
  rr=sorted(matching,key=lambda x:x["sample_id"])
  observed=metrics(rr)
  for name,v in observed.items():
   req(abs(float(z[name])-v)<1e-9,
      f"replay mismatch {ds}/{z['margin_mode']}/{z['method']}/{name}")
  summary.append(z)
out=Path(arg.out);out.mkdir(parents=True,exist_ok=True)
with (out/"savs3_honest11_method_V_metrics.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
print("SAVS3_HONEST_EIGHTY_EIGHT_METHODS_ALL_PATIENTS_REPLAYED=PASS")
print("SAVS3_S_POOL_H_SELECTION_AND_FINAL_V_PARTITION=PASS")
for mode in modes:
 for arm in arms:
  v=[x for x in summary if x["margin_mode"]==mode and x["method"]==arm]
  print("SAVS3_DATASET_MEANS",mode,arm," ".join(
   key+"="+format(sum(float(x[key]) for x in v)/11,".6f")
   for key in ("positive_f1","macro_f1","ba","logloss","mean_active","mean_comparisons")))
print("SAVS3_NOT_ORIGINAL_E_OR_UNTOUCHED_CONFIRMATION")
