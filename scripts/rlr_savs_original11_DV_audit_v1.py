#!/usr/bin/env python3
"""Independent patientwise replay of 11 original cohort D/V SAVS arms."""
import argparse,csv,math
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True)
p.add_argument("--out",required=True)
a=p.parse_args()
datasets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
modes=("RLR10","RELATIVE10PCT")
arms=("CONDITIONAL_ANCHOR_ABSTAIN","GLOBAL_SAME_RULES",
      "GLOBAL_FOUR_RULES","ANCHOR_ROOT_ONLY")
def need(x,msg):
 if not x:raise SystemExit("SAVS_DV_AUDIT_STOP "+msg)
def read(path):
 with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def f1(y,h,c):
 tp=sum(z==c and a==c for z,a in zip(y,h))
 fp=sum(z!=c and a==c for z,a in zip(y,h))
 fn=sum(z==c and a!=c for z,a in zip(y,h))
 return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
def recompute(rows):
 y=[int(z["truth"]) for z in rows]
 h=[int(z["label"]) for z in rows]
 p=[float(z["prob"]) for z in rows]
 active=[int(z["n_active"]) for z in rows]
 considered=[int(z["n_considered"]) for z in rows]
 root=[int(z["n_root"]) for z in rows]
 n=len(rows)
 need(min(y.count(0),y.count(1))>0 and all(
      yi in (0,1) and hi in (0,1) and math.isfinite(pi) and 0<=pi<=1
      for yi,hi,pi in zip(y,h,p)),"invalid patient labels/p")
 need(all(0<=ac<=co<=8 and rt in (0,3,5)
          for ac,co,rt in zip(active,considered,root)),
      "rule activity inconsistent")
 need(all(hi==int(pi>=.5) for hi,pi in zip(h,p)),
      "wrong decision threshold")
 ids=sorted(range(n),key=lambda i:p[i])
 ranks=[0.]*n;k=0
 while k<n:
  j=k+1
  while j<n and p[ids[j]]==p[ids[k]]:j+=1
  r=(k+1+j)/2
  for z in range(k,j):ranks[ids[z]]=r
  k=j
 n1=y.count(1);n0=y.count(0)
 x=sorted(active)
 median=x[n//2] if n%2 else (x[n//2-1]+x[n//2])/2
 clipped=[max(1e-15,min(1-1e-15,z)) for z in p]
 return dict(positive_f1=f1(y,h,1),
  macro_f1=(f1(y,h,0)+f1(y,h,1))/2,
  ba=(sum(yy==0 and hh==0 for yy,hh in zip(y,h))/n0+
       sum(yy==1 and hh==1 for yy,hh in zip(y,h))/n1)/2,
  auc=(sum(ranks[i] for i,z in enumerate(y) if z==1)-
        n1*(n1+1)/2)/(n1*n0),
  logloss=sum(-(yy*math.log(pp)+(1-yy)*math.log1p(-pp))
              for yy,pp in zip(y,clipped))/n,
  brier=sum((yy-pp)**2 for yy,pp in zip(y,p))/n,
  mean_pairs_active=sum(active)/n,
  median_pairs_active=median,
  mean_pairs_examined=sum(considered)/n,
  root_comparisons=root[0],
  pct_no_vote=sum(z==0 for z in active)/n,
  pct_le2_active=sum(z<=2 for z in active)/n)
allrows=[]
for ds in datasets:
 found=list(Path(a.input).rglob(ds+"/metrics_V.csv"))
 need(len(found)==1,ds+" missing or multiple metrics")
 metrics=read(found[0]);pred=read(found[0].parent/"patient_V.csv")
 shapes=read(found[0].parent/"structure_D.csv")
 need(len(metrics)==8 and len(shapes)==2,ds+" incomplete arms")
 need(set((x["margin_mode"],x["method"]) for x in metrics)==
      {(z,w) for z in modes for w in arms},ds+" unexpected models")
 for row in metrics:
  rr=[z for z in pred if z["margin_mode"]==row["margin_mode"]
      and z["method"]==row["method"]]
  need(len(rr)==int(row["n_V"]) and
       len({z["sample_id"] for z in rr})==len(rr),
       ds+" misaligned V patient")
  calc=recompute(sorted(rr,key=lambda z:z["sample_id"]))
  for name,value in calc.items():
   need(abs(float(row[name])-value)<1e-9,
        ds+"/"+row["margin_mode"]+"/"+row["method"]+"/"+name)
  need(row["evidence"]==
       "PUBLIC_SOURCE_ORIGINAL_OUTER_T_ONLY_INNER_V_DEVELOPMENT",
       ds+" unjustified E claim")
  allrows.append(row)
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_savs_inner_V_methods.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=allrows[0].keys())
 w.writeheader();w.writerows(allrows)
print("SAVS_ALL11_EIGHTY_EIGHT_DV_METHOD_METRICS_RECALCULATED=PASS")
print("SAVS_ALL11_PATIENTWISE_ABSTENTION_AND_COMPLEXITY_REPLAYED=PASS")
print("SAVS_ORIGINAL_OUTER_E_NEVER_USED=PASS")
for mode in modes:
 for arm in arms:
  v=[z for z in allrows if z["margin_mode"]==mode and z["method"]==arm]
  line=" ".join(
    field+"="+format(sum(float(z[field]) for z in v)/11,".6f")
    for field in ("positive_f1","macro_f1","ba","auc","logloss",
      "mean_pairs_active","mean_pairs_examined","pct_no_vote"))
  print("SAVS_ALL11_DV_MEAN",mode,arm,line)
print("SAVS_DV_IS_EXPLORATORY_NOT_OUTER_E_CONFIRMATION")
