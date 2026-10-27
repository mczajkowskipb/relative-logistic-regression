#!/usr/bin/env python3
"""Independent 11 cohort q-margin RLR inner V patient-level score replay."""
import csv,math,argparse,collections,statistics
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True)
p.add_argument("--out",required=True)
q=p.parse_args()
datasets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
methods=tuple(f"QRLR_K{k:02d}_Q{q:02d}" for k in (10,20) for q in (0,25,50))
def need(value,message):
 if not value:raise SystemExit("QRLR_DV_AUDIT_FAIL "+message)
def read(path):
 with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def stats(y,h,p):
 def f1(cls):
  tp=sum(a==cls and b==cls for a,b in zip(y,h))
  fp=sum(a!=cls and b==cls for a,b in zip(y,h))
  fn=sum(a==cls and b!=cls for a,b in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 n1=y.count(1);n0=y.count(0);need(n1>0 and n0>0,"one V class")
 order=sorted(range(len(p)),key=lambda i:p[i])
 ranks=[0.]*len(p);i=0
 while i<len(p):
  j=i+1
  while j<len(p) and p[order[j]]==p[order[i]]:j+=1
  rank=(i+1+j)/2
  for k in range(i,j):ranks[order[k]]=rank
  i=j
 clip=[max(1e-15,min(1-1e-15,x)) for x in p]
 return {"positive_f1":f1(1),"macro_f1":(f1(0)+f1(1))/2,
  "ba":(sum(a==0 and b==0 for a,b in zip(y,h))/n0+
        sum(a==1 and b==1 for a,b in zip(y,h))/n1)/2,
  "auc":(sum(ranks[i] for i,a in enumerate(y) if a==1)-
         n1*(n1+1)/2)/(n1*n0),
  "logloss":sum(-(a*math.log(b)+(1-a)*math.log1p(-b))
                 for a,b in zip(y,clip))/len(y),
  "brier":sum((a-b)**2 for a,b in zip(y,p))/len(y)}
outrows=[]
for ds in datasets:
 paths=list(Path(q.input).rglob(ds+"/metrics_V.csv"))
 need(len(paths)==1,ds+" missing/duplicate metrics")
 folder=paths[0].parent
 ms=read(folder/"metrics_V.csv")
 ps=read(folder/"patient_V.csv")
 frozen=read(folder/"frozen_D_margin_panels.csv")
 man=read(folder/"manifest.csv")
 need(len(ms)==len(frozen)==6 and len(man)==1,ds+" wrong 6 panel count")
 need({x["method"] for x in ms}==set(methods),ds+" methods mismatch")
 D=set(man[0]["D_ids"].split("|"));V=set(man[0]["V_ids"].split("|"))
 need(len(D)==int(man[0]["n_D"]) and len(V)==int(man[0]["n_V"])
       and not(D & V),ds+" train/validation contamination")
 for model in methods:
  m=next(x for x in ms if x["method"]==model)
  rr=[x for x in ps if x["method"]==model]
  need({x["sample_id"] for x in rr}==V and len(rr)==len(V),
       ds+"/"+model+" patient V coverage")
  y=[int(x["truth"]) for x in rr]
  h=[int(x["label"]) for x in rr]
  probs=[float(x["prob"]) for x in rr]
  activities=[int(x["n_active"]) for x in rr]
  tested=[int(x["n_tested"]) for x in rr]
  need(all(a in (0,1) and b in (0,1) for a,b in zip(y,h)) and
        all(0<z<1 and math.isfinite(z) for z in probs) and
        all(l==int(p>=.5) for l,p in zip(h,probs)),
       ds+"/"+model+" labels/probabilities invalid")
  need(all(0<=v<=u<=int(m["nominal_K"]) for v,u in
           zip(activities,tested)) and
       len(set(tested))==1 and tested[0]==int(m["n_tested"]),
       ds+"/"+model+" invalid active ternary profile")
  vals=stats(y,h,probs)
  for key,val in vals.items():
   need(abs(float(m[key])-val)<1e-9,
        ds+"/"+model+"/"+key+" differs from patient metrics")
  need(abs(float(m["mean_active_V"])-statistics.mean(activities))<1e-9,
      ds+"/"+model+" mean active mismatch")
  f=next(x for x in frozen if x["method"]==model)
  margins=list(map(float,f["margins"].split("|")))
  need(len(margins)==int(m["nominal_K"]) and min(margins)>=10,
      ds+"/"+model+" thresholds incorrect")
  need((all(x==10 for x in margins) if int(model[-2:])==0 else True),
       ds+"/"+model+" q0 original threshold not 10")
  need(m["evidence"]=="D_FIT_UNSUPERVISED_THRESHOLDS_V_TEST_NO_ORIGINAL_OUTER_E",
      ds+"/"+model+" evidence mismatch")
  outrows.append(m)
out=Path(q.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_quantile_rank_DV_metrics.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=outrows[0].keys());w.writeheader();w.writerows(outrows)
print("QRLR_ALL11_SIXTY_SIX_METHOD_PATIENT_METRICS_REPLAYED=PASS")
print("QRLR_NO_OUTER_E_OR_SAMPLE_OVERLAP=PASS")
for name in methods:
 rows=[x for x in outrows if x["method"]==name]
 def mean(key):return statistics.mean(float(x[key]) for x in rows)
 print("QRLR_MEAN",name,
  f"F1={mean('positive_f1'):.6f}",
  f"MacroF1={mean('macro_f1'):.6f}",
  f"AUC={mean('auc'):.6f}",
  f"LL={mean('logloss'):.6f}",
  f"active={mean('mean_active_V'):.3f}",
  f"tested={mean('n_tested'):.3f}")
print("QRLR_DV_EXPLORATORY_PREVIOUSLY_EXPOSED_COHORTS")
