#!/usr/bin/env python3
"""Independent held-out patient replay for bounded RLR teacher capacity."""
import argparse,csv,math,statistics,itertools,random
from pathlib import Path
ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True)
ap.add_argument("--out",required=True)
args=ap.parse_args()
COHORTS=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
MODELS=tuple(f"{goal}_CAP{k}" for goal in ("TEACHER","LABEL")
              for k in range(4))+("RLR20_REFERENCE",)
def need(ok,msg):
 if not ok:raise SystemExit("RLR_CAPACITY_AUDIT_FAILURE "+msg)
def read(path):
 with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def metric(y,h,p,pt,ref,c):
 n=len(y);n0=y.count(0);n1=y.count(1)
 need(n0>0 and n1>0,"V class support")
 def f1(cl):
  tp=sum(a==cl and b==cl for a,b in zip(y,h))
  fp=sum(a!=cl and b==cl for a,b in zip(y,h))
  fn=sum(a==cl and b!=cl for a,b in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.0
 order=sorted(range(n),key=lambda i:p[i]);rank=[0.]*n
 i=0
 while i<n:
  j=i+1
  while j<n and p[order[j]]==p[order[i]]:j+=1
  for t in range(i,j):rank[order[t]]=(i+1+j)/2
  i=j
 cap=lambda v:max(1e-15,min(1-1e-15,v))
 ll=lambda a,b:-a*math.log(cap(b))-(1-a)*math.log1p(-cap(b))
 return dict(positive_f1=f1(1),macro_f1=(f1(0)+f1(1))/2,
  ba=.5*(sum(a==b==0 for a,b in zip(y,h))/n0+
         sum(a==b==1 for a,b in zip(y,h))/n1),
  auc=(sum(rank[i] for i in range(n) if y[i]==1)-n1*(n1+1)/2)/(n0*n1),
  logloss=statistics.mean(ll(a,b) for a,b in zip(y,p)),
  brier=statistics.mean((a-b)**2 for a,b in zip(y,p)),
  teacher_fidelity_brier=statistics.mean((a-b)**2 for a,b in zip(ref,pt)),
  teacher_soft_cross_entropy=statistics.mean(ll(a,b) for a,b in zip(ref,pt)),
  teacher_agreement=statistics.mean((a>=.5)==(b>=.5) for a,b in zip(ref,pt)),
  mean_checks=statistics.mean(c),max_checks=max(c))
allrows=[]
for ds in COHORTS:
 files=list(Path(args.input).rglob(ds+"/metrics_V.csv"))
 need(len(files)==1,ds+" unique source results")
 folder=files[0].parent
 m=read(folder/"manifest.csv");records=read(files[0]);raw=read(folder/"patient_V.csv")
 need(len(m)==1 and len(records)==9 and
      {x["method"] for x in records}==set(MODELS),ds+" missing 9 methods")
 mm=m[0]
 need(mm["source"]=="SCHLIEP_RECONSTRUCTED_ARFF_SHA_UNVERIFIED" and
      mm["teacher"]=="ORIGINAL_RLR20_D_OOF_3FOLD" and
      mm["candidate_pool_D_only"]=="SWITCHBOX_TOP20" and
      int(mm["n_outer_T"])==int(mm["n_D"])+int(mm["n_V"]) and
      int(mm["n_outer_E_UNREAD"])>0 and
      0<=int(mm["surrogate_OOF_folds"])<=3,ds+" provenance")
 D=set(mm["D_ids"].split("|"));V=set(mm["V_ids"].split("|"))
 need(len(D)==int(mm["n_D"]) and len(V)==int(mm["n_V"]) and not D&V,
      ds+" disjoint D/V")
 baseline={}
 for name in MODELS:
  z=next(x for x in records if x["method"]==name)
  row=sorted((x for x in raw if x["method"]==name),key=lambda x:x["sample_id"])
  need(len(row)==len(V) and {r["sample_id"] for r in row}==V,
       ds+"/"+name+" heldout V coverage")
  y=[int(q["truth"]) for q in row]
  h=[int(q["label"]) for q in row]
  pl=[float(q["prob_label"]) for q in row]
  pt=[float(q["prob_teacher"]) for q in row]
  ref=[float(q["teacherV"]) for q in row]
  c=[int(q["checks"]) for q in row]
  need(all(x in (0,1) for x in y+h) and
       all(math.isfinite(x) and 0<x<1 for x in pl+pt+ref),
       ds+"/"+name+" nonfinite probabilities/classes")
  need(all(int(p>=.5)==v for p,v in zip(pl,h)),ds+"/"+name+" label map")
  if name!="RLR20_REFERENCE":
   k=int(z["cap"]);used=int(z["actual_nodes"])
   need(0<=used<=k<=3 and int(z["n_unique_pairs"])<=used and
        int(z["n_unique_genes"])<=2*used and
        all(0<=q<=min(k,2) for q in c),
        ds+"/"+name+" real tree capacity violation")
   need(all(len([t for t in q["path"].split(" / ") if t])==qchecks
            for q,qchecks in zip(row,c)),
        ds+"/"+name+" path trace/check count mismatch")
   if k==0:need(used==0 and set(c)=={0},ds+" cap0 nonconstant")
   baseline[name]=(float(z["train_optimized_loss"]),
       float(z["train_teacher_sse_fraction"]),
       float(z["train_true_label_ce_fraction"]))
  else:
   need(len(set(c))==1 and c[0]==int(mm["parent_nonzero_beta"]),
        ds+" parent measured comparison count")
  score=metric(y,h,pl,pt,ref,c)
  for key,value in score.items():
   need(abs(float(z[key])-value)<1e-9,
        ds+"/"+name+"/"+key+" patient metric replay mismatch")
  need(z["scope"]=="ORIGINAL_OUTER_FIRST_T_ONLY_D_V_NO_E",
       ds+"/"+name+" invalid training scope")
  allrows.append(z)
 for goal in ("TEACHER","LABEL"):
  loss=[baseline[f"{goal}_CAP{k}"][0] for k in range(4)]
  need(all(loss[i]>=loss[i+1]-1e-8 for i in range(3)),
       ds+"/"+goal+" training oracle frontier NON_MONOTONE")
 for k in range(4):
  teacher=baseline[f"TEACHER_CAP{k}"]
  label=baseline[f"LABEL_CAP{k}"]
  need(teacher[1]<=label[1]+1e-8 and label[2]<=teacher[2]+1e-8,
       ds+f" cap{k} exact objective does NOT dominate comparator on D")
out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_capacity_metrics.csv").open("w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=allrows[0].keys())
 w.writeheader();w.writerows(allrows)
print("RLR_CAPACITY_ALL11_NINETYNINE_PATIENT_METHOD_SUMMARIES_REPLAYED=PASS")
print("RLR_CAPACITY_EXACT_0_TO_3_FRONTIER_DOMINANCE_D=PASS")
print("RLR_CAPACITY_FIXED_DV_NO_OUTER_E_AND_MAX2=PASS")
by={(z["dataset"],z["method"]):z for z in allrows}
def means(key,m):return statistics.mean(float(by[d,m][key]) for d in COHORTS)
for goal in ("TEACHER","LABEL"):
 for k in range(4):
  m=f"{goal}_CAP{k}"
  print("RLR_CAPACITY_11MEAN",m,
        "F1=%.6f"%means("positive_f1",m),
        "AUC=%.6f"%means("auc",m),
        "LL=%.6f"%means("logloss",m),
        "FIDELITY_BRIER=%.6f"%means("teacher_fidelity_brier",m),
        "TRAIN_TEACHER_SSE_RATIO=%.6f"%means("train_teacher_sse_fraction",m),
        "CHECKS=%.4f"%means("mean_checks",m))
tests=[]
rng=random.Random(20261029)
contrasts=[
 (f"TEACHER_CAP{k}",f"LABEL_CAP{k}") for k in range(4)]+[
 ("TEACHER_CAP3","TEACHER_CAP0"),("LABEL_CAP3","LABEL_CAP0")]
for A,B in contrasts:
 for metric_name in ("positive_f1","teacher_fidelity_brier","mean_checks"):
  dif=[float(by[d,A][metric_name])-float(by[d,B][metric_name]) for d in COHORTS]
  obs=statistics.mean(dif)
  p=sum(abs(sum(a*(1 if bits>>j&1 else -1) for j,a in enumerate(dif))/11)>=abs(obs)-1e-12
        for bits in range(2048))/2048
  boots=sorted(statistics.mean(rng.choices(dif,k=11)) for _ in range(4000))
  tests.append(dict(A=A,B=B,metric=metric_name,delta=obs,
      exact_signflip_p=p,bootstrap95_lower=boots[100],bootstrap95_upper=boots[3899]))
with (out/"paired_cohort_signflips.csv").open("w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=tests[0].keys())
 w.writeheader();w.writerows(tests)
print("RLR_CAPACITY_EXACT_ELEVEN_COHORT_SIGNFLIPS=PASS")
print("RLR_CAPACITY_EXPLORATORY_SCHLIEP_RECONSTRUCTION_ORIGINAL_ARFF_UNVERIFIED")
