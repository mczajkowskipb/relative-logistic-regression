#!/usr/bin/env python3
"""RLR-GLVT v5 independently replay actual inner-V patient scores and complexity."""
import argparse,csv,collections,math,statistics,itertools,random
from pathlib import Path
q=argparse.ArgumentParser()
q.add_argument("--input",required=True)
q.add_argument("--out",required=True)
a=q.parse_args()
cohorts=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
models=("GLOBAL_LABEL_ONLY","GLOBAL_RLR_LOGIT_OOF",
        "GREEDY_LABEL_TSP_TOP20","RLR20_TEACHER","KTSP5_NATIVE")
tree=set(models[:3]);global_models=set(models[:2])
def check(cond,msg):
 if not cond:raise SystemExit("RLRGLVT_AUDIT_FAIL "+msg)
def read(path):
 with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def calc(y,h,p,teacher,checks):
 n=len(y);n0=y.count(0);n1=y.count(1)
 check(n0>0 and n1>0,"invalid class support in validation")
 def f1(cls):
  tp=sum(x==cls and v==cls for x,v in zip(y,h))
  fp=sum(x!=cls and v==cls for x,v in zip(y,h))
  fn=sum(x==cls and v!=cls for x,v in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 order=sorted(range(n),key=lambda i:p[i]);ranks=[0.]*n
 i=0
 while i<n:
  j=i+1
  while j<n and p[order[j]]==p[order[i]]:j+=1
  mid=(i+1+j)/2
  for k in range(i,j):ranks[order[k]]=mid
  i=j
 cc=[max(1e-15,min(1-1e-15,v)) for v in p]
 ct=[max(1e-8,min(1-1e-8,v)) for v in p]
 return dict(
  positive_f1=f1(1),macro_f1=(f1(0)+f1(1))/2,
  ba=.5*(sum(x==0 and z==0 for x,z in zip(y,h))/n0+
         sum(x==1 and z==1 for x,z in zip(y,h))/n1),
  auc=(sum(ranks[i] for i,x in enumerate(y) if x==1)-
       n1*(n1+1)/2)/(n0*n1),
  logloss=sum(-x*math.log(v)-(1-x)*math.log1p(-v)
              for x,v in zip(y,cc))/n,
  brier=sum((x-v)**2 for x,v in zip(y,p))/n,
  mean_evaluated_pairs=statistics.mean(checks),
  max_evaluated_pairs=max(checks),
  teacher_decision_agreement=sum((x>=.5)==(v>=.5)
       for x,v in zip(teacher,p))/n,
  teacher_soft_cross_entropy=sum(-x*math.log(v)-
       (1-x)*math.log1p(-v) for x,v in zip(teacher,ct))/n,
  teacher_fidelity_brier=sum((x-v)**2 for x,v in zip(teacher,p))/n)
allrows=[]
for ds in cohorts:
 paths=list(Path(a.input).rglob(ds+"/metrics_V.csv"))
 check(len(paths)==1,ds+" missing or duplicate cohort result")
 folder=paths[0].parent
 m=read(folder/"manifest.csv")
 metr=read(folder/"metrics_V.csv")
 pred=read(folder/"patient_V.csv")
 nodes=read(folder/"nodes_D.csv")
 check(len(m)==1 and len(metr)==5 and
       {x["method"] for x in metr}==set(models),
       ds+" frozen 5 arms or manifest mismatch")
 mm=m[0]
 check(mm["source"]=="SCHLIEP_RECONSTRUCTED_ORIGINAL_ARFF_SHA_UNVERIFIED" and
       int(mm["n_outer_T"])==int(mm["n_D"])+int(mm["n_V"]) and
       int(mm["RLR20_K"])==20 and
       int(mm["teacher_oof_folds"])==3 and
       0<=int(mm["teacher_ridge_fallback_folds"])<=3 and
       abs(float(mm["teacher_weight_fixed"])-.25)<1e-9 and
       int(mm["max_depth"])==2 and int(mm["max_internal"])==3,
       ds+" honest D/V provenance or fixed specification invalid")
 D=set(mm["D_ids"].split("|"));V=set(mm["V_ids"].split("|"))
 check(len(D)==int(mm["n_D"]) and len(V)==int(mm["n_V"]) and
       not(D&V) and int(mm["n_outer_E_UNREAD"])>0,
       ds+" D/V patient overlap/outer E")
 for model in models:
  z=next(x for x in metr if x["method"]==model)
  patient=[x for x in pred if x["method"]==model]
  check(len(patient)==len(V) and
        {x["sample_id"] for x in patient}==V,
        ds+"/"+model+" V sample coverage")
  patient=sorted(patient,key=lambda x:x["sample_id"])
  yy=[int(x["truth"]) for x in patient]
  hh=[int(x["label"]) for x in patient]
  pp=[float(x["prob"]) for x in patient]
  tt=[float(x["teacher_prob"]) for x in patient]
  c=[int(x["n_evaluated"]) for x in patient]
  check(all(x in (0,1) for x in yy+hh) and
        all(math.isfinite(x) and 0<=x<=1 for x in pp) and
        all(math.isfinite(x) and 0<x<1 for x in tt),
        ds+"/"+model+" bad patient probabilities")
  if model!="KTSP5_NATIVE":
   check(all(x==int(p>=.5) for x,p in zip(hh,pp)),
        ds+"/"+model+" train-only threshold != reported native label")
  if model in global_models:
   check(all(x in (0,1,2) for x in c) and
         int(z["n_internal_nodes"])<=3 and
         int(z["n_leaves"])<=4 and
         int(z["n_unique_relations"])<=3 and
         int(z["n_unique_model_genes"])<=6,
         ds+"/"+model+" exact global tree exceeds complexity")
  if model=="GREEDY_LABEL_TSP_TOP20":
   check(all(0<=x<=2 for x in c) and
         int(z["n_internal_nodes"])<=3 and
         int(z["n_unique_relations"])<=3,
         ds+" label tree max node count invalid")
  if model=="KTSP5_NATIVE":
   check(len(set(c))==1 and 1<=c[0]<=5,
         ds+" native switchBox realized kTSP K invalid")
  if model in tree:
   nd=[x for x in nodes if x["method"]==model]
   n_internal=int(z["n_internal_nodes"])
   n_leaf=int(z["n_leaves"])
   check(len(nd)==n_internal+n_leaf and
         sum(bool(x["node_predicate"]) for x in nd)==n_internal,
         ds+"/"+model+" tree structure mismatch")
   if model in global_models:
    check(all(">" in x["node_predicate"]
              for x in nd if x["node_predicate"]),
          ds+"/"+model+" missing genuine relational node definition")
  vals=calc(yy,hh,pp,tt,c)
  for key,value in vals.items():
   check(abs(float(z[key])-value)<1e-9,
         ds+"/"+model+"/"+key+" patient replay mismatch")
  check(z["training_scope"]=="ORIGINAL_OUTER_T_TO_D_V_NO_E",
        ds+"/"+model+" invalid scope")
  allrows.append(z)
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_global_logit_variance_tree_DV_metrics.csv").open(
     "w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=allrows[0].keys())
 w.writeheader();w.writerows(allrows)
print("RLRGLVT_ALL11_FIFTY_FIVE_PATIENT_SCORES_INDEPENDENTLY_REPLAYED=PASS")
print("RLRGLVT_ALL11_REAL_MAX_TWO_PAIR_PATH_BUDGET_AND_DV=PASS")
def mean(model,key):
 return statistics.mean(float(z[key]) for z in allrows if z["method"]==model)
for model in models:
 print("RLRGLVT_11_COHORT_MEAN",model,
  "F1=%.6f"%mean(model,"positive_f1"),
  "MacroF1=%.6f"%mean(model,"macro_f1"),
  "AUC=%.6f"%mean(model,"auc"),
  "LL=%.6f"%mean(model,"logloss"),
  "NPAIR=%.4f"%mean(model,"mean_evaluated_pairs"))
# No falsely independent 55 model cohorts. EXACT 11 biological cohort test
# units, all predeclared 2^11 signflips pre-fixed for 4 contrasts.
by={(z["dataset"],z["method"]):z for z in allrows}
contrasts=[
 ("GLOBAL_RLR_LOGIT_OOF","GLOBAL_LABEL_ONLY"),
 ("GLOBAL_RLR_LOGIT_OOF","GREEDY_LABEL_TSP_TOP20"),
 ("GLOBAL_LABEL_ONLY","GREEDY_LABEL_TSP_TOP20"),
 ("GLOBAL_RLR_LOGIT_OOF","RLR20_TEACHER"),
 ("GLOBAL_RLR_LOGIT_OOF","KTSP5_NATIVE")]
tests=[]
rng=random.Random(20261029)
for A,B in contrasts:
 for metric in ("positive_f1","macro_f1","auc","logloss","mean_evaluated_pairs"):
  dif=[float(by[d,A][metric])-float(by[d,B][metric]) for d in cohorts]
  obs=statistics.mean(dif)
  cnt=sum(abs(sum(v*t for v,t in zip(dif,signs))/11)>=abs(obs)-1e-12
    for signs in itertools.product((-1,1),repeat=11))
  boots=sorted(statistics.mean(rng.choices(dif,k=11))
              for _ in range(8000))
  tests.append(dict(A=A,B=B,metric=metric,delta_mean=obs,
    exact_p=cnt/2048,bootstrap95_lo=boots[200],
    bootstrap95_hi=boots[7799]))
with (out/"paired_11dataset_exact_signflips.csv").open(
     "w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=tests[0].keys())
 w.writeheader();w.writerows(tests)
print("RLRGLVT_11_COHORT_PAIRED_EXACT_SIGNFLIP=PASS")
print("RLRGLVT_EVIDENCE_EXPLORATORY_SOURCE_ARFF_PARITY_UNVERIFIED")
