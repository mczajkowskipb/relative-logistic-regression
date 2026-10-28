#!/usr/bin/env python3
"""Independent patientwise replay: strict train-only 11-cohort RLR-distilled trees.

No claim of original TSPDT reproduction, genuine new tree novelty, original
ARFF parity, or independent cohort validation without direct evidence.
"""
import argparse,csv,math,statistics,collections
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--input",required=True)
p.add_argument("--out",required=True)
args=p.parse_args()
sets=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
arms=("LABEL_BINARY_SAME_RLR20","LABEL_BINARY_TSP_TOP20",
      "DISTILL_BINARY_RLR20","DISTILL_TERNARY_RLR20",
      "RLR20_TEACHER","KTSP5_NATIVE")
tree_arms=set(arms[:4])
def need(yes,msg):
 if not yes:raise SystemExit("RLRTREE_INDEPENDENT_AUDIT_STOP "+msg)
def read(path):
 with path.open(newline="",encoding="utf-8") as f:
  return list(csv.DictReader(f))
def calc(y,pr,cl,depth,teacher):
 n=len(y);n0=y.count(0);n1=y.count(1)
 need(min(n0,n1)>0,"validation single class")
 def f1(label):
  tp=sum(yi==label and pred==label for yi,pred in zip(y,cl))
  fp=sum(yi!=label and pred==label for yi,pred in zip(y,cl))
  fn=sum(yi==label and pred!=label for yi,pred in zip(y,cl))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 ordering=sorted(range(n),key=lambda i:pr[i])
 ranks=[0.]*n;i=0
 while i<n:
  j=i+1
  while j<n and pr[ordering[j]]==pr[ordering[i]]:j+=1
  mid=(i+1+j)/2
  for k in range(i,j):ranks[ordering[k]]=mid
  i=j
 cp=[max(1e-15,min(1-1e-15,z)) for z in pr]
 return dict(
  positive_f1=f1(1),macro_f1=(f1(0)+f1(1))/2,
  ba=(sum(yi==0 and hi==0 for yi,hi in zip(y,cl))/n0+
      sum(yi==1 and hi==1 for yi,hi in zip(y,cl))/n1)/2,
  auc=(sum(ranks[i] for i,yi in enumerate(y) if yi==1)-
      n1*(n1+1)/2)/(n0*n1),
  logloss=sum(-(yi*math.log(pi)+(1-yi)*math.log1p(-pi))
              for yi,pi in zip(y,cp))/n,
  brier=sum((yi-pi)**2 for yi,pi in zip(y,pr))/n,
  mean_evaluated_pairs=statistics.mean(depth),
  max_evaluated_pairs=max(depth),
  student_teacher_agreement=sum(int(z>=.5)==int(t>=.5)
                             for z,t in zip(pr,teacher))/n,
  teacher_soft_cross_entropy=sum(
   -t*math.log(x)-(1-t)*math.log1p(-x)
   for t,x in zip(teacher,cp))/n,
  teacher_fidelity_brier=sum((t-x)**2
                            for t,x in zip(teacher,pr))/n)
summ=[]
for ds in sets:
 paths=list(Path(args.input).rglob(ds+"/metrics_V.csv"))
 need(len(paths)==1,f"exactly one metrics table {ds}")
 folder=paths[0].parent
 met=read(folder/"metrics_V.csv")
 patient=read(folder/"patient_V.csv")
 nodes=read(folder/"tree_nodes_D.csv")
 man=read(folder/"manifest.csv")
 need(len(met)==6 and len(man)==1 and
      {z["method"] for z in met}==set(arms),
      f"{ds}: six frozen methods and manifest required")
 m=man[0]
 need(m["source"]=="PUBLIC_SCHLIEP_RECONSTRUCTION_ARFF_PARITY_UNVERIFIED" and
      int(m["teacher_K"])==20 and int(m["teacher_OOF_folds"])==3 and
      int(m["max_depth"])==2 and int(m["max_internal"])==3 and
      int(m["n_D"])+int(m["n_V"])==int(m["n_outer_T"]),
      "original T/D/V provenance misrecorded "+ds)
 D=set(m["D_ids"].split("|"));V=set(m["V_ids"].split("|"))
 need(not(D&V) and len(D)==int(m["n_D"]) and len(V)==int(m["n_V"]),
      "D/V overlap or patient accounting "+ds)
 for method in arms:
  mm=next(z for z in met if z["method"]==method)
  pp=[z for z in patient if z["method"]==method]
  need(len(pp)==len(V) and {z["sample_id"] for z in pp}==V,
       "per-model patient coverage "+ds+"/"+method)
  rr=sorted(pp,key=lambda z:z["sample_id"])
  y=[int(z["truth"]) for z in rr]
  prob=[float(z["prob"]) for z in rr]
  pred=[int(z["label"]) for z in rr]
  depth=[int(z["evaluated_pairs"]) for z in rr]
  teacher=[float(z["teacher_prob"]) for z in rr]
  need(all(yi in (0,1) and hi in (0,1) for yi,hi in zip(y,pred)) and
       all(math.isfinite(z) and 0<z<1 for z in prob+teacher),
       "invalid scores "+ds+"/"+method)
  if method!="KTSP5_NATIVE":
   need(all(hi==int(q>=.5) for hi,q in zip(pred,prob)),
        "student/RLR label threshold inconsistency "+ds+"/"+method)
  if method in tree_arms:
   need(all(0<=d<=2 for d in depth) and
        int(mm["n_internal_nodes"])<=3 and
        int(mm["n_leaves"])<=7 and
        int(mm["n_unique_model_genes"])<=6 and
        int(mm["n_tree_pairs"])<=3,
        "autonomous compact tree size invalid "+ds+"/"+method)
   ns=[z for z in nodes if z["method"]==method]
   need(len(ns)==int(mm["n_internal_nodes"])+int(mm["n_leaves"]),
        "node table size mismatched "+ds+"/"+method)
   need(sum(bool(z["split_pair"]) for z in ns)==int(mm["n_internal_nodes"]),
        "internal nodes mismatched "+ds+"/"+method)
   if int(mm["n_internal_nodes"])>0:
    need(max(depth)>0,"nontrivial trained tree never used")
  elif method=="KTSP5_NATIVE":
   need(len(set(depth))==1 and 1<=depth[0]<=5,
        "native kTSP5 actual K invalid "+ds)
  values=calc(y,prob,pred,depth,teacher)
  for key,v in values.items():
   need(abs(float(mm[key])-v)<1e-9,
        "patientwise metric mismatch "+ds+"/"+method+"/"+key)
  need(mm["evidence"]=="ORIGINAL_OUTER_T_D_ONLY_TRAIN_V_ONLY_COMPARISON_NO_E",
       "incorrect evidential scope")
  summ.append(mm)
out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_distilled_trees_V_metrics.csv").open("w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=summ[0].keys())
 w.writeheader();w.writerows(summ)
print("RLRTREE_ALL11_SIXTY_SIX_PATIENTWISE_V_METHODS_REPLAYED=PASS")
print("RLRTREE_STRICT_D_V_COVERAGE_AND_TREE_PATH_BUDGET_REPLAYED=PASS")
for method in arms:
 z=[a for a in summ if a["method"]==method]
 fields=("positive_f1","macro_f1","auc","logloss",
         "mean_evaluated_pairs","student_teacher_agreement")
 print("RLRTREE_11_COHORT_MEAN",method," ".join(
   field+"="+format(statistics.mean(float(x[field]) for x in z),".6f")
   for field in fields))
print("RLRTREE_DV_PREVIOUSLY_EXPOSED_COHORTS_NOT_OUTER_E_CONFIRMATION")
