#!/usr/bin/env python3
"""Independent RXA multi-test cascade 11-cohort patient-level metric and routing audit."""
import argparse,collections,csv,glob,itertools,math,pathlib,statistics
ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True)
ap.add_argument("--out",required=True)
args=ap.parse_args()
out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
arms=("ORIGINAL_RLR20","BINARY_MT2_LEAF","TERNARY_MT2_LEAF",
      "BINARY_MT2_RLR_BACKSTOP","TERNARY_MT2_RLR_BACKSTOP",
      "TERNARY_SAME6_ENET")
def read(path):
 with open(path,newline="",encoding="utf8") as fh:return list(csv.DictReader(fh))
def save(name,rows):
 assert rows
 with (out/name).open("w",newline="",encoding="utf8") as fh:
  w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def close(a,b):
 return math.isclose(float(a),float(b),rel_tol=2e-8,abs_tol=2e-9)
def auc(y,p):
 n1=sum(y);n0=len(y)-n1
 assert n1>0 and n0>0
 order=sorted(range(len(y)),key=lambda i:p[i])
 ranks=[0.0]*len(y);j=0
 while j<len(y):
  k=j
  while k+1<len(y) and p[order[k+1]]==p[order[j]]:k+=1
  mean_rank=(j+k+2)/2
  for m in range(j,k+1):ranks[order[m]]=mean_rank
  j=k+1
 return (sum(r for r,t in zip(ranks,y) if t==1)-n1*(n1+1)/2)/(n1*n0)
def measure(y,p):
 q=[min(1-1e-15,max(1e-15,v)) for v in p]
 h=[int(v>=.5) for v in p]
 tp=sum(t==1 and z==1 for t,z in zip(y,h))
 fp=sum(t==0 and z==1 for t,z in zip(y,h))
 fn=sum(t==1 and z==0 for t,z in zip(y,h))
 bins=[min(9,int(v*10)) for v in p]
 ece=sum(len(ix)/len(y)*abs(statistics.mean(p[i] for i in ix)-
                                  statistics.mean(y[i] for i in ix))
         for b in range(10) if (ix:=[i for i,k in enumerate(bins) if k==b]))
 return dict(logloss=statistics.mean(-t*math.log(v)-(1-t)*math.log1p(-v)
                    for t,v in zip(y,q)),
   brier=statistics.mean((t-v)**2 for t,v in zip(y,p)),
   auc=auc(y,p),positive_f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0,
   ece=ece)
files=sorted(glob.glob(args.input+"/**/fold_metrics.csv",recursive=True))
assert len(files)==11,("expected 11 datasets, not 55 fake independent studies",len(files))
allmetrics=[];allpatients=[];allfolds=[];studies=[];by_study={}
for mfile in files:
 directory=pathlib.Path(mfile).parent
 metrics=read(mfile);pts=read(directory/"patients.csv");folds=read(directory/"fold_manifest.csv")
 assert len(metrics)==30 and len(folds)==5
 ds={q["dataset"] for q in metrics}
 assert len(ds)==1
 ds=next(iter(ds))
 assert set(x["method"] for x in metrics)==set(arms)
 assert {int(x["fold"]) for x in metrics}==set(range(5))
 assert all(x["source"]=="SCHLIEP_RECONSTRUCTED_NOT_ORIGINAL_ARFF_SHA_VERIFIED"
            for x in folds),"V1 provenance must not claim original verified"
 for f in folds:
  tr=f["train_ids"].split("|")
  te=f["test_ids"].split("|")
  assert len(tr)==int(f["n_train"]) and len(te)==int(f["n_test"])
  assert set(tr).isdisjoint(set(te))
  assert len(f["pair_6"].split("|"))==6
  assert len(f["orientation"].split("|"))==6
  assert int(f["min_leaf_support"])==5 and int(f["fixed_prior_weight"])==10
 for arm in arms:
  ids=[x["sample_id"] for x in pts if x["method"]==arm]
  assert len(ids)==len(set(ids)),(ds,arm,"repeated outer E patient")
 assert len(pts)==sum(int(z["n_test"]) for z in metrics)
 counts=collections.defaultdict(lambda:collections.defaultdict(dict))
 for x in pts:
  key=(x["sample_id"],int(x["fold"]))
  assert x["method"] not in counts[key]["row"]
  counts[key]["row"][x["method"]]=x
 assert all(set(x["row"])==set(arms) for x in counts.values())
 for f in folds:
  fold=int(f["fold"])
  fold_rows={z["sample_id"] for z in pts
              if z["method"]=="ORIGINAL_RLR20" and int(z["fold"])==fold}
  assert fold_rows==set(f["test_ids"].split("|"))
 for m in metrics:
  fold=int(m["fold"]);arm=m["method"]
  patients=[x for x in pts if int(x["fold"])==fold and x["method"]==arm]
  assert len(patients)==int(m["n_test"])
  y=[int(z["truth"]) for z in patients]
  p=[float(z["prob"]) for z in patients]
  cost=[int(z["pair_checks"]) for z in patients]
  fb=[int(z["backstop"]) for z in patients]
  route=[z["exit_route"] for z in patients]
  assert set(y)=={0,1} and all(0<=z<=1 and math.isfinite(z) for z in p)
  assert all(int(z["pred"])==int(q>=.5) for z,q in zip(patients,p))
  assert all(z in ("NODE1","NODE2","DEFER","RLR","ENET") for z in route)
  measure_0=measure(y,p)
  for name,val in measure_0.items():
   assert close(val,m[name]),(ds,fold,arm,name,val,m[name])
  assert close(statistics.mean(cost),m["mean_pair_checks"])
  assert close(statistics.mean(fb),m["mean_backstop_rate"])
  for label,field in (("NODE1","node1_exit_rate"),("NODE2","node2_exit_rate"),
                      ("DEFER","unresolved_rate")):
   assert close(sum(z==label for z in route)/len(route),m[field])
  nparent=int(m["n_pairs_parent"])
  for z in patients:
   route=z["exit_route"];k=int(z["pair_checks"]);backstop=int(z["backstop"])
   if arm=="ORIGINAL_RLR20":assert route=="RLR" and k==nparent and backstop==0
   elif arm=="TERNARY_SAME6_ENET":assert route=="ENET" and k==6 and backstop==0
   else:
    assert route in ("NODE1","NODE2","DEFER")
    assert k==(3 if route=="NODE1" else (6+nparent if route=="DEFER" and "BACKSTOP" in arm else 6))
    assert backstop==int(route=="DEFER" and "BACKSTOP" in arm)
 for key,pack in counts.items():
  x=pack["row"]
  rlr=float(x["ORIGINAL_RLR20"]["prob"])
  for prefix in ("BINARY","TERNARY"):
   leaf=x[prefix+"_MT2_LEAF"]
   rescue=x[prefix+"_MT2_RLR_BACKSTOP"]
   assert leaf["exit_route"]==rescue["exit_route"],(ds,key,prefix,"route mismatch")
   if rescue["exit_route"]=="DEFER":
    assert close(rescue["prob"],rlr),("non-global RLR logistic rescue",ds,key)
   else:
    assert close(rescue["prob"],leaf["prob"]),("confident leaf changed",ds,key)
 cohort=[]
 for method in arms:
  t=[z for z in metrics if z["method"]==method]
  assert len(t)==5
  cohort.append(dict(dataset=ds,method=method,**{
   field:statistics.mean(float(z[field]) for z in t)
   for field in ("logloss","brier","auc","positive_f1","ece","mean_pair_checks",
                 "mean_backstop_rate","node1_exit_rate","node2_exit_rate","unresolved_rate")
  }))
 studies+=cohort
 by_study[ds]={z["method"]:z for z in cohort}
 allmetrics+=metrics;allpatients+=pts;allfolds+=folds
assert len(allmetrics)==330 and len(allfolds)==55 and len(studies)==66
save("all11_330_fold_metrics.csv",allmetrics)
save("all11_55_fold_manifests.csv",allfolds)
save("all11_all_patient_probabilities.csv",allpatients)
save("all11_study_equal_metrics.csv",studies)
summary=[]
for arm in arms:
 q=[by_study[ds][arm] for ds in sorted(by_study)]
 summary.append(dict(method=arm,n_studies=11,**{
  field:statistics.mean(float(z[field]) for z in q)
  for field in ("logloss","brier","auc","positive_f1","ece",
                "mean_pair_checks","mean_backstop_rate",
                "node1_exit_rate","node2_exit_rate","unresolved_rate")
 }))
save("summary_all_methods.csv",summary)
contrasts=[
 ("TERNARY_MT2_LEAF","BINARY_MT2_LEAF","logloss"),
 ("TERNARY_MT2_LEAF","BINARY_MT2_LEAF","positive_f1"),
 ("TERNARY_MT2_RLR_BACKSTOP","BINARY_MT2_RLR_BACKSTOP","logloss"),
 ("TERNARY_MT2_RLR_BACKSTOP","BINARY_MT2_RLR_BACKSTOP","positive_f1"),
]
inference=[]
for a,b,metric in contrasts:
 delta=[float(by_study[ds][a][metric])-float(by_study[ds][b][metric])
        for ds in sorted(by_study)]
 avg=statistics.mean(delta)
 vals=[abs(sum(s*d for s,d in zip(signs,delta))/11)
       for signs in itertools.product((-1,1),repeat=11)]
 p=sum(z>=abs(avg)-1e-13 for z in vals)/2048
 inference.append(dict(arm_A=a,arm_B=b,metric=metric,
  mean_A_minus_B=avg,exact_two_sided_signflip_p=p,
  n_independent_cohorts=11,source_sha256_parity="NOT_VERIFIED",
  hypothesis_status="EXPLORATORY_11_PREVIOUSLY_EXPOSED_DATASETS"))
# BH across four pre-specified comparisons, including null results.
order=sorted(range(4),key=lambda i:inference[i]["exact_two_sided_signflip_p"])
qvals=[min(1,inference[idx]["exact_two_sided_signflip_p"]*4/(j+1))
       for j,idx in enumerate(order)]
minimum=1.0
for j in reversed(range(4)):
 minimum=min(minimum,qvals[j]);inference[order[j]]["exploratory_BH_q"]=minimum
save("paired_predefined_comparisons.csv",inference)
for z in summary:print("RXA_MTDEF_SUMMARY",z)
for z in inference:print("RXA_MTDEF_PAIRED",z)
print("RXA_MTDEF_V1_INDEPENDENT_330_METHOD_FOLD_METRIC_REPLAY=PASS")
print("RXA_MTDEF_V1_SAME_PATIENT_ROUTING_AND_RLR_BACKSTOP=PASS")
print("RXA_MTDEF_V1_11_COHORT_SIGNFLIP_BH4=PASS")
print("RXA_MTDEF_V1_ORIGINAL_ARFF_SHA256_PARITY=NOT_VERIFIED")
