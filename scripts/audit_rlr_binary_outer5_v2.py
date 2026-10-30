#!/usr/bin/env python3
"""Independent replay for 11 cohorts × 5 folds × 3 arms; dataset-level exact sign flips."""
import argparse,collections,csv,glob,itertools,math,pathlib,statistics
arg=argparse.ArgumentParser()
arg.add_argument("--input",required=True)
arg.add_argument("--out",required=True)
args=arg.parse_args();out=pathlib.Path(args.out);out.mkdir(parents=True,exist_ok=True)
methods=("RLR20_ORIGINAL","BINARY_ENET20_PURIST_INSPIRED","BINARY_MCP20_DECAF_INSPIRED")
def read(f):
 with open(f,newline="") as h:return list(csv.DictReader(h))
def dump(name,rows):
 assert rows,(name,"empty")
 with (out/name).open("w",newline="") as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def close(a,b):
 return math.isclose(float(a),float(b),rel_tol=2e-8,abs_tol=2e-9)
def f1(y,h,cls):
 tp=sum(t==cls and z==cls for t,z in zip(y,h))
 fp=sum(t!=cls and z==cls for t,z in zip(y,h))
 fn=sum(t==cls and z!=cls for t,z in zip(y,h))
 return (2*tp/(2*tp+fp+fn)) if 2*tp+fp+fn else 0
def auc(y,prob):
 n1=sum(y);n0=len(y)-n1
 assert n1>0 and n0>0
 ranks=[0.0]*len(prob)
 for k,grp in itertools.groupby(sorted(range(len(prob)),key=lambda i:prob[i]),key=lambda i:prob[i]):
  idx=list(grp);avg=(min([j for j,v in enumerate(sorted(prob)) if v==k])+1+max([j for j,v in enumerate(sorted(prob)) if v==k])+1)/2
  for i in idx:ranks[i]=avg
 return (sum(r for r,t in zip(ranks,y) if t==1)-n1*(n1+1)/2)/(n1*n0)
def ece(y,p):
 bins=[min(9,int(q*10)) for q in p]
 return sum(len(ix)/len(y)*abs(statistics.mean(p[i] for i in ix)-statistics.mean(y[i] for i in ix))
  for b in range(10) if (ix:=[i for i,z in enumerate(bins) if z==b]))
files=sorted(glob.glob(args.input+"/**/fold_metrics.csv",recursive=True))
assert len(files)==11,("expected 11 cohorts",len(files))
allmetrics=[];allpatients=[];dsrows=[]
for metricsfile in files:
 directory=pathlib.Path(metricsfile).parent
 m=read(metricsfile);p=read(directory/"patients.csv");s=read(directory/"fold_manifest.csv")
 assert len(m)==15 and len(s)==5
 ds={z["dataset"] for z in m}
 assert len(ds)==1
 ds=next(iter(ds))
 assert set(z["method"] for z in m)==set(methods)
 assert len(p)==sum(int(z["n_test"]) for z in m)
 # All five folds provide disjoint test patients covering each dataset once per model.
 for name in methods:
  ids=[z["sample_id"] for z in p if z["method"]==name]
  assert len(ids)==len(set(ids)),(ds,name,"test overlap")
 assert len(set(z["sample_id"] for z in p))==len(p)//3
 for z in s:
  assert set(z["train_sample_ids"].split("|")).isdisjoint(set(z["test_sample_ids"].split("|")))
  assert int(z["n_train"])+int(z["n_test"])==len(p)//3
 for z in m:
  pts=[q for q in p if q["method"]==z["method"] and q["fold"]==z["fold"]]
  assert len(pts)==int(z["n_test"])
  y=[int(q["truth"]) for q in pts]
  raw=[float(q["prob"]) for q in pts]
  pred=[int(q["pred"]) for q in pts]
  assert set(y)=={0,1} and all(0<=v<=1 and math.isfinite(v) for v in raw)
  assert all((v>=.5)==bool(z2) for v,z2 in zip(raw,pred))
  pp=[max(1e-15,min(1-1e-15,v)) for v in raw]
  ll=statistics.mean(-t*math.log(q)-(1-t)*math.log1p(-q) for t,q in zip(y,pp))
  br=statistics.mean((t-q)**2 for t,q in zip(y,raw))
  assert close(ll,z["logloss"]),(ds,z["method"],z["fold"],"ll",ll,z["logloss"])
  assert close(br,z["brier"]),(ds,z["method"],z["fold"],"brier",br,z["brier"])
  assert close(f1(y,pred,1),z["positive_f1"])
  assert close(auc(y,raw),z["auc"]),(ds,z["method"],z["fold"],"auc")
  assert close(ece(y,raw),z["ece"]),(ds,z["method"],z["fold"],"ece")
 for name in methods:
  mm=[z for z in m if z["method"]==name]
  assert len(mm)==5 and {int(z["fold"]) for z in mm}==set(range(5))
  dsrows.append(dict(dataset=ds,method=name,**{
   k:statistics.mean(float(z[k]) for z in mm)
   for k in ("logloss","brier","auc","positive_f1","ece","active_coefficients")}))
 allmetrics+=m;allpatients+=p
assert len(allmetrics)==165 and len(dsrows)==33
dump("all11_165_outer_fold_metrics.csv",allmetrics)
dump("all11_outer_patient_predictions.csv",allpatients)
dump("all11_dataset_equal_metrics.csv",dsrows)
by=collections.defaultdict(dict)
for z in dsrows:by[z["dataset"]][z["method"]]=z
assert len(by)==11 and all(set(g)==set(methods) for g in by.values())
comparisons=[];summaries=[]
for name in methods:
 z=[by[ds][name] for ds in sorted(by)]
 summaries.append(dict(method=name,n_datasets=11,
  **{k:statistics.mean(float(q[k]) for q in z)
  for k in ("logloss","brier","auc","positive_f1","ece","active_coefficients")}))
for target in methods[1:]:
 for metric in ("logloss","brier","auc","positive_f1"):
  diff=[float(by[ds][methods[0]][metric])-float(by[ds][target][metric]) for ds in sorted(by)]
  delta=statistics.mean(diff)
  signs=[abs(sum(t*s for t,s in zip(diff,sign))/len(diff))
         for sign in itertools.product((-1,1),repeat=11)]
  p=sum(q>=abs(delta)-1e-13 for q in signs)/2048
  comparisons.append(dict(comparator=target,metric=metric,
      mean_RLR_minus_binary=delta,exact_two_sided_signflip_p=p,
      n_independent_dataset_units=11,n_outer_folds_per_dataset=5,
      trained_clinical_purist_or_decaf="NO",original_paper_ARFF_SHA256_parity="NOT_VERIFIED"))
# exploratory BH adjustment across all eight tests
order=sorted(range(len(comparisons)),key=lambda i:comparisons[i]["exact_two_sided_signflip_p"])
qraw=[min(1,comparisons[ix]["exact_two_sided_signflip_p"]*len(comparisons)/(i+1))
      for i,ix in enumerate(order)]
minq=1.0
for pos in reversed(range(len(order))):
 minq=min(minq,qraw[pos]);comparisons[order[pos]]["exploratory_BH_q"]=minq
dump("all11_summary.csv",summaries)
dump("all11_paired_exact_inference.csv",comparisons)
for z in summaries:print("OUTER5_SUMMARY",z)
for z in comparisons:print("OUTER5_PAIRED",z)
print("RLR_BINARY_OUTER5_165_FOLD_METRICS_INDEPENDENT_REPLAY=PASS")
print("RLR_BINARY_OUTER5_EVERY_PATIENT_EXACTLY_ONCE_PER_METHOD=PASS")
print("RLR_BINARY_OUTER5_EXACT_2048_SIGNFLIPS_AND_BH8=PASS")
