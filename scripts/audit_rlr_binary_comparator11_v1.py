#!/usr/bin/env python3
"""Independent replay, 11 cohort equal-weight summaries, exact paired 2^11 inference."""
import argparse,csv,glob,math,statistics,itertools,pathlib,collections
p=argparse.ArgumentParser()
p.add_argument("--input",required=True);p.add_argument("--out",required=True);a=p.parse_args()
out=pathlib.Path(a.out);out.mkdir(parents=True,exist_ok=True)
files=sorted(glob.glob(a.input+"/**/metrics.csv",recursive=True))
assert len(files)==11,(len(files),files)
MET=[];PAT=[]
arms={"RLR20_ORIGINAL","BINARY_ENET20_PURIST_INSPIRED","BINARY_MCP20_DECAF_INSPIRED"}
for file in files:
 d=pathlib.Path(file).parent
 with open(file,newline="") as h: m=list(csv.DictReader(h))
 with open(d/"patients.csv",newline="") as h: pt=list(csv.DictReader(h))
 with open(d/"manifest.csv",newline="") as h: manifest=list(csv.DictReader(h))
 assert len(m)==3 and len(manifest)==1 and set(z["method"] for z in m)==arms
 assert int(manifest[0]["n_outer_E_UNREAD"])>0 and len(pt)==3*int(manifest[0]["n_V"])
 ids=set(z["sample_id"] for z in pt)
 assert len(ids)==int(manifest[0]["n_V"])
 for z in m:
  ps=[q for q in pt if q["method"]==z["method"]]
  y=[int(q["truth"]) for q in ps]; prob=[float(q["prob"]) for q in ps]
  assert len(y)==int(z["n_V"]) and set(y)=={0,1}
  clipped=[min(1-1e-15,max(1e-15,x)) for x in prob]
  ll=sum(-t*math.log(x)-(1-t)*math.log1p(-x) for t,x in zip(y,clipped))/len(y)
  brier=sum((t-x)**2 for t,x in zip(y,clipped))/len(y)
  assert math.isclose(ll,float(z["logloss"]),rel_tol=1e-8,abs_tol=1e-9)
  assert math.isclose(brier,float(z["brier"]),rel_tol=1e-8,abs_tol=1e-9)
  assert all(int(q["pred"])==int(q["prob"] and float(q["prob"])>=.5) for q in ps)
 MET.extend(m);PAT.extend(pt)
assert len(MET)==33
by=collections.defaultdict(dict)
for z in MET:by[z["dataset"]][z["method"]]=z
assert len(by)==11 and all(set(z)==arms for z in by.values())
def dump(path,rows):
 with open(path,"w",newline="") as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
dump(out/"all11_metrics.csv",MET)
dump(out/"all11_patients.csv",PAT)
rows=[]
for target in ("BINARY_ENET20_PURIST_INSPIRED","BINARY_MCP20_DECAF_INSPIRED"):
 for metric in ("logloss","brier","auc","positive_f1"):
  diffs=[float(z["RLR20_ORIGINAL"][metric])-float(z[target][metric]) for _,z in sorted(by.items())]
  mean=statistics.mean(diffs)
  vals=[abs(sum(d*s for d,s in zip(diffs,sign))/11) for sign in itertools.product((-1,1),repeat=11)]
  pv=sum(x>=abs(mean)-1e-13 for x in vals)/len(vals)
  rows.append(dict(comparator=target,metric=metric,mean_RLR_minus_binary=mean,
   exact_two_sided_signflip_p=pv,n_independent_cohorts=11,exploratory="YES",original_ARFF_sha_parity="NOT_VERIFIED"))
dump(out/"paired_inference.csv",rows)
for arm in sorted(arms):
 vals=[z[arm] for z in by.values()]
 print("METHOD",arm," ".join(k+"="+str(round(statistics.mean(float(z[k]) for z in vals),6)) for k in ("positive_f1","auc","logloss","brier","ece","active_coefficients")))
for z in rows:print("PAIR",z["comparator"],z["metric"],z["mean_RLR_minus_binary"],z["exact_two_sided_signflip_p"])
print("RLR_BINARY_ALL11_THREE_METHODS_PATIENT_REPLAY=PASS")
print("RLR_BINARY_EXACT_SIGNFLIP_2048=PASS")
