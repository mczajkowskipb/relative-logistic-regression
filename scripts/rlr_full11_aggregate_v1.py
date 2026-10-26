#!/usr/bin/env python3
"""Aggregate 11-dataset 470-outer-fold same-pair RLR component comparison.
Dataset, not repeat/fold, is unit for any exploratory Wilcoxon comparison.
Strictly labels source reconstructions NOT exact original ARFF provenance.
"""
from __future__ import annotations
import argparse,csv,json,math,itertools
from collections import defaultdict
from pathlib import Path
SETS=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
METHODS=("ORIGINAL_RLR","REGULARIZED_EQUAL_VOTE")
K=(10,20,50)
def need(flag,msg):
 if not flag:raise SystemExit("FULL11_AGGREGATE_STOP: "+msg)
def read_csv(file):
 with file.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))
def write_csv(file,rows):
 with file.open("w",newline="",encoding="utf8") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]))
  w.writeheader();w.writerows(rows)
def exact_signed_rank_p(diffs):
  """Exact two-sided signed-rank sign-permutation p for n<=11 datasets.
  Drops exact zeros, assigns midranks to tied absolute effects, enumerates
  all 2^n possible sign assignments (no asymptotic approximation).
  """
  x=[float(v) for v in diffs if abs(v)>0]
  n=len(x)
  if n==0:return 1.0
  assert n<=11
  magnitudes=[abs(v) for v in x]
  order=sorted(range(n),key=lambda i:magnitudes[i])
  ranks=[0.]*n
  j=0
  while j<n:
    end=j+1
    while end<n and magnitudes[order[end]]==magnitudes[order[j]]:
      end+=1
    mid=(j+1+end)/2
    for k in range(j,end):ranks[order[k]]=mid
    j=end
  observed=sum(ranks[i] for i in range(n) if x[i]>0)
  center=sum(ranks)/2
  extreme=0
  for signs in itertools.product((0,1),repeat=n):
    total=sum(ranks[i] for i in range(n) if signs[i])
    if abs(total-center)>=abs(observed-center)-1e-12:
      extreme+=1
  return extreme/(2**n)

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--input",default="downloaded")
 p.add_argument("--out",default="results/rlr_full11_aggregate_v1")
 p.add_argument("--mode",choices=("full","smoke"),default="full")
 a=p.parse_args()
 f=list(Path(a.input).rglob("method_metrics.csv"))
 need(len(f)==11,f"exactly 11 dataset metric sets required; found {len(f)}")
 per={}
 for file in f:
  rows=read_csv(file)
  need(bool(rows),"empty study results")
  ds=rows[0]["dataset"]
  need(ds in SETS and ds not in per,"dataset duplication or off-paper name")
  need(all(x["dataset"]==ds and x["scientific_stage"]==
    "SOURCE_RECONSTRUCTION_NOT_ARFF_PARITY" for x in rows),
    "provenance or dataset record inconsistent")
  need(all(x["method"] in METHODS and int(x["K"]) in K for x in rows),
    "unknown compared arm")
  expect_splits=(10 if ds in ("west-2001","yeoh-2002-v1") else 50)
  if a.mode=="smoke":expect_splits=1
  pairs={(x["split_id"],int(x["K"]),x["method"]) for x in rows}
  need(len(rows)==len(pairs)==expect_splits*3*2,
       "missing/duplicate experiment splits")
  rr=sorted({x["split_id"] for x in rows})
  expected=[f"rep{r:02d}_fold{j:02d}" for r in range(expect_splits//5)
              for j in range(5)] if a.mode=="full" else ["rep00_fold00"]
  need(rr==expected,f"not original per-cohort CV schedule for {ds}")
  per[ds]=rows
 need(set(per)==set(SETS),"publication inventory not fully represented")
 output=Path(a.out); output.mkdir(parents=True,exist_ok=True)
 item=[]
 for ds in SETS:
  rows=per[ds]
  for k in K:
   baseline=[x for x in rows if x["method"]=="ORIGINAL_RLR" and int(x["K"])==k]
   equal=[x for x in rows if x["method"]=="REGULARIZED_EQUAL_VOTE" and int(x["K"])==k]
   need(len(equal)==len(baseline)>0,"missing matched fold method")
   dd=[float(x["delta_logloss"]) for x in equal]
   bb={x["split_id"]:float(x["brier"]) for x in baseline}
   bdelta=[float(x["brier"])-bb[x["split_id"]] for x in equal]
   need(all(math.isfinite(v) for v in dd+bdelta),"nonfinite comparison")
   item.append(dict(dataset=ds,K=k,n_outer_folds=len(equal),
       mean_delta_logloss=sum(dd)/len(dd),
       mean_delta_brier=sum(bdelta)/len(bdelta),
       n_RLR_better=sum(z>0 for z in dd),
       n_equal_vote_better=sum(z<0 for z in dd),
       source_parity="ORIGINAL_ARFF_SHA256_NOT_VERIFIED"))
 write_csv(output/"dataset_level.csv",item)
 summary=[]
 for k in K:
  vals=[next(x for x in item if x["dataset"]==ds and x["K"]==k) for ds in SETS]
  d=[x["mean_delta_logloss"] for x in vals]
  summary.append(dict(K=k,n_datasets=11,
    n_total_outer_folds=sum(x["n_outer_folds"] for x in vals),
    dataset_mean_delta_logloss=sum(d)/len(d),
    dataset_median_delta_logloss=sorted(d)[len(d)//2],
    n_RLR_better_datasets=sum(x>0 for x in d),
    n_equal_vote_better_datasets=sum(x<0 for x in d),
    wilcoxon_two_sided_p_raw=exact_signed_rank_p(d),
    scientific_stage="EXPLORATORY_RECONSTRUCTED_SOURCE_NO_ARFF_PARITY"))
 ps=[x["wilcoxon_two_sided_p_raw"] for x in summary]
 order=sorted(range(len(ps)),key=lambda j:ps[j])
 adjusted=[0.0]*len(ps); running=0.
 for rank,idx in enumerate(order):
  running=max(running,min(1.,(len(ps)-rank)*ps[idx]))
  adjusted[idx]=running
 for i in range(len(ps)):summary[i]["holm_p_across_3_exploratory_K"]=adjusted[i]
 write_csv(output/"summary_exploratory.csv",summary)
 (output/"manifest.json").write_text(json.dumps(dict(
    datasets=list(SETS),K=list(K),source_type="PUBLIC_SCHLIEP_FILTERED_TEXT",
    original_ARFF_sha256_match=False,mode=a.mode,
    total_matched_folds=sum(x["n_outer_folds"] for x in item)//3,
    paired_fits=sum(len(v) for v in per.values())),indent=2)+"\n")
 for r in summary:print("FULL11_AGGREGATE_RESULT",json.dumps(r))
 print("FULL11_AGGREGATE_11_BINARY_PUBLICATION_DATASETS=PASS")
 if a.mode=="full":
  need(all(x["n_total_outer_folds"]==470 for x in summary),
       "not full frozen CV count")
  print("FULL11_AGGREGATE_470_CV_FOLDS_PER_MODEL=PASS")
 print("FULL11_AGGREGATE_ORIGINAL_ARFF_SHA256_PARITY=NOT_VERIFIED")
if __name__=="__main__":main()
