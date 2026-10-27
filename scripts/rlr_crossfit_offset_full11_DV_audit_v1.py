#!/usr/bin/env python3
"""Independent patient-level D/V checker of 11-compCancer CFO model trees."""
from __future__ import annotations
import argparse,csv,math,collections
from pathlib import Path
P=argparse.ArgumentParser()
P.add_argument("--input",required=True)
P.add_argument("--out",required=True)
o=P.parse_args()
COHORTS=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ARMS=("ORIGINAL_RLR50","RESIDUAL_LOCAL_INSAMPLE_OFFSET",
 "RESIDUAL_LOCAL_CROSSFIT_OFFSET",
 "SAMEPAIRS_GLOBAL_INSAMPLE_OFFSET",
 "SAMEPAIRS_GLOBAL_CROSSFIT_OFFSET")
def stop(s):raise SystemExit("CFO_FULL11_AUDIT_STOP "+s)
def read(p):
 with p.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def ll(y,p):
 return sum(-(a*math.log(max(1e-15,min(1-1e-15,b)))+
               (1-a)*math.log1p(-max(1e-15,min(1-1e-15,b))))
            for a,b in zip(y,p))/len(y)
def scores(y,p):
 h=[int(x>=.5) for x in p]
 def f(cl):
  tp=sum(a==cl and b==cl for a,b in zip(y,h))
  fp=sum(a!=cl and b==cl for a,b in zip(y,h))
  fn=sum(a==cl and b!=cl for a,b in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.
 if not y.count(0) or not y.count(1):stop("one V class absent")
 return dict(LL=ll(y,p),macroF1=(f(0)+f(1))/2,positiveF1=f(1),
    BA=(sum(a==0 and b==0 for a,b in zip(y,h))/y.count(0)+
        sum(a==1 and b==1 for a,b in zip(y,h))/y.count(1))/2)
results=[]
for ds in COHORTS:
 locations=list(Path(o.input).rglob(f"{ds}/selection_DV.csv"))
 if len(locations)!=1:stop(f"cohort selection missing or duplicate {ds}:{len(locations)}")
 meta=read(locations[0])
 if len(meta)!=1:stop("multiple metadata rows "+ds)
 z=meta[0]
 if z["dataset"]!=ds or z["original_outer_split"]!="rep00_fold00" or z["evidence"]!="ORIGINAL_NAME_PUBLIC_SOURCE_D_V_ONLY_NO_OUTER_E":stop("wrong study stage "+ds)
 nT=int(z["n_T"]);nD=int(z["n_D"]);nV=int(z["n_V"]);nE=int(z["n_outer_E_UNREAD"])
 if nD+nV!=nT or nT+nE!=int(z["n_total"]):stop("split inventory "+ds)
 accepted=int(z["accepted"])
 if accepted not in (0,1) or bool(accepted)!=(z["state"]=="INNER_V_ACCEPTED"):stop("accept state "+ds)
 rows=list(locations[0].parent.glob("patient_V.csv"))
 if z["state"] in ("SAFE_FALLBACK_SMALL_D","SAFE_FALLBACK_SMALL_OOF","SAFE_FALLBACK_NO_ROOT","ERROR_NEEDS_DIAGNOSIS"):
  if accepted or rows or int(z["n_roots_evaluated"])!=0:stop("fake skip "+ds)
 else:
  if z["state"] not in ("INNER_V_ACCEPTED","INNER_V_REJECTED"):stop("unexpected state "+ds)
  if len(rows)!=1:stop("missing patient scores "+ds)
  d=read(rows[0])
  g=collections.defaultdict(list)
  for q in d:g[(int(q["root_id"]),q["method"])].append(q)
  roots=sorted(set(k[0] for k in g))
  if roots!=list(range(1,len(roots)+1)) or len(roots)!=int(z["n_roots_evaluated"]):stop("root grid "+ds)
  if set(a for _,a in g)!=set(ARMS):stop("model set "+ds)
  allmetrics={}
  for root in roots:
   yref=None;idsref=None
   for arm in ARMS:
    rr=g[root,arm]
    ids=[q["sample_id"] for q in rr]
    truth=[int(q["truth"]) for q in rr]
    pp=[float(q["prob"]) for q in rr]
    if len(rr)!=nV or len(set(ids))!=nV:stop("V patient duplication "+ds)
    if any(not math.isfinite(x) or x<0 or x>1 for x in pp):stop("invalid p "+ds)
    if yref is None:yref=truth;idsref=ids
    elif yref!=truth or idsref!=ids:stop("method patient misalignment "+ds)
    allmetrics[root,arm]=scores(truth,pp)
  j=int(z["best_root_id"])
  if j not in roots:stop("chosen root id "+ds)
  best=min(roots,key=lambda k:allmetrics[k,"RESIDUAL_LOCAL_CROSSFIT_OFFSET"]["LL"])
  if j!=best:stop("not min V tree "+ds)
  orig=allmetrics[j,"ORIGINAL_RLR50"]
  tree=allmetrics[j,"RESIDUAL_LOCAL_CROSSFIT_OFFSET"]
  glob=allmetrics[j,"SAMEPAIRS_GLOBAL_CROSSFIT_OFFSET"]
  cmp={"V_original_LL":orig["LL"],
       "V_tree_LL":tree["LL"],
       "V_samepairs_global_LL":glob["LL"],
       "V_gain_vs_original":orig["LL"]-tree["LL"],
       "V_gain_vs_global":glob["LL"]-tree["LL"],
       "V_original_MacroF1":orig["macroF1"],
       "V_tree_MacroF1":tree["macroF1"],
       "V_tree_minus_original_BA":tree["BA"]-orig["BA"]}
  for name,val in cmp.items():
   if abs(float(z[name])-val)>1e-9:stop(f"metric mismatch {ds}/{name}")
  decision=(cmp["V_gain_vs_original"]>=.005 and
            cmp["V_gain_vs_global"]>=.002 and
            cmp["V_tree_minus_original_BA"]>=-.02 and
            cmp["V_tree_MacroF1"]>=cmp["V_original_MacroF1"]-.01)
  if decision!=bool(accepted):stop("V gate not reproduced "+ds)
 results.append(z)
path=Path(o.out);path.mkdir(parents=True,exist_ok=True)
with (path/"all11_training_only_V_gate.csv").open("w",newline="",encoding="utf-8") as f:
 w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
state=collections.Counter(z["state"] for z in results)
print("CFO_FULL11_ALL_11_ORIGINAL_NAMED_DV_COHORTS=PASS")
print("CFO_FULL11_EXACT_DV_PATIENT_PREDICTIONS_AND_GATES_REPLAYED=PASS")
print("CFO_FULL11_ALL_ORIGINAL_OUTER_E_OUTCOMES_UNREAD=PASS")
print("CFO_FULL11_STATES "+str(dict(state)))
print("CFO_FULL11_V_ACCEPTED="+str(state.get("INNER_V_ACCEPTED",0))+"/11")
for r in results:
 print("CFO_FULL11_CASE "+r["dataset"]+
       " state="+r["state"]+
       " V_gain_RLR="+r["V_gain_vs_original"]+
       " V_gain_samepairs="+r["V_gain_vs_global"]+
       " V_macroF1_orig="+r["V_original_MacroF1"]+
       " V_macroF1_tree="+r["V_tree_MacroF1"])
print("CFO_FULL11_PAPER_ARFF_SHA_PARITY_UNVERIFIED_AND_NOT_FULL470")
