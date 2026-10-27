#!/usr/bin/env python3
"""Independent 11-cohort COT RLR vs real KTSP50 first-fold E audit."""
from __future__ import annotations
import csv,math,argparse
from pathlib import Path
P=argparse.ArgumentParser()
P.add_argument("--input",required=True)
P.add_argument("--out",required=True)
A=P.parse_args()
DATASETS=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
ARMS=("KTSP50_NATIVE","ORIGINAL_RLR50_THRESHOLD05",
      "RLR50_DOOF_GLOBAL_THRESHOLD",
      "RLR50_ORDINAL_TREE_DIAGNOSTIC",
      "RLR50_ORDINAL_TREE_V_GATED")
def fail(s):raise SystemExit("COT_AUDIT_STOP "+s)
def rd(path):
  with path.open(encoding="utf-8",newline="") as f:
    return list(csv.DictReader(f))
def f1(y,h,c):
  tp=sum(a==c and b==c for a,b in zip(y,h))
  fp=sum(a!=c and b==c for a,b in zip(y,h))
  fn=sum(a==c and b!=c for a,b in zip(y,h))
  return 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0
def roc(y,p):
  k1=y.count(1);k0=y.count(0)
  idx=sorted(range(len(p)),key=lambda i:p[i])
  ranks=[0.0]*len(p)
  i=0
  while i<len(p):
    j=i+1
    while j<len(p) and p[idx[j]]==p[idx[i]]:j+=1
    rank=(i+1+j)/2
    for q in range(i,j):ranks[idx[q]]=rank
    i=j
  return (sum(ranks[i] for i,a in enumerate(y) if a==1)-
          k1*(k1+1)/2)/(k0*k1)
def metrics(y,p,h):
  n1=y.count(1);n0=y.count(0)
  if not n1 or not n0:fail("one class missing E")
  return dict(
    positive_f1=f1(y,h,1),
    macro_f1=.5*(f1(y,h,1)+f1(y,h,0)),
    ba=.5*(sum(a==1 and b==1 for a,b in zip(y,h))/n1+
           sum(a==0 and b==0 for a,b in zip(y,h))/n0),
    auc=roc(y,p),brier=sum((a-b)**2 for a,b in zip(y,p))/len(y),
    logloss=sum(-(a*math.log(max(1e-15,min(1-1e-15,b)))+
                  (1-a)*math.log1p(-max(1e-15,min(1-1e-15,b))))
                for a,b in zip(y,p))/len(y))
rows=[]
summary=[]
for ds in DATASETS:
  files=list(Path(A.input).rglob(f"{ds}/method_metrics.csv"))
  if len(files)!=1:fail(f"expected one method file for {ds}; got {len(files)}")
  origin=files[0].parent
  m=rd(origin/"method_metrics.csv")
  pred=rd(origin/"patient_predictions.csv")
  manifest=rd(origin/"split_manifest.csv")
  sel=rd(origin/"decision_selection.csv")
  if len(m)!=5 or len(manifest)!=1 or len(sel)!=1:fail(f"bad metadata {ds}")
  fold=manifest[0]
  if fold["split_id"]!="rep00_fold00":fail("wrong original split "+ds)
  T=set(fold["T_ids"].split("|"));E=set(fold["E_ids"].split("|"))
  nE=int(fold["n_E"])
  if T&E or len(E)!=nE or len(T)!=int(fold["n_train"]):
    fail("train/E leakage or invalid N "+ds)
  if len(pred)!=5*nE:fail("patient prediction count "+ds)
  byarm={r["method"]:r for r in m}
  if set(byarm)!=set(ARMS):fail("method set mismatch "+ds)
  arr={}
  for arm in ARMS:
    pp=[x for x in pred if x["method"]==arm]
    data={x["sample_id"]:x for x in pp}
    if len(pp)!=len(data) or set(data)!=E:fail("patient alignment "+ds+"/"+arm)
    y=[];p=[];h=[]
    for sid in sorted(E):
      row=data[sid]
      yi=int(row["truth"]);pi=float(row["prob"]);hi=int(row["label"])
      if yi not in (0,1) or hi not in (0,1) or not math.isfinite(pi) or not 0<=pi<=1:
        fail("patient data invalid "+ds+"/"+arm)
      y.append(yi);p.append(pi);h.append(hi)
    arr[arm]={"y":y,"p":p,"h":h,"data":data}
    v=metrics(y,p,h)
    for name,num in v.items():
      if abs(float(byarm[arm][name])-num)>1e-9:fail("metric mismatch "+ds+"/"+arm+"/"+name)
    rows.append(byarm[arm])
  base=arr["ORIGINAL_RLR50_THRESHOLD05"]
  for arm in ARMS:
    x=arr[arm]
    if x["y"]!=base["y"]:fail("truth mismatch "+ds+"/"+arm)
    if arm=="KTSP50_NATIVE":continue
    if x["p"]!=base["p"]:fail("changed RLR probabilities "+ds+"/"+arm)
    if arm=="ORIGINAL_RLR50_THRESHOLD05":
      if x["h"]!=[int(v>=.5) for v in x["p"]]:fail("native label mismatch")
    if arm=="RLR50_DOOF_GLOBAL_THRESHOLD":
      t=float(sel[0]["global_threshold"])
      if x["h"]!=[int(v>=t) for v in x["p"]]:fail("global threshold mismatch")
    for sid,z in x["data"].items():
      if float(z["RLR50_original_prob"])!=float(base["data"][sid]["prob"]):
        fail("parent p mismatch "+ds+"/"+arm)
  acc=int(sel[0]["accepted"])
  if acc not in (0,1):fail("invalid gate")
  expected=arr["RLR50_ORDINAL_TREE_DIAGNOSTIC"] if acc else arr["RLR50_DOOF_GLOBAL_THRESHOLD"]
  if arr["RLR50_ORDINAL_TREE_V_GATED"]["h"]!=expected["h"]:
    fail("tree acceptance fallback not exact "+ds)
  if acc:
    if float(sel[0]["V_F1_gain"])<.02 or float(sel[0]["V_BA_gain"])<-.02:
      fail("accepted tree with invalid V gate "+ds)
  if any(x["auc"]!=byarm["ORIGINAL_RLR50_THRESHOLD05"]["auc"] or
         x["logloss"]!=byarm["ORIGINAL_RLR50_THRESHOLD05"]["logloss"] or
         x["brier"]!=byarm["ORIGINAL_RLR50_THRESHOLD05"]["brier"]
         for a,x in byarm.items() if a!="KTSP50_NATIVE"):
    fail("probability score metrics changed "+ds)
  summary.append(dict(dataset=ds,n_E=nE,accepted=acc,
    inner_V_F1_gain=sel[0]["V_F1_gain"],
    F1_kTSP50=float(byarm["KTSP50_NATIVE"]["positive_f1"]),
    F1_RLR05=float(byarm["ORIGINAL_RLR50_THRESHOLD05"]["positive_f1"]),
    F1_RLRglobal=float(byarm["RLR50_DOOF_GLOBAL_THRESHOLD"]["positive_f1"]),
    F1_RLRtree=float(byarm["RLR50_ORDINAL_TREE_V_GATED"]["positive_f1"]),
    Macro_RLRtree=float(byarm["RLR50_ORDINAL_TREE_V_GATED"]["macro_f1"]),
    AUC_RLR=float(byarm["ORIGINAL_RLR50_THRESHOLD05"]["auc"])))
out=Path(A.out);out.mkdir(parents=True,exist_ok=True)
with (out/"all11_firstfold_independent_audit.csv").open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
with (out/"all11_method_fold_metrics.csv").open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
print("COT_11_COHORTS_MATCHED_CANONICAL_FIRST_E_FOLD=PASS")
print("COT_55_METHOD_FOLD_METRICS_INDEPENDENTLY_REPLAYED=PASS")
print("COT_RLR50_PARENT_PROBABILITY_AUC_LOGLOSS_BRIER_IDENTICAL=PASS")
print("COT_INNER_V_GATE_AND_GLOBAL_THRESHOLD_FALLBACK=PASS")
print(f"COT_V_ACCEPTED_TREE_COHORTS={sum(v['accepted'] for v in summary)}/11")
for arm in ARMS:
  z=[float(x["positive_f1"]) for x in rows if x["method"]==arm]
  print(f"COT_SUMMARY {arm} mean_positive_F1={sum(z)/len(z):.6f}")
for r in summary:
  print(f"COT_CASE {r['dataset']} accept={r['accepted']} kTSP={r['F1_kTSP50']:.3f} RLR05={r['F1_RLR05']:.3f} Global={r['F1_RLRglobal']:.3f} Gated={r['F1_RLRtree']:.3f}")
print("COT_EXPOSED_PUBLIC_SOURCE_FIRST_FOLD_ONLY_NOT_ORIGINAL470")
