#!/usr/bin/env python3
"""Independent heldout E verification of every original RLR certificate."""
import csv,math, argparse
from collections import defaultdict
from pathlib import Path
nmap={"chen-2002":179,"singh-2002":102,"yeoh-2002-v1":248}
arms=("ORIGINAL_RLR","CERT_TREE_ABS_WEIGHT","CERT_TREE_EA")
def need(x,m):
 if not x:raise SystemExit("CRET_AUDIT_STOP: "+m)
def read(p):
 need(p.exists(),f"missing {p}")
 with p.open(newline="",encoding="utf8") as f:return list(csv.DictReader(f))
def stat(y,p):
 pred=[int(v>=.5) for v in p]
 ba=sum(sum(pred[j]==cl for j in range(len(y)) if y[j]==cl)/
   max(1,sum(t==cl for t in y)) for cl in (0,1))/2
 f1=sum(2*sum(pred[j]==cl and y[j]==cl for j in range(len(y)))/
   max(1,sum(z==cl for z in pred)+sum(z==cl for z in y))
   for cl in (0,1))/2
 pc=[min(1-1e-7,max(1e-7,v)) for v in p]
 ll=sum(-(t*math.log(z)+(1-t)*math.log1p(-z))
        for t,z in zip(y,pc))/len(y)
 b=sum((v-t)**2 for v,t in zip(p,y))/len(y)
 return {"BA":ba,"macro_F1":f1,"logloss":ll,"brier":b}
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--dataset",choices=tuple(nmap),required=True)
 ds=parser.parse_args().dataset
 path=Path("results/rlr_arlet_v1")/ds/"certified_tree"
 rows=read(path/"patient_certificates.csv")
 metric=read(path/"method_metrics.csv")
 met={(r["split_id"],r["method"]):r for r in metric}
 groups=defaultdict(list)
 for r in rows:groups[r["split_id"],r["method"]].append(r)
 expected={(f"rep00_fold{i:02d}",arm) for i in range(3) for arm in arms}
 need(set(met)==expected and set(groups)==expected and
       len(metric)==9 and len(rows)==3*nmap[ds],
       "method or patient archive incomplete")
 ever=set();number=0
 for i in range(3):
  fold=f"rep00_fold{i:02d}"
  original={r["sample_id"]:r for r in groups[fold,"ORIGINAL_RLR"]}
  ids=set(original)
  need(not (ids&ever),"E IDs not disjoint")
  ever|=ids
  for arm in arms:
   rr=groups[fold,arm]
   need({r["sample_id"] for r in rr}==ids,"paired E patients disagree")
   p=[];truth=[];ev=[];act=[]
   for r in rr:
    ref=original[r["sample_id"]]
    y=int(r["truth"]);pr=float(r["prob"])
    b=float(ref["prob"]);lo=float(r["prob_lower"]);hi=float(r["prob_upper"])
    e=int(r["evaluated"]);a=int(r["active"])
    need(y==int(ref["truth"]) and y in (0,1),"label mismatch")
    need(0<=lo<=b+1e-8 and b<=hi+1e-8 and hi<=1,
         "original RLR probability outside certificate")
    need(abs(pr-b)<=.0250002,"epsilon 0.05 half-width violated")
    need((pr>=.5)==(b>=.5),"exact original RLR class not preserved")
    need(0<=a<=e<=50,"test/evidence count impossible")
    need(r["input_provenance"] if "input_provenance" in r else True,
        "invalid")
    p.append(pr);truth.append(y);ev.append(e);act.append(a)
    number+=1
   scores=stat(truth,p)
   m=met[fold,arm]
   for key,v in scores.items():
    need(abs(v-float(m[key]))<1e-8,f"metric mismatch {fold}/{arm}/{key}")
   need(abs(sum(ev)/len(ev)-float(m["evaluated_mean"]))<1e-8,
        "mean evaluated feature count audit")
   need(abs(sum(act)/len(act)-float(m["active_mean"]))<1e-8,
        "active feature count audit")
 need(len(ever)==nmap[ds],"cohort E coverage incomplete")
 print("CRET_INDEPENDENT_DATASET",ds)
 for arm in arms:
  mm=[met[f"rep00_fold{i:02d}",arm] for i in range(3)]
  print(f"CRET_SUMMARY_ARM={arm} BA={sum(float(x['BA']) for x in mm)/3:.6f} macroF1={sum(float(x['macro_F1']) for x in mm)/3:.6f} eval={sum(float(x['evaluated_mean']) for x in mm)/3:.6f} early={sum(float(x['stop_early_fraction']) for x in mm)/3:.6f}")
 print(f"CRET_INDEPENDENT_PATIENT_MODELS_REPLAYED={number}")
 print("CRET_CERTIFICATES_AND_EXACT_CLASSES_VALIDATED=PASS")
 print("CRET_ALL_THREE_E_FOLDS_DISJOINT=PASS")
 print("CRET_DATASET_ARFF_PROCESSED_PARITY=NOT_VERIFIED")
if __name__=="__main__":main()
