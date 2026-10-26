#!/usr/bin/env python3
"""Independent heldout patient/BA/macro-F1/logloss Brier replay, T-GRLR V1."""
import csv,math,argparse,json
from pathlib import Path
from collections import defaultdict
N={"chen-2002":179,"gordon-2002":181,"yeoh-2002-v1":248}
METHODS=("ORIGINAL_RLR_K10","ORIGINAL_RLR_K50","TSP_1","kTSP_9",
         "TSPDT_GREEDY_PROXY","RGT50_GATED_RESIDUAL_TREE")
def need(cond,msg):
 if not cond:raise SystemExit("RGT50_AUDIT_STOP: "+msg)
def read(p):
 need(p.is_file(),"missing "+str(p))
 with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def stats(y,p):
 # BA, macro-F1 and Brier are computed from ORIGINAL raw probabilities.
 # Only logarithmic loss is clipped, exactly as in the R implementation.
 pc=[min(1-1e-8,max(1e-8,v)) for v in p]
 z=[int(v>=.5) for v in p]
 ba=sum(sum(y[i]==c and z[i]==c for i in range(len(y)))/
      max(1,sum(v==c for v in y)) for c in (0,1))/2
 f1=sum(2*sum(y[i]==c and z[i]==c for i in range(len(y)))/
        max(1,sum(v==c for v in y)+sum(v==c for v in z))
        for c in (0,1))/2
 ll=sum(-(v*math.log(t)+(1-v)*math.log1p(-t))
         for v,t in zip(y,pc))/len(y)
 br=sum((v-t)**2 for v,t in zip(y,p))/len(y)
 return dict(BA=ba,macro_F1=f1,logloss=ll,brier=br)
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--dataset",required=True,choices=list(N))
 ds=ap.parse_args().dataset
 root=Path("results/rlr_rgt50_v1")/ds
 m=read(root/"method_metrics.csv")
 p=read(root/"patient_predictions.csv")
 folds=read(root/"split_manifest.csv")
 struct=json.loads((root/"tree_structure.json").read_text())
 expected={(f"rep00_fold{j:02d}",arm) for j in range(3) for arm in METHODS}
 idx={(r["split_id"],r["method"]):r for r in m}
 need(len(idx)==len(m)==18 and set(idx)==expected,"method-fold inventory")
 groups=defaultdict(list)
 for r in p:
   key=(r["split_id"],r["method"])
   need(key in expected,"unknown arm")
   groups[key].append(r)
 need(set(groups)==expected,"missing E patient groups")
 all_e=set();npred=0
 for fold in folds:
   fid=fold["split_id"];E=set(fold["E_ids"].split("|"))
   need(len(E)==int(fold["n_E"]) and
         len(E & all_e)==0,"duplicate E patients")
   all_e.update(E)
   orig={q["sample_id"]:q for q in groups[fid,"ORIGINAL_RLR_K50"]}
   need(set(orig)==E,"original K50 E partition discrepancy")
   for arm in METHODS:
     rr=groups[fid,arm]
     need(len(rr)==len(E) and {q["sample_id"] for q in rr}==E,
          "missing or repeated matched E patient")
     truth=[];probs=[]
     for q in rr:
       y=int(q["truth"]);prob=float(q["prob"])
       need(y in (0,1) and math.isfinite(prob) and 0<=prob<=1,
            "nonbinary E or bad probability")
       need(int(orig[q["sample_id"]]["truth"])==y,"positive-class mismatch")
       truth.append(y);probs.append(prob)
       npred+=1
     computed=stats(truth,probs)
     met=idx[fid,arm]
     for col,value in computed.items():
       need(abs(value-float(met[col]))<1e-8,
          f"patient metric replay mismatch {fid}/{arm}/{col}")
     need(met["scientific_stage"]==
          "RGT50_ORIGINAL_PAPER_RECONSTRUCTED_SOURCE_ARFF_NOT_VERIFIED",
          "incorrect claim of original ARFF parity")
   need(fid in struct and "V_gate_accepted" in struct[fid],
        "missing D/V structure and gate")
   accepted=str(struct[fid]["V_gate_accepted"]).lower()=="true"
   if not accepted:
     pa={q["sample_id"]:q for q in groups[fid,"RGT50_GATED_RESIDUAL_TREE"]}
     need(all(abs(float(pa[z]["prob"])-float(orig[z]["prob"]))<1e-8
              for z in E),"fallback does not equal original K50")
 need(len(folds)==3 and len(all_e)==N[ds] and npred==6*N[ds],
     "E partition not complete / predicted patient models")
 print("RGT50_INDEPENDENT_DATASET",ds)
 for arm in METHODS:
   arr=[idx[f"rep00_fold{i:02d}",arm] for i in range(3)]
   print(f"RGT50_AUDIT_ARM={arm} BA={sum(float(r['BA']) for r in arr)/3:.5f} macroF1={sum(float(r['macro_F1']) for r in arr)/3:.5f} logloss={sum(float(r['logloss']) for r in arr)/3:.5f}")
 print("RGT50_AUDIT_E_PREDICTION_PAIRS",npred)
 print("RGT50_INDEPENDENT_3_FOLD_6_METHOD_SCORE_PARITY=PASS")
 print("RGT50_PATIENT_TEST_SETS_DISJOINT_AND_COMPLETE=PASS")
 print("RGT50_ORIGINAL_PROCESSED_ARFF_SHA=NOT_VERIFIED")
if __name__=="__main__":main()
