#!/usr/bin/env python3
"""Independent E-patient replay of C-RET-L2 mathematical class/posterior bounds."""
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
N={"alizadeh-2000-v1":42,"pomeroy-2002-v1":34}
ARMS=("ORIGINAL_RLR","CERT_STATIC_ABS","CERT_MYOPIC_CONTEXT",
      "CERT_TREE_DEPTH2_LOOKAHEAD")
def check(ok,msg):
 if not ok:raise SystemExit("CL2_INDEPENDENT_AUDIT_STOP: "+msg)
def read(p):
 check(p.is_file(),"missing "+str(p))
 with p.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def metric(y,p):
 pred=[int(v>=.5) for v in p]
 ba=sum(sum(pred[i]==c for i in range(len(y)) if y[i]==c)/
          max(1,sum(t==c for t in y)) for c in (0,1))/2
 f1=sum(2*sum(pred[i]==c and y[i]==c for i in range(len(y)))/
        max(1,sum(pred[i]==c for i in range(len(y)))+
              sum(y[i]==c for i in range(len(y)))) for c in (0,1))/2
 pc=[max(1e-7,min(1-1e-7,v)) for v in p]
 ll=sum(-(yy*math.log(v)+(1-yy)*math.log1p(-v))
        for yy,v in zip(y,pc))/len(y)
 br=sum((v-yy)**2 for yy,v in zip(y,p))/len(y)
 return dict(BA=ba,macro_F1=f1,logloss=ll,brier=br)
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--dataset",required=True,choices=tuple(N))
 ds=ap.parse_args().dataset
 root=Path("results/rlr_cl2_v1")/ds
 mets=read(root/"metrics.csv")
 patients=read(root/"patient_certificates.csv")
 structures=json.loads((root/"tree_structures.json").read_text())
 indexes={(m["split_id"],m["method"]):m for m in mets}
 expected={(f"rep00_fold{i:02d}",arm) for i in range(3) for arm in ARMS}
 check(len(mets)==len(indexes)==12 and set(indexes)==expected,
       "method by fold metrics not complete")
 grouped=defaultdict(list)
 for row in patients:
  key=(row["split_id"],row["method"])
  check(key in expected,"unknown arm")
  grouped[key].append(row)
 check(set(grouped)==expected and len(patients)==4*N[ds],
       "E patient score inventory invalid")
 all_E=set()
 for f in range(3):
  fid=f"rep00_fold{f:02d}"
  ref={r["sample_id"]:r for r in grouped[fid,"ORIGINAL_RLR"]}
  check(len(ref)==len(grouped[fid,"ORIGINAL_RLR"]),
        "original E sample IDs duplicated")
  ids=set(ref)
  check(len(ids & all_E)==0,"outer E patient reused across folds")
  all_E|=ids
  check(fid in structures,"tree structure missing")
  struct=structures[fid]
  chosen=struct["chosen_type"]
  check(chosen in ("lookahead","static"),"unknown policy mode")
  if chosen=="lookahead":
   sub=struct["chosen_branch_second"]
   # Three explicitly named ternary branches must be represented by
   # distinct second pair tests for a genuinely conditional policy.
   if isinstance(sub,dict):
    check(set(sub)=={"-1","0","1"},"missing ternary branch")
    check(len(set(sub.values()))>=2,"not truly branched")
   else: check(len(set(sub))>=2,"not truly branched")
   check(float(struct["V_gain"])>=.25-1e-8,
         "validation acceptance gate bypassed")
  for arm in ARMS:
   rr=grouped[fid,arm]
   check(len(rr)==len(ids) and {r["sample_id"] for r in rr}==ids,
         "E score IDs mismatch")
   y=[];p=[];tested=[];active=[]
   for r in rr:
    source=ref[r["sample_id"]]
    truth=int(r["truth"]);full=float(source["prob"])
    approx=float(r["prob"]);original=float(r["original_prob"])
    lo=float(r["lower"]);hi=float(r["upper"])
    cnt=int(r["evaluated"]);act=int(r["active"])
    check(truth==int(source["truth"]) and truth in (0,1),
          "truth mismatch")
    check(math.isfinite(approx) and 0<=lo<=hi<=1 and
          math.isfinite(full) and 0<=full<=1,
          "certificate score invalid")
    check(abs(original-full)<1e-12,
          "frozen original RLR patient probability mismatch")
    check(lo-1e-8<=full<=hi+1e-8,
          "full RLR posterior outside certified interval")
    check(hi-lo<=.0500002 and
          abs(approx-full)<=.0250002,
          "posterior error bound violated")
    check((approx>=.5)==(full>=.5),
          "RLR classification changed")
    check(0<=act<=cnt<=50,
          "invalid evaluated/active relational test count")
    if arm=="CERT_TREE_DEPTH2_LOOKAHEAD":
     if chosen=="lookahead":
      check(r["route"] in ("-1","0","1"),
            "tree branch route not recorded")
      check(r["path"]!="",
            "patient-specific tree path missing")
     else:
      check(r["route"]=="STATIC",
            "unaccepted tree still routing conditionally")
    y.append(truth);p.append(approx)
    tested.append(cnt);active.append(act)
   independent=metric(y,p)
   met=indexes[fid,arm]
   for col,v in independent.items():
    check(abs(v-float(met[col]))<1e-8,
          f"independent patient metric mismatch {fid}/{arm}/{col}")
   check(abs(sum(tested)/len(tested)-
              float(met["mean_evaluated"]))<1e-8,
         "evaluated count mismatch")
   check(abs(sum(active)/len(active)-
              float(met["mean_active"]))<1e-8,
         "active evidence count mismatch")
   check(met["training_provenance"]==
      "ORIGINAL_MANUSCRIPT_PUBLIC_SOURCE_NOT_ARFF_PARITY",
      "incorrect ARFF reproduction claim")
 check(len(all_E)==N[ds],"3 outer test folds do not cover cohort")
 print("CL2_INDEPENDENT_DATASET="+ds)
 for arm in ARMS:
  results=[indexes[f"rep00_fold{i:02d}",arm] for i in range(3)]
  print("CL2_INDEPENDENT_ARM="+arm,
        "BA="+format(sum(float(a["BA"]) for a in results)/3,".5f"),
        "macroF1="+format(sum(float(a["macro_F1"]) for a in results)/3,".5f"),
        "LL="+format(sum(float(a["logloss"]) for a in results)/3,".5f"),
        "mean_tests="+format(sum(float(a["mean_evaluated"]) for a in results)/3,".3f"))
 print("CL2_LOOKAHEAD_ACCEPTED_FOLDS="+str(sum(
     structures[f"rep00_fold{i:02d}"]["chosen_type"]=="lookahead" for i in range(3))))
 print("CL2_INDEPENDENT_ALL_CLASS_AND_POSTERIOR_CERTIFICATES=PASS")
 print("CL2_INDEPENDENT_3_FOLDS_COVER_PATIENTS_ONCE=PASS")
 print("CL2_DATASET_ORIGINAL_PROCESSED_ARFF_PARITY=NOT_VERIFIED")
if __name__=="__main__":main()
