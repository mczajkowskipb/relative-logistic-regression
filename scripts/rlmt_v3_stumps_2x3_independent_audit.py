#!/usr/bin/env python3
"""Independent read-only 18-fold-method RLR v3 pilot prediction audit."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OUT="results/rlmt_v3_stumps_2x3_pilot_v1"
OLD="results/rlmt_dc_published_data_5fold_v1"
DS=("chen-2002","gordon-2002")
FOLD=tuple(f"rep00_fold{i:02d}" for i in range(3))
MODES=("STUMP_TSP_BETA","STUMP_RIDGE_OFFSET","TREE_D2_PRUNED_BETA")

def need(cond,msg):
    if not cond:raise SystemExit("V3_PILOT_AUDIT_STOP: "+msg)

def read(p):
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:
        r=list(csv.DictReader(f))
    need(bool(r),"empty "+str(p))
    return r

def num(v):
    try:v=float(v)
    except (TypeError,ValueError) as e:
        raise SystemExit("V3_PILOT_AUDIT_STOP: nonnumeric") from e
    need(math.isfinite(v),"nonfinite")
    return v

def scores(y,probs,positive):
    need(len(y)==len(probs)>0 and
         all(0<=p<=1 for p in probs),
         "invalid probabilities")
    binary=[int(z==positive) for z in y]
    need(set(binary)=={0,1},"one-class test fold")
    q=[max(1e-15,min(1-1e-15,p)) for p in probs]
    return (
      -sum(a*math.log(p)+(1-a)*math.log1p(-p)
            for a,p in zip(binary,q))/len(y),
      sum((p-a)**2 for p,a in zip(probs,binary))/len(y)
    )

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--original-ref",required=True)
    v=a.parse_args()
    need(Path.cwd().resolve()==ROOT,"run from v3 pilot worktree")
    orig=Path(v.original_ref).resolve(strict=True)
    out=(ROOT/OUT).resolve(strict=True)
    old=(orig/OLD).resolve(strict=True)
    need((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "original source incomplete")
    m=read(out/"run_manifest.csv")
    oldm=read(old/"run_manifest.csv")
    need(len(m)==18 and len(oldm)==285,"unexpected record counts")
    original={(r["world"],r["dataset"],r["split_id"],r["method"]):r
              for r in oldm}
    entries={(r["dataset"],r["split_id"],r["method"]):r for r in m}
    need(len(entries)==18,"duplicate v3 rows")
    need(set(entries)=={(d,f,s) for d in DS for f in FOLD for s in MODES},
         "not predeclared 2x3x3 design")
    warnfile=out/"warning_inventory.csv"
    counts=defaultdict(int)
    if warnfile.exists():
        for r in read(warnfile):
            counts[(r["dataset"],r["split_id"],r["method"])]+=(int(r["frequency"]))
    summary=defaultdict(list)
    verified=0
    for dataset in DS:
      for sid in FOLD:
        parent=None
        expected=None
        for method in MODES:
            row=entries[(dataset,sid,method)]
            need(row["status"]=="OK" and
                 row["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEV_ONLY",
                 "model status/stage invalid")
            model=(ROOT/row["model_dir"]).resolve(strict=True)
            need(model.is_relative_to(out),"model outside new results root")
            pred=read(model/"predictions_test.csv")
            key=[(x["sample_id"],x["y_true"]) for x in pred]
            need(len(key)==len(set(key)),"duplicate test sample")
            if expected is None: expected=key
            need(key==expected,"different heldout sample grids")
            oldrow=original[(row["world"],dataset,sid,"RLR")]
            frozen=(orig/oldrow["model_dir"]).resolve(strict=True)
            need(frozen.is_relative_to(old),"old path escapes original")
            frozenpred=read(frozen/"predictions_test.csv")
            need([(x["sample_id"],x["y_true"]) for x in frozenpred]==key,
                 "frozen original RLR sample mismatch")
            support=read(frozen.parent/"class_support.csv")
            need(len(support)==2,"frozen classes absent")
            pos=support[1]["class"]
            truth=[x["y_true"] for x in pred]
            pp=[num(x["parent_prob"]) for x in pred]
            px=[num(x["prob"]) for x in pred]
            if parent is None:parent=pp
            need(len(parent)==len(pp) and all(
              abs(a-b)<1e-12 for a,b in zip(parent,pp)),
              "aligned parent differs among v3 models")
            ll,brier=scores(truth,px,pos)
            ll0,_=scores(truth,pp,pos)
            llold,_=scores(truth,[num(x["prob"]) for x in frozenpred],pos)
            for k,value in [
              ("logloss",ll),("brier",brier),
              ("aligned_logloss",ll0),("original_RLR_logloss",llold),
              ("delta_vs_aligned",ll-ll0),
              ("delta_vs_original",ll-llold)]:
                need(abs(num(row[k])-value)<1e-9,
                     "saved metric mismatch: "+str((dataset,sid,method,k)))
            need(int(row["warning_count"])==
                 counts[(dataset,sid,method)],
                 "warning inventory does not match count")
            nodes=int(row["realized_internal_nodes"])
            depth=int(row["chosen_depth"])
            need(nodes>=0 and nodes<=3,"invalid branch count")
            need(depth in (1,2),"invalid tree depth")
            if method=="STUMP_TSP_BETA":
                need(nodes in (0,1),"stump is not a stump")
            elif method=="STUMP_RIDGE_OFFSET":
                need(nodes in (0,1),"ridge-offset stump >1 node")
            else:
                need(nodes<=3 and depth in (1,2),
                     "depth2 pruned model invalid")
            summary[method].append((dataset,sid,ll-ll0,ll-llold,nodes))
            verified+=1
    need(verified==18,"not all 18 method-folds audited")
    for method in MODES:
        z=summary[method]
        gain=sum(x[2]<-1e-10 for x in z)
        harm=sum(x[2]>1e-10 for x in z)
        print(f"{method}: real_tree={sum(x[4]>0 for x in z)}/6 "
              f"LLvsAligned win/loss/tie={gain}/{harm}/{6-gain-harm} "
              f"meanDeltaAligned={sum(x[2] for x in z)/6:+.8f} "
              f"meanDeltaOriginal={sum(x[3] for x in z)/6:+.8f}")
    print("RLMT_V3_PILOT_INDEPENDENT_AUDIT=PASS")
    print("RLMT_V3_PILOT_AUDITED=18")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")

if __name__=="__main__":main()
