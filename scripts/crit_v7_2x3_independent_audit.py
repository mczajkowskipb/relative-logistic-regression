#!/usr/bin/env python3
"""CRIT-v7 read-only heldout patient-level independent pilot audit."""
import argparse
import csv
import math
from pathlib import Path
from statistics import mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v7_original_rlr_2x3_v1"
OLD="results/rlmt_dc_published_data_5fold_v1"
DS=("chen-2002","gordon-2002")
FOLDS=("rep00_fold00","rep00_fold01","rep00_fold02")
METHODS=("ORIGINAL_RLR","GLOBAL_CAL",
         "GLOBAL_INTERACTION","CRIT_TREE")
def need(ok,msg):
    if not ok:raise SystemExit("CRIT7_AUDIT_STOP: "+msg)
def csvrows(p):
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    need(rows,"empty "+str(p))
    return rows
def num(s):
    try:v=float(s)
    except (TypeError,ValueError) as e:
        raise SystemExit("CRIT7_AUDIT_STOP: nonnumeric") from e
    need(math.isfinite(v),"nonfinite")
    return v
def ll(y,p,pos):
    need(len(y)==len(p)>0 and all(0<=z<=1 for z in p),
         "invalid probability")
    yy=[int(t==pos) for t in y]
    need(set(yy)=={0,1},"both classes required for test")
    z=[max(1e-15,min(1-1e-15,a)) for a in p]
    return (-sum(t*math.log(q)+(1-t)*math.log1p(-q)
             for t,q in zip(yy,z))/len(y),
            sum((t-q)**2 for t,q in zip(yy,p))/len(y))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-ref",required=True)
    args=ap.parse_args()
    ref=Path(args.original_ref).resolve(strict=True)
    old=(ref/OLD).resolve(strict=True)
    need(Path.cwd().resolve()==ROOT,"run inside CRIT v7 worktree")
    need((old/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "old original 285 records incomplete")
    m=csvrows(OUT/"run_manifest.csv")
    prior=csvrows(old/"run_manifest.csv")
    need(len(m)==24 and len(prior)==285,
         "incorrect pilot and original model inventory")
    k={(r["dataset"],r["split_id"],r["method"]):r for r in m}
    o={(r["dataset"],r["split_id"],r["method"]):r for r in prior}
    need(len(k)==24,"duplicate method-fold keys")
    need(set(k)=={(d,s,method) for d in DS for s in FOLDS
                  for method in METHODS},
         "not exact 2x3x4 frozen design")
    scores={method:[] for method in METHODS}
    for ds in DS:
      for sid in FOLDS:
        baseline=o[(ds,sid,"RLR")]
        olddir=(ref/baseline["model_dir"]).resolve(strict=True)
        need(olddir.is_relative_to(old),
             "original frozen prediction path escapes archive")
        original=csvrows(olddir/"predictions_test.csv")
        classes=csvrows(olddir.parent/"class_support.csv")
        need(len(classes)==2,"positive label reference unavailable")
        positive=classes[1]["class"]
        ids=[(r["sample_id"],r["y_true"]) for r in original]
        need(len(ids)==len(set(ids)),"duplicate original patient ids")
        y=[r["y_true"] for r in original]
        origprob=[num(r["prob"]) for r in original]
        origll,_=ll(y,origprob,positive)
        for method in METHODS:
          row=k[(ds,sid,method)]
          need(row["status"]=="OK" and
               row["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEV_ONLY",
               "model status or honesty marker invalid")
          p=(ROOT/row["model_dir"]).resolve(strict=True)
          need(p.is_relative_to(OUT.resolve(strict=True)),
               "new model escaped output root")
          saved=csvrows(p/"predictions_test.csv")
          need([(r["sample_id"],r["y_true"]) for r in saved]==ids,
               "patient/sample label alignment mismatch")
          probs=[num(r["prob"]) for r in saved]
          root=[num(r["original_RLR_prob"]) for r in saved]
          need(max(abs(x-z) for x,z in zip(root,origprob))<1e-9,
               "CRIT7 ORIGINAL_RLR parent not original paper")
          v,b=ll(y,probs,positive)
          need(abs(v-num(row["logloss"]))<1e-9 and
               abs(b-num(row["brier"]))<1e-9 and
               abs(origll-num(row["original_logloss"]))<1e-9 and
               abs((v-origll)-num(row["delta_vs_original"]))<1e-9,
               "saved metrics disagree with patient predictions")
          need(int(row["n_test"])==len(saved),
               "unexpected test sample count")
          n=int(row["tree_accepted"])
          need(n in (0,1),"actual two-leaf tree status invalid")
          if method!="CRIT_TREE":
              need(n==0,"control must never report tree")
          if n:
              names=[row[z] for z in (
                  "gate_gene_a","gate_gene_b",
                  "leaf_false_gene_a","leaf_false_gene_b",
                  "leaf_true_gene_a","leaf_true_gene_b")]
              need(all(names) and
                   names[0]!=names[1] and
                   names[2]!=names[3] and
                   names[4]!=names[5],
                   "accepted genuine tree lacks two named leaf experts")
          leaf=[int(a["leaf_id"]) for a in saved]
          need(all(x in (0,1,2) for x in leaf),"invalid saved leaf ID")
          if method=="CRIT_TREE":
              need((all(x==0 for x in leaf) if n==0 else
                     all(x in (1,2) for x in leaf)),
                   "saved true leaf route does not match split state")
          else:
              need(all(x==0 for x in leaf),
                   "non-tree control stores leaf decision")
          if method=="CRIT_TREE" and n==0:
              need(max(abs(q-z) for q,z in zip(probs,origprob))<1e-12,
                  "pruned tree must be exactly original paper RLR")
          scores[method].append((ds,sid,v-origll,n))
    need(all(len(v)==6 for v in scores.values()),
         "all fold-method predictions must be checked")
    for method in METHODS:
        rows=scores[method]
        wins=sum(z[2]<-1e-10 for z in rows)
        losses=sum(z[2]>1e-10 for z in rows)
        print(f"{method}: accepted_tree_folds={sum(x[3] for x in rows)}/6 "
              f"foldW/L/T_vs_original={wins}/{losses}/{6-wins-losses} "
              f"mean_deltaLL_original={mean(x[2] for x in rows):+.9f}")
    print("CRIT7_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS")
    print("CRIT7_30_METHOD_FOLD_PREDICTION_VECTORS_AUDITED=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")
if __name__=="__main__":main()
