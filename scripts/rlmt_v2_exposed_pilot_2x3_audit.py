#!/usr/bin/env python3
"""Independent audit for EXACT 2×3×2 exposed-paper H1/H2 pilot results.

No refitting and no editing of archived models. Checks that old original
RLR sample predictions, new aligned parent and both v2 methods are paired,
reproduces logloss/Brier and enforces exact no-split fallback.
"""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
FOLDER="results/rlmt_v2_exposed_pilot_2x3_v1"
OLD="results/rlmt_dc_published_data_5fold_v1"
DATASETS=("chen-2002","gordon-2002")
SPLITS=tuple(f"rep00_fold{i:02d}" for i in range(3))
MODES=("H1","H2")
TOL=1e-9

def need(ok,message):
    if not ok:raise SystemExit("V2_PILOT_AUDIT_STOP: "+message)

def readcsv(p):
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def value(x):
    try: v=float(x)
    except (ValueError,TypeError):
        raise SystemExit("V2_PILOT_AUDIT_STOP: nonnumeric probability")
    need(math.isfinite(v),"nonfinite probability")
    return v

def score(ys,ps,positive):
    need(len(ys)==len(ps)>0,"missing predictions")
    vals=[int(y==positive) for y in ys]
    need(set(vals)=={0,1},"test fold missing class")
    need(all(0<=p<=1 for p in ps),"probability out of bounds")
    def clipped(p):return min(1-1e-15,max(1e-15,p))
    ll=-sum(y*math.log(clipped(p))+
        (1-y)*math.log1p(-clipped(p))
        for y,p in zip(vals,ps))/len(ys)
    brier=sum((p-y)**2 for p,y in zip(ps,vals))/len(vals)
    return ll,brier

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--reference-repo",required=True)
    args=parser.parse_args()
    need(ROOT==Path.cwd().resolve(),"must run from pilot worktree")
    ref=Path(args.reference_repo).resolve(strict=True)
    pilot=(ROOT/FOLDER).resolve(strict=True)
    old=(ref/OLD).resolve(strict=True)
    entries=readcsv(pilot/"run_manifest.csv")
    oldentries=readcsv(old/"run_manifest.csv")
    need(len(entries)==12 and len(oldentries)==285,
         "expected 12 pilot fits and 285 archival fit records")
    orig={(x["world"],x["dataset"],x["split_id"],x["method"]):x
          for x in oldentries}
    need(len(orig)==285,"duplicates in archival fit index")
    values={}
    summary=defaultdict(list)
    for e in entries:
        world,dataset,sid,mode=(e[k] for k in ("world","dataset","split_id","mode"))
        need(dataset in DATASETS and sid in SPLITS and mode in MODES,
             "unexpected cohort/split/method")
        k=dataset,sid,mode
        need(k not in values,"duplicate method fold")
        values[k]=e
        need(e["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY" and
             e["status"]=="OK","incomplete fit/evidence status")
        model=(ROOT/e["model_dir"]).resolve(strict=True)
        need(model.is_relative_to(pilot),"model path outside new result root")
        arch=orig.get((world,dataset,sid,"RLR"))
        need(arch is not None,"missing exact historical RLR comparator")
        archival=(ref/arch["model_dir"]).resolve(strict=True)
        need(archival.is_relative_to(old),"archived model path outside old root")
        obs=readcsv(archival/"predictions_test.csv")
        new=readcsv(model/"predictions_test.csv")
        need(len(obs)==len(new)==int(e["n_test"]),"paired sample count mismatch")
        need([(r["sample_id"],r["y_true"]) for r in obs] ==
             [(r["sample_id"],r["y_true"]) for r in new],
             "heldout sample ID or label mismatch")
        support=readcsv(archival.parent/"class_support.csv")
        need(len(support)==2,"original classes missing")
        positive=support[1]["class"]
        yy=[r["y_true"] for r in new]
        p=[value(r["prob"]) for r in new]
        rootp=[value(r["parent_prob"]) for r in new]
        origp=[value(r["prob"]) for r in obs]
        got=score(yy,p,positive)
        parent=score(yy,rootp,positive)
        original=score(yy,origp,positive)
        for field,v in (
            ("new_logloss",got[0]),("new_brier",got[1]),
            ("aligned_parent_logloss",parent[0]),
            ("aligned_brier",parent[1]),
            ("original_RLR_logloss",original[0]),
            ("original_RLR_brier",original[1]),
            ("delta_vs_aligned",got[0]-parent[0]),
            ("delta_vs_original",got[0]-original[0])):
            need(abs(value(e[field])-v)<TOL,
                 "saved metric differs "+str((k,field)))
        accepted=e["accepted"].lower()=="true"
        need(e["accepted"].lower() in ("true","false"),
             "invalid acceptance Boolean")
        if not accepted:
            need(p==rootp,"no-split not exactly aligned parent")
        summary[mode].append(
            (dataset,sid,accepted,got[0]-parent[0],
             got[0]-original[0]))
    need(set(values)=={
      (dataset,split,mode)
      for dataset in DATASETS for split in SPLITS for mode in MODES
    },"not the exact 2×3×2 predeclared study")
    for dataset in DATASETS:
        for sid in SPLITS:
            h=readcsv(ROOT/values[(dataset,sid,"H1")]["model_dir"]/
                      "predictions_test.csv")
            o=readcsv(ROOT/values[(dataset,sid,"H2")]["model_dir"]/
                      "predictions_test.csv")
            need([x["parent_prob"] for x in h] ==
                 [x["parent_prob"] for x in o],
                 "H1 and H2 aligned parent differ; model comparison invalid")
    print("RLMT_V2_PILOT_AUDITED_FITS=12")
    print("RLMT_V2_PILOT_PAIRED_FOLDS=6")
    for mode in MODES:
        group=summary[mode]
        wins=sum(row[3]<-1e-10 for row in group)
        losses=sum(row[3]>1e-10 for row in group)
        print(f"{mode}: accepted={sum(x[2] for x in group)}/6 "
              f"vs_aligned_win/loss/tie="
              f"{wins}/{losses}/{6-wins-losses} "
              f"mean_delta_logloss_vs_aligned="
              f"{sum(x[3] for x in group)/6:+.7f}")
    print("RLR_ALIGNED_DISTINCT_FROM_ARCHIVED_RLR=TRUE")
    print("RLMT_V2_PILOT_INDEPENDENT_AUDIT=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY")
if __name__=="__main__": main()
