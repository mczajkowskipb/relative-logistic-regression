#!/usr/bin/env python3
"""Read-only Stage D patient gain, archived V score and Monte Carlo inventory auditor."""
from __future__ import annotations
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v8_stageD_null_simulation_v1"
DATASETS=("chen-2002","gordon-2002")
FOLDS=("rep00_fold00","rep00_fold01","rep00_fold02")
METHODS=range(1,7)

def check(ok,msg):
    if not ok: raise SystemExit("STAGED_AUDIT_STOP: "+msg)

def read(path):
    check(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        a=list(csv.DictReader(f))
    check(bool(a),"empty "+str(path))
    return a

def num(x):
    f=float(x)
    check(math.isfinite(f),"nonfinite")
    return f

def main():
    check(len(sys.argv)==3,"usage audit.py CRIT7_WORKTREE ORIGINAL_REF")
    v7=Path(sys.argv[1]).resolve(strict=True)
    v7manifest=read(v7/"results/crit_v7_original_rlr_2x3_v1"/"run_manifest.csv")
    check(len(v7manifest)==24,"v7 prior manifest incomplete")
    original_ref=Path(sys.argv[2]).resolve(strict=True)
    original_archive=original_ref/"results/rlmt_dc_published_data_5fold_v1"
    original_manifest=read(original_archive/"run_manifest.csv")
    original_keys={(r["dataset"],r["split_id"],r["method"]):r for r in original_manifest}
    check(len(original_manifest)==285,"original archive inventory incomplete")
    gains=read(OUT/"candidate_V_gains.csv")
    scores=read(OUT/"candidate_score_manifest.csv")
    sim=read(OUT/"null_draws.csv")
    summary=read(OUT/"selection_null_summary.csv")
    check(len(scores)==36 and len(sim)==30000 and len(summary)==6,
          "unexpected simulation file dimensions")
    expected={(ds,f) for ds in DATASETS for f in FOLDS}
    groups=defaultdict(dict)
    for r in gains:
        key=(r["dataset"],r["split_id"])
        rank=int(r["D_rank"])
        check(key in expected and rank in METHODS,
              "unknown candidate V record")
        groups[key].setdefault(rank,[]).append(r)
    check(set(groups)==expected,"missing fold V vectors")
    scores_by={(r["dataset"],r["split_id"],int(r["D_rank"])):r
               for r in scores}
    check(len(scores_by)==36,"score duplicates")
    summary_by={(r["dataset"],r["split_id"]):r for r in summary}
    check(len(summary_by)==6 and set(summary_by)==expected,
          "missing simulation summary")
    previous={(r["dataset"],r["split_id"],r["method"]):r
              for r in v7manifest}
    sim_groups=defaultdict(list)
    for x in sim:
        key=(x["dataset"],x["split_id"])
        check(key in expected,"unknown null replication fold")
        sim_groups[key].append(x)
    check(set(sim_groups)==expected,"missing null simulation fold")
    for ds,sid in sorted(expected):
        candidate=groups[ds,sid]
        check(set(candidate)==set(METHODS),"not all 6 candidate ranks")
        ref_ids=[r["sample_id"] for r in candidate[1]]
        n=len(ref_ids)
        check(n>=20 and len(set(ref_ids))==n,"invalid V patient support")
        orig=original_keys[ds,sid,"RLR"]
        olddir=(original_ref/orig["model_dir"]).resolve(strict=True)
        check(olddir.is_relative_to(original_archive.resolve(strict=True)),
              "archived model path traversal")
        class_rows=read(olddir.parent/"class_support.csv")
        check(len(class_rows)==2,"original factor class levels not found")
        positive=class_rows[1]["class"]
        lowers=[]
        for rank in METHODS:
            rr=candidate[rank]
            check(len(rr)==n and [r["sample_id"] for r in rr]==ref_ids,
                  "V patient mismatch across candidates")
            losses=[]
            for row in rr:
                y=row["y_true"]
                p0=num(row["baseline_prob"])
                p1=num(row["candidate_prob"])
                check(0<=p0<=1 and 0<=p1<=1,
                      "candidate probability outside [0,1]")
                gain=num(row["patient_gain"])
                y=int(row["y_true"]==positive)
                p0=max(1e-15,min(1-1e-15,p0))
                p1=max(1e-15,min(1-1e-15,p1))
                recomputed=-(y*math.log(p0)+(1-y)*math.log1p(-p0))+(
                            y*math.log(p1)+(1-y)*math.log1p(-p1))
                check(abs(gain-recomputed)<1e-9,
                      "saved patient logloss gain incorrect")
                losses.append(gain)
            avg=sum(losses)/n
            se=math.sqrt(sum((g-avg)**2 for g in losses)/(n-1)/n)
            lower=avg-se-.002
            s=scores_by[ds,sid,rank]
            check(int(s["n_V"])==n and
                  abs(avg-num(s["confirmation_mean"]))<1e-8 and
                  abs(se-num(s["confirmation_se"]))<1e-8 and
                  abs(lower-num(s["confirmation_lower"]))<1e-8 and
                  num(s["saved_score_parity_max_abs"])<1e-7,
                  "recomputed V candidate score differs from archive")
            lowers.append(lower)
        prior=previous[ds,sid,"CRIT_TREE"]
        check(int(prior["tree_accepted"]) ==
              int(max(lowers)>.002),
              "v7 archival accept policy differs")
        s=summary_by[ds,sid]
        check(int(s["n_V"])==n and int(s["n_candidates"])==6 and
              int(s["repetitions"])==5000 and
              int(s["observed_best_D_rank"])==1+max(range(6),key=lambda i:lowers[i]) and
              abs(num(s["observed_first_lower"])-lowers[0])<1e-8 and
              abs(num(s["observed_best_lower"])-max(lowers))<1e-8,
              "archived simulation summary initial scores disagree")
        draws=sorted(sim_groups[ds,sid],key=lambda x:int(x["replicate"]))
        check(len(draws)==5000 and
              [int(x["replicate"]) for x in draws]==list(range(5000)),
              "missing or duplicated pseudo-null trials")
        sums=[0,0,0]
        premiums=[]
        for draw in draws:
            first=num(draw["first_D_lower"])
            best=num(draw["best_of_six_V_lower"])
            bonf=num(draw["max_six_bonf_lower"])
            flags=[int(draw[k]) for k in
                   ("first_accepted","best_of_six_accepted","bonf_accepted")]
            check(all(i in (0,1) for i in flags) and
                  flags==[int(first>.002),int(best>.002),int(bonf>.002)] and
                  best>=first-1e-10 and bonf<=best+1e-10,
                  "inconsistent null replication statistic")
            sums=[sums[i]+flags[i] for i in range(3)]
            premiums.append(best-first)
        names=("null_first_accept_rate","null_best_six_accept_rate",
               "null_bonf_six_accept_rate")
        check(all(abs(sums[i]/5000-num(s[names[i]]))<1e-7 for i in range(3)) and
              abs(sum(premiums)/5000 -
                  num(s["null_mean_selection_premium"]))<1e-7,
              "stored Monte Carlo summary not reproducible from draws")
        print(f"STAGED_AUDIT {ds}:{sid}: nV={n} "
              f"first_accept={sums[0]/5000:.4f} "
              f"max_accept={sums[1]/5000:.4f} "
              f"bonf_accept={sums[2]/5000:.4f}")
    print("STAGED_ARCHIVED_V_GAIN_METRIC_PARITY_AUDIT=PASS")
    print("STAGED_30000_NULL_DRAW_SUMMARY_AUDIT=PASS")
    print("STAGED_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")

if __name__=="__main__":
    main()
