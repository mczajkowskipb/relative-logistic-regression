#!/usr/bin/env python3
"""Stage D: paired sign-flip *sensitivity* simulation for six frozen V matrices.

Centered observed per-patient candidate gains define an artificial
zero-mean null. One Rademacher multiplier per patient is shared across
six candidates. Repeated 5000 times/fold. Descriptive selection-risk
stress test ONLY, not a calibrated FWER/p-value for original CRIT.
"""
from __future__ import annotations
import csv
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist, mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v8_stageD_null_simulation_v1"
REPS=5000
SEED=20261025
COST=0.002
MIN_GAIN=0.002
BONF_Z=NormalDist().inv_cdf(1-0.05/6)
DATASETS=("chen-2002","gordon-2002")
FOLDS=("rep00_fold00","rep00_fold01","rep00_fold02")

def need(ok,msg):
    if not ok:raise SystemExit("STAGED_SIM_STOP: "+msg)

def read_csv(file):
    need(file.is_file(),"missing "+str(file))
    with file.open(newline="",encoding="utf-8") as f:
        vals=list(csv.DictReader(f))
    need(bool(vals),"empty "+str(file))
    return vals

def val(s):
    f=float(s)
    need(math.isfinite(f),"non-finite number")
    return f

def statistic(a,z=1.):
    n=len(a)
    need(n>=3,"validation sample insufficient")
    avg=sum(a)/n
    variance=sum((q-avg)**2 for q in a)/(n-1)
    return avg-z*math.sqrt(variance/n)-COST

def run():
    gains=read_csv(OUT/"candidate_V_gains.csv")
    manifest=read_csv(OUT/"candidate_score_manifest.csv")
    need(len(manifest)==36,"missing archived six candidate scores per fold")
    fold_candidates=defaultdict(dict)
    for x in gains:
        key=(x["dataset"],x["split_id"])
        rank=int(x["D_rank"])
        need(rank in range(1,7),"invalid candidate rank")
        fold_candidates[key].setdefault(rank,[]).append(x)
    scores={(r["dataset"],r["split_id"],int(r["D_rank"])):r
            for r in manifest}
    need(len(scores)==36,"candidate score duplicates")
    required={(d,f) for d in DATASETS for f in FOLDS}
    need(set(fold_candidates)==required,"incorrect fold inventory")

    detail=[]
    summ=[]
    for ix,(ds,fold) in enumerate(sorted(required)):
        candidates=fold_candidates[ds,fold]
        need(set(candidates)==set(range(1,7)),"incomplete candidate matrix")
        sample_ids=[r["sample_id"] for r in candidates[1]]
        n=len(sample_ids)
        need(len(set(sample_ids))==n,"V sample ids duplicated")
        arr=[]
        for k in range(1,7):
            records=candidates[k]
            need([r["sample_id"] for r in records]==sample_ids,
                 "candidate V patient identity/order differs")
            col=[val(r["patient_gain"]) for r in records]
            observed=statistic(col)
            archived=val(scores[ds,fold,k]["confirmation_lower"])
            need(abs(observed-archived)<1e-7,
                 "recomputed observed V score differs from RDS")
            arr.append(col)
        means=[mean(col) for col in arr]
        # One paired sign per patient keeps cross-candidate correlation.
        centered=[[arr[k][i]-means[k] for k in range(6)] for i in range(n)]
        rng=random.Random(SEED+ix*1009)
        count1=countmax=countbonf=0
        selection_premium=[]
        observed_lower=[statistic(col) for col in arr]
        for rep in range(REPS):
            signs=[1 if rng.getrandbits(1) else -1 for _ in range(n)]
            samples=[
                [signs[i]*centered[i][k] for i in range(n)]
                for k in range(6)
            ]
            lowers=[statistic(col) for col in samples]
            bonf=[statistic(col,BONF_Z) for col in samples]
            first=lowers[0]
            best=max(lowers)
            best_bonf=max(bonf)
            a1=int(first>MIN_GAIN)
            am=int(best>MIN_GAIN)
            ab=int(best_bonf>MIN_GAIN)
            count1+=a1
            countmax+=am
            countbonf+=ab
            selection_premium.append(best-first)
            detail.append(dict(
                dataset=ds,split_id=fold,replicate=rep,
                first_D_lower=round(first,12),
                best_of_six_V_lower=round(best,12),
                max_six_bonf_lower=round(best_bonf,12),
                first_accepted=a1,
                best_of_six_accepted=am,
                bonf_accepted=ab))
        observed_first=observed_lower[0]
        observed_best=max(observed_lower)
        observed_best_rank=1+observed_lower.index(observed_best)
        summ.append(dict(
            dataset=ds,split_id=fold,n_V=n,n_candidates=6,
            repetitions=REPS,seed=SEED+ix*1009,
            observed_first_lower=round(observed_first,12),
            observed_best_lower=round(observed_best,12),
            observed_best_D_rank=observed_best_rank,
            null_first_accept_rate=round(count1/REPS,8),
            null_best_six_accept_rate=round(countmax/REPS,8),
            null_bonf_six_accept_rate=round(countbonf/REPS,8),
            null_mean_selection_premium=round(mean(selection_premium),10),
            bonferroni_z=round(BONF_Z,9)))
        print(f"STAGED_SIM {ds} {fold} nV={n} observedVbestRank={observed_best_rank} "
              f"paired-null accept: rank1={count1/REPS:.4f} "
              f"max6={countmax/REPS:.4f} "
              f"bonf6={countbonf/REPS:.4f} "
              f"mean_null_selection_premium={mean(selection_premium):.6f}",flush=True)
    with (OUT/"null_draws.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=detail[0].keys())
        w.writeheader()
        w.writerows(detail)
    with (OUT/"selection_null_summary.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=summ[0].keys())
        w.writeheader()
        w.writerows(summ)
    need(len(detail)==6*REPS and len(summ)==6,"incomplete simulation")
    print("STAGED_NULL_30000_PAIRED_SIMULATION_DRAWS=PASS")
    print("STAGED_ZERO_MEAN_SENSITIVITY_ONLY_NOT_CALIBRATED_FWER")
    print("STAGED_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")

if __name__=="__main__":
    run()
