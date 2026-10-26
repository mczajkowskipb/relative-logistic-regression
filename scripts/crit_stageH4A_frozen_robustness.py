#!/usr/bin/env python3
"""Stage H4A: predeclared independent synthetic stress test of FROZEN H3 LSH.

Crucial: import and call the unchanged H3 blind_rank_pair_lsh algorithm.
Only its global N is set to the current D size; HASH_BITS=18 and HASH_TABLES=24
remain unchanged. No oracle anchor, original-RLR fit, V or heldout E.
Failure to return six results is an OBSERVED metric, not a program failure.
"""
from __future__ import annotations

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import crit_stageH3_blind_lsh_pilot as frozen

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_stageH4A_blind_robustness_v1"

N_SIZES=(60,120,240)
DIMENSIONS=(32,128,512)
AMPLITUDES=(0.0,1.0,1.6,2.4)
REPS=8
START_SEED=20280617
ALGORITHM_TAG="FROZEN_H3_LSH_18BITS_24TABLES"
PRIMARY="TRUE_PAIR_PRODUCT_IN_TOP6_D_ONLY"
N_CASES=len(N_SIZES)*len(DIMENSIONS)*len(AMPLITUDES)*REPS

def need(ok,msg):
    if not ok:
        raise RuntimeError("H4A_STOP: "+msg)

def case_seed(n,p,amp_index,rep):
    return START_SEED + n*101 + p*313 + amp_index*100003 + rep*47

def make_case(n,p,amp,seed):
    rng=np.random.default_rng(seed)
    X=rng.standard_normal((n,p)).astype(np.float32)
    z=np.where(X[:,0]>X[:,1],1.0,-1.0)
    r=np.where(X[:,2]>X[:,3],1.0,-1.0)
    eta=amp*z*r
    prob=1/(1+np.exp(-eta))
    y=(rng.random(n)<prob).astype(np.uint8)
    return X,y

def write_csv(path,rows,cols):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

def main():
    need(OUT.is_dir() and (OUT/"SOURCE_SHA.txt").is_file()
         and not (OUT/"study_summary.csv").exists(),
         "new pinned source archive required; refuse overwrite")
    need((frozen.HASH_BITS,frozen.HASH_TABLES,frozen.SCORE_BATCH)==
         (18,24,4096),"frozen H3 algorithm settings changed")
    summaries=[]
    top6=[]
    for n in N_SIZES:
      for p in DIMENSIONS:
        for ai,amp in enumerate(AMPLITUDES):
          for rep in range(1,REPS+1):
            seed=case_seed(n,p,ai,rep)
            X,y=make_case(n,p,amp,seed)
            need(len(np.unique(y))==2,
                 f"unsupported one-class sample seed {seed}")
            frozen.N=n
            res=frozen.blind_rank_pair_lsh(X,y,seed)
            no=len(res["top"])
            need(0<=no<=6 and
                 res["candidate_pairs_scored"]>=0,
                 "frozen H3 returned impossible shortlist")
            collision=int(res["true_joint_hash_collision"])
            hit=int(res["true_joint_in_top6"])
            need(hit<=collision and (hit==1)==(res["true_joint_top6_rank"]>0),
                 "true pair hash collision / shortlist contradiction")
            row=dict(
              n_D=n,p=p,amplitude=amp,amplitude_index=ai,
              replicate=rep,seed=seed,
              label_positive=int(y.sum()),
              method=ALGORITHM_TAG,
              endpoint=PRIMARY,
              hash_bits=frozen.HASH_BITS,
              hash_tables=frozen.HASH_TABLES,
              pair_feature_count=int(res["pair_feature_count"]),
              candidate_pairs_proposed=int(res["candidate_pairs_proposed"]),
              candidate_pairs_scored=int(res["candidate_pairs_scored"]),
              shortlist_size=no,
              true_joint_hash_collision=collision,
              true_joint_in_top6=hit,
              true_joint_top6_rank=int(res["true_joint_top6_rank"]),
              best_score=float(res["best_score"]),
              wall_seconds=float(res["seconds"]),
              process_peak_rss_kb=int(res["peak_rss_kb"]),
              scientific_stage="INDEPENDENT_SIMULATION_DEVELOPMENT_ONLY"
            )
            summaries.append(row)
            for c in res["top"]:
                top6.append(dict(
                   n_D=n,p=p,amplitude_index=ai,replicate=rep,
                   seed=seed,rank=c["rank"],pairA_a=c["pairA_a"],
                   pairA_b=c["pairA_b"],pairB_a=c["pairB_a"],
                   pairB_b=c["pairB_b"],score=c["score"],
                   true_pair_product=c["true_pair_product"]
                ))
            print(f"H4A_CASE n={n} p={p} signal={amp:.1f} rep={rep} "
                  f"collision={collision} recovered={hit} "
                  f"top_n={no} proposed={res['candidate_pairs_scored']} "
                  f"elapsed_s={res['seconds']:.3f}",flush=True)
    need(len(summaries)==N_CASES and len(set(
         (x["n_D"],x["p"],x["amplitude_index"],x["replicate"])
         for x in summaries))==N_CASES,
         "missing/duplicate 288 study IDs")
    write_csv(OUT/"study_summary.csv",summaries,list(summaries[0].keys()))
    need(all(0<=x["shortlist_size"]<=6 for x in summaries) and
         len(top6)==sum(x["shortlist_size"] for x in summaries),
         "shortlist row inventory inconsistent")
    write_csv(OUT/"top6_candidates.csv",top6,
      ("n_D","p","amplitude_index","replicate","seed","rank",
       "pairA_a","pairA_b","pairB_a","pairB_b","score","true_pair_product"))
    for amp in AMPLITUDES:
      print(f"=== H4A SIGNAL {amp:.1f} ===")
      for n in N_SIZES:
       for p in DIMENSIONS:
        rr=[x for x in summaries if x["amplitude"]==amp
            and x["n_D"]==n and x["p"]==p]
        hits=sum(x["true_joint_in_top6"] for x in rr)
        collisions=sum(x["true_joint_hash_collision"] for x in rr)
        means=float(np.mean([x["wall_seconds"] for x in rr]))
        print(f"H4A_CELL signal={amp:.1f} n={n} p={p} "
              f"hit={hits}/8 collision={collisions}/8 "
              f"mean_elapsed_s={means:.3f}")
    print("STAGEH4A_288_FROZEN_LSH_STUDIES=PASS")
    print("STAGEH4A_36_PREDECLARED_CELLS=PASS")
    print("STAGEH4A_TRUE_RANK_PAIR_RESULTS_ARCHIVED=PASS")
    print("STAGEH4A_SCIENTIFIC_STAGE=SIMULATION_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
