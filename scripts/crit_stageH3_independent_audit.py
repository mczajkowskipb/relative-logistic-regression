#!/usr/bin/env python3
"""Independent Stage H3 checksum-free replay of 240 saved top-six scores.

Regenerates every D synthetic sample from the pinned seed, recomputes
per-candidate signed pair-product correlation, checks rank/true pair IDs,
source inventory and reporting. Does not reproduce every LSH bucket join.
"""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

OUT=Path("results/crit_stageH3_blind_lsh_fix1_v1")
PS=(32,64,128,256,512)
SCENARIOS=("null","xor")
REPS=4
N=120

def need(x,msg):
    if not x: raise SystemExit("STAGEH3_AUDIT_STOP: "+msg)

def load(p):
    need(p.is_file(),"missing "+str(p))
    with p.open(newline="",encoding="utf-8") as f: return list(csv.DictReader(f))

def regenerate(p,seed,scenario):
    rng=np.random.default_rng(seed)
    X=rng.standard_normal((N,p)).astype(np.float32)
    z=np.where(X[:,0]>X[:,1],1.0,-1.0)
    r=np.where(X[:,2]>X[:,3],1.0,-1.0)
    eta=2.4*z*r if scenario=="xor" else np.zeros(N,dtype=float)
    y=(rng.random(N)<1/(1+np.exp(-eta))).astype(np.uint8)
    return X,y

def score(X,y,a,b,c,d):
    r1=(X[:,a-1]>X[:,b-1])
    r2=(X[:,c-1]>X[:,d-1])
    agreement=np.count_nonzero((r1^r2)!=y)
    return abs(2*agreement/N-1)

def main():
    summaries=load(OUT/"study_summary.csv")
    top=load(OUT/"top6_candidates.csv")
    expected={(s,p,r) for s in SCENARIOS for p in PS for r in range(1,REPS+1)}
    meta={(x["scenario"],int(x["p"]),int(x["replicate"])):x for x in summaries}
    need(len(meta)==len(expected)==40 and set(meta)==expected,
         "duplicated/missing study IDs")
    group=defaultdict(list)
    for x in top:
        k=(x["scenario"],int(x["p"]),int(x["replicate"]))
        need(k in expected,"unexpected shortlist ID")
        group[k].append(x)
    need(len(top)==240 and set(group)==expected,
         "top-six shortlist archive incomplete")
    for scenario,p,rep in sorted(expected):
        row=meta[scenario,p,rep]
        seed=20270211+(SCENARIOS.index(scenario)+1)*100003+p*317+rep*43
        X,y=regenerate(p,seed,scenario)
        M=p*(p-1)//2
        need(int(row["seed"])==seed and int(row["n_D"])==N and
             row["method"]=="BLIND_LSH_PAIR_PRODUCT" and
             row["scientific_stage"]=="D_ONLY_SYNTHETIC_DEVELOPMENT" and
             int(row["pair_feature_count"])==M and
             int(row["hash_bits"])==18 and
             int(row["hash_tables"])==24,
             "seed/hash or study design drift")
        n=int(row["candidate_pairs_scored"])
        need(n>0 and n==int(row["candidate_pairs_proposed"]) and
             math.isfinite(float(row["seconds"])) and
             float(row["seconds"])>=0 and
             float(row["peak_rss_kb"])>0,
             "invalid candidate budget/run metrics")
        cand=sorted(group[scenario,p,rep],key=lambda x:int(x["rank"]))
        need(len(cand)==6 and
             [int(x["rank"]) for x in cand]==list(range(1,7)),
             "wrong candidate rank inventory")
        earlier=float("inf")
        seen=set()
        true_ix=[]
        for x in cand:
            a,b,c,d=[int(x[k]) for k in
                     ("pairA_a","pairA_b","pairB_a","pairB_b")]
            need(1<=a<b<=p and 1<=c<d<=p and
                 not ({a,b}&{c,d}),"pair orientation/disjoint violation")
            pair=((a,b),(c,d))
            need(pair not in seen,"duplicate shortlisted interaction")
            seen.add(pair)
            observed=score(X,y,a,b,c,d)
            printed=float(x["score"])
            need(abs(observed-printed)<2e-7,
                 "RNG re-generation / saved D interaction score mismatch")
            need(printed<=earlier+1e-7 and 0<=printed<=1,
                 "top-six score order violation")
            earlier=printed
            t=int({(a,b),(c,d)}=={(1,2),(3,4)})
            need(int(x["true_pair_product"])==t,
                 "planted pair flags inconsistent")
            if t:true_ix.append(int(x["rank"]))
            need(int(x["seed"])==seed and
                 x["method"]=="BLIND_LSH_PAIR_PRODUCT",
                 "shortlisted patient cohort/method drift")
        need(len(true_ix)<=1 and
             int(row["true_joint_in_top6"])==int(bool(true_ix)) and
             int(row["true_joint_top6_rank"])==
                 (true_ix[0] if true_ix else 0) and
             int(row["true_joint_hash_collision"])>=int(bool(true_ix)) and
             abs(float(row["best_score"])-float(cand[0]["score"]))<1e-7,
             "true pair recovery / hash collision flag discrepancy")
    print("=== H3 BLIND XOR PAIR-PRODUCT RECOVERY ===")
    for scenario in SCENARIOS:
        for p in PS:
            rr=[meta[scenario,p,k] for k in range(1,REPS+1)]
            collisions=sum(int(x["true_joint_hash_collision"]) for x in rr)
            hits=sum(int(x["true_joint_in_top6"]) for x in rr)
            print(f"{scenario:4s} p={p:3d} "
                  f"true_hash_collision={collisions}/{REPS} "
                  f"true_top6={hits}/{REPS}")
    print("STAGEH3_240_TOP6_SCORES_RECOMPUTED_FROM_FRESH_D=PASS")
    print("STAGEH3_40_STUDY_BLIND_SOURCE_INVENTORY_AUDITED=PASS")
    print("STAGEH3_SCIENTIFIC_STAGE=D_ONLY_SYNTHETIC_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
