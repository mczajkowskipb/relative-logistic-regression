#!/usr/bin/env python3
"""Independent Stage H4A archive audit: regenerate D, verify every top-six score.

No imports of the original H3 estimator; independent NumPy generator and
score reconstruction. Verify predeclared seeds, design, truth orientation,
shortlists, candidate counts, cumulative process RSS and source identity.
"""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_stageH4A_blind_robustness_v1"
NS=(60,120,240)
PS=(32,128,512)
AMPS=(0.0,1.0,1.6,2.4)
REPS=8
SEED0=20280617

def need(ok,what):
    if not ok:raise SystemExit("STAGEH4A_AUDIT_STOP: "+what)

def read_csv(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def regenerate(n,p,amp,seed):
    rng=np.random.default_rng(seed)
    X=rng.standard_normal((n,p)).astype(np.float32)
    z=np.where(X[:,0]>X[:,1],1.0,-1.0)
    r=np.where(X[:,2]>X[:,3],1.0,-1.0)
    eta=amp*z*r
    y=(rng.random(n)<1/(1+np.exp(-eta))).astype(np.uint8)
    return X,y

def independent_score(X,y,pairs):
    (a,b),(c,d)=pairs
    A=X[:,a-1]>X[:,b-1]
    B=X[:,c-1]>X[:,d-1]
    # Independent expression of H3: abs(2 * misagreement/n - 1).
    disagree=np.count_nonzero((A^B)!=y)
    return abs(2.0*disagree/len(y)-1.0)

def wilson(k,n,z=1.95996398454):
    x=k/n
    den=1+z*z/n
    center=(x+z*z/(2*n))/den
    half=z*math.sqrt(x*(1-x)/n+z*z/(4*n*n))/den
    return (center-half,center+half)

def main():
    m=read_csv(OUT/"study_summary.csv")
    t=read_csv(OUT/"top6_candidates.csv")
    expected={(n,p,ai,rep) for n in NS for p in PS
              for ai in range(len(AMPS)) for rep in range(1,REPS+1)}
    study={(int(x["n_D"]),int(x["p"]),
            int(x["amplitude_index"]),int(x["replicate"])):x for x in m}
    need(len(m)==len(expected)==288 and set(study)==expected,
         "study inventory missing or duplicated")
    groups=defaultdict(list)
    for x in t:
        key=(int(x["n_D"]),int(x["p"]),
             int(x["amplitude_index"]),int(x["replicate"]))
        need(key in expected,"unknown finalist cohort")
        groups[key].append(x)
    checked=0
    for key in sorted(expected):
        n,p,ai,rep=key
        row=study[key]
        amp=AMPS[ai]
        seed=SEED0+n*101+p*313+ai*100003+rep*47
        X,y=regenerate(n,p,amp,seed)
        need(int(row["seed"])==seed and
             abs(float(row["amplitude"])-amp)<1e-12 and
             int(row["label_positive"])==int(y.sum()) and
             0<int(y.sum())<n,
             "sample generator or seed drift")
        need(row["method"]=="FROZEN_H3_LSH_18BITS_24TABLES" and
             row["endpoint"]=="TRUE_PAIR_PRODUCT_IN_TOP6_D_ONLY" and
             row["scientific_stage"]==
               "INDEPENDENT_SIMULATION_DEVELOPMENT_ONLY" and
             int(row["hash_bits"])==18 and
             int(row["hash_tables"])==24 and
             int(row["pair_feature_count"])==p*(p-1)//2,
             "parameter, method identity or stage drift")
        proposed=int(row["candidate_pairs_proposed"])
        scored=int(row["candidate_pairs_scored"])
        need(proposed==scored and proposed>=0 and
             int(row["process_peak_rss_kb"])>0 and
             math.isfinite(float(row["wall_seconds"])) and
             float(row["wall_seconds"])>=0,
             "invalid candidate/running inventory")
        block=sorted(groups.get(key,[]),key=lambda x:int(x["rank"]))
        k=int(row["shortlist_size"])
        need(len(block)==k and 0<=k<=6 and
             [int(x["rank"]) for x in block]==list(range(1,k+1)),
             "shortlist cardinality/ranks")
        need((k>0)==(proposed>0),"non-empty shortlist/proposal mismatch")
        previously=float("inf")
        seen=set()
        truthrank=0
        for candidate in block:
            a,b,c,d=(int(candidate[name]) for name in
                ("pairA_a","pairA_b","pairB_a","pairB_b"))
            pairs=((a,b),(c,d))
            need(1<=a<b<=p and 1<=c<d<=p and
                 not ({a,b}&{c,d}) and pairs not in seen,
                 "pair IDs noncanonical, overlapping or duplicated")
            seen.add(pairs)
            raw=independent_score(X,y,pairs)
            val=float(candidate["score"])
            need(math.isfinite(val) and 0<=val<=1.000001 and
                 abs(raw-val)<3e-7 and
                 val<=previously+3e-7,
                 "independent D score or rank ordering incorrect")
            previously=val
            wastrue=int({(a,b),(c,d)}=={(1,2),(3,4)})
            need(int(candidate["true_pair_product"])==wastrue and
                 int(candidate["seed"])==seed,
                 "true-pair or seed label mismatch")
            if wastrue:
                need(not truthrank,"two true pair products in shortlist")
                truthrank=int(candidate["rank"])
            checked+=1
        hit=int(row["true_joint_in_top6"])
        collision=int(row["true_joint_hash_collision"])
        need((collision,hit) in ((0,0),(1,0),(1,1)) and
             hit==int(truthrank>0) and
             int(row["true_joint_top6_rank"])==truthrank,
             "collision/hit/rank inconsistency")
        if k:
            need(abs(float(row["best_score"]) -
                     float(block[0]["score"]))<3e-7,
                 "best score disagrees with rank1")
        else:
            need(math.isnan(float(row["best_score"])),
                 "empty shortlist should have no score")
    print("=== FROZEN BLIND LSH: INDEPENDENT SYNTHETIC ROBUSTNESS ===")
    for ai,amp in enumerate(AMPS):
      for n in NS:
        for p in PS:
            rows=[study[n,p,ai,rep] for rep in range(1,REPS+1)]
            k=sum(int(x["true_joint_in_top6"]) for x in rows)
            coll=sum(int(x["true_joint_hash_collision"]) for x in rows)
            lo,hi=wilson(k,REPS)
            runtime=mean(float(x["wall_seconds"]) for x in rows)
            print(f"H4A_FINAL signal={amp:.1f} n={n} p={p} "
                  f"top6={k}/{REPS} collision={coll}/{REPS} "
                  f"wilson95=[{lo:.3f},{hi:.3f}] mean_sec={runtime:.3f}")
    print(f"STAGEH4A_{checked}_TOPK_SCORES_INDEPENDENTLY_REPLAYED=PASS")
    print("STAGEH4A_288_SOURCE_SEED_AND_FEATURE_AUDITS=PASS")
    print("STAGEH4A_36_CELLS_AUDITED=PASS")
    print("STAGEH4A_NULL_REFERENCE_NOT_TYPE1_ERROR")

if __name__=="__main__":
    main()
