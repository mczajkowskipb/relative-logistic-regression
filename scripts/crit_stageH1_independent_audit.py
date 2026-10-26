#!/usr/bin/env python3
"""Independent H1 screening archive structural / provenance auditor.

Rechecks all 256 method-study inventories, pair-budget/canonical orientation,
D-only pair membership and top-six independent unordered interaction indices.
No claim to reconstruct R RNG or evaluate original-RLR or held-out models.
"""
from __future__ import annotations
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_stageH1_bounded_screening_v1"
PS=(32,64,128,256,512)
SCENARIOS=("null","xor")
STRATEGIES=("FIXED_RANDOM","SEED_RANDOM","AWARE_GRAPH")
REFS={(1,2),(3,4)}

def stop(ok,msg):
    if not ok:
        raise SystemExit("STAGEH1_AUDIT_STOP: "+msg)

def load(path):
    stop(path.is_file(),"missing "+str(path))
    with path.open("r",newline="",encoding="utf-8") as f:
        a=list(csv.DictReader(f))
    stop(bool(a),"empty "+str(path))
    return a

def main():
    summ=load(OUT/"screening_summary.csv")
    pools=load(OUT/"candidate_pools.csv")
    top6=load(OUT/"top6_pair_interactions.csv")
    expected={(s,p,r,m)
              for s in SCENARIOS for p in PS for r in range(1,9)
              for m in (STRATEGIES+("EXHAUSTIVE32",) if p==32
                        else STRATEGIES)}
    summary={(r["scenario"],int(r["p"]),int(r["replicate"]),r["method"]):r
             for r in summ}
    stop(len(summary)==len(expected)==256 and set(summary)==expected,
         "summary duplicate/missing study")
    pg=defaultdict(list)
    for r in pools:
        k=(r["scenario"],int(r["p"]),int(r["replicate"]),r["method"])
        stop(k in expected,"unexpected candidate pool id")
        pg[k].append(r)
    tg=defaultdict(list)
    for r in top6:
        k=(r["scenario"],int(r["p"]),int(r["replicate"]),r["method"])
        stop(k in expected,"unexpected interaction shortlist ID")
        tg[k].append(r)
    stop(set(pg)==expected and set(tg)==expected,"missing pool/top6 IDs")
    stop(len(top6)==256*6 and len(pools)==256*320+16*(496-320),
         "incorrect total candidate/interaction counts")
    for key in sorted(expected):
        s,p,rep,method=key
        row=summary[key]
        c=sorted(pg[key],key=lambda x:int(x["candidate_index"]))
        K=496 if method=="EXHAUSTIVE32" else 320
        stop(len(c)==K and int(row["pool_budget"])==K and
             [int(x["candidate_index"]) for x in c]==list(range(1,K+1)),
             "candidate pool budget / indices invalid")
        pairs=[]
        for x in c:
            a=int(x["gene_a"])
            b=int(x["gene_b"])
            stop(1<=a<b<=p,"candidate unordered canonical pair invalid")
            pairs.append((a,b))
        stop(len(set(pairs))==K,"candidate pool duplicated")
        if method=="EXHAUSTIVE32":
            stop(p==32 and
                 set(pairs)=={(i,j) for i in range(1,33)
                             for j in range(i+1,33)},
                 "exhaustive32 pool is not exact")
        first=(1,2) in pairs
        second=(3,4) in pairs
        stop(int(row["reference_gate_in_pool"])==int(first) and
             int(row["reference_pair_in_pool"])==int(second) and
             int(row["both_reference_pairs_in_pool"])==int(first and second),
             "reference pair coverage flags not matched")
        matches=[]
        rs=sorted(tg[key],key=lambda x:int(x["joint_rank"]))
        stop(len(rs)==6 and
             [int(x["joint_rank"]) for x in rs]==list(range(1,7)),
             "not six joint-ranked results")
        past=float("inf")
        used=set()
        for x in rs:
            a=int(x["pairA_index"])
            b=int(x["pairB_index"])
            score=float(x["score"])
            stop(math.isfinite(score) and score>=0 and score<=past+1e-7,
                 "not descending finite interaction score")
            past=score
            stop(1<=a<b<=K and frozenset((a,b)) not in used,
                 "joint candidate index order invalid")
            used.add(frozenset((a,b)))
            pairA=pairs[a-1]
            pairB=pairs[b-1]
            stop(not (set(pairA)&set(pairB)),
                 "joint candidate gene sets overlap")
            matches.append({pairA,pairB}==REFS)
        stop(int(row["reference_joint_in_top6"])==int(any(matches)),
             "orientation-invariant truth joint top6 doesn't replay")
        stop(int(row["seed"])==
               20261026+(SCENARIOS.index(s)+1)*100000+p*300+rep*19 and
             int(row["n_D"])==120,
             "predeclared seed/support mismatch")
        stop(row["scientific_stage"]=="D_ONLY_SYNTHETIC_SCREENING_DEV",
             "wrong stage label")
        t=float(row["screening_seconds"])
        stop(math.isfinite(t) and t>=0,"runtime missing")
        if method=="AWARE_GRAPH":
            ta=float(row["aware_graph_seconds"])
            stop(math.isfinite(ta) and ta>=0,
                 "aware graph resource time missing")
        else:
            stop(row["aware_graph_seconds"] in ("","NA"),
                 "unexpected graph run on non-aware strategy")
    print("=== STAGE H1 D-ONLY SCREENING RESULTS ===")
    for s in SCENARIOS:
      for p in PS:
       for method in STRATEGIES+("EXHAUSTIVE32",) if p==32 else STRATEGIES:
        vals=[summary[s,p,r,method] for r in range(1,9)]
        covered=sum(int(r["both_reference_pairs_in_pool"]) for r in vals)
        joint=sum(int(r["reference_joint_in_top6"]) for r in vals)
        cost=mean(float(r["screening_seconds"]) for r in vals)
        graph=(mean(float(r["aware_graph_seconds"]) for r in vals)
               if method=="AWARE_GRAPH" else 0.0)
        print(f"{s:4s} p={p:3d} {method:15s} "
              f"both_pairs={covered}/8 pair_interaction_top6={joint}/8 "
              f"score_s={cost:.4f} graph_s={graph:.4f}")
    print("STAGEH1_256_INDEPENDENT_PROVENANCE_CHECKS=PASS")
    print("STAGEH1_D_ONLY_PAIR_BUDGETS_AND_SYMMETRIC_JOINT_AUDITED=PASS")
    print("STAGEH1_SCIENTIFIC_STAGE=SCREENING_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
