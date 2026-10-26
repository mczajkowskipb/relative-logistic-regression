#!/usr/bin/env python3
"""Independent artifact auditor for Stage H2 ORACLE anchor diagnostic."""
from __future__ import annotations
import csv
import math
from pathlib import Path
from collections import defaultdict
from statistics import mean

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_stageH2_anchored_completion_v1"
SCENARIOS=("null","xor")
PS=(32,64,128,256,512,1024)
ARMS={
    "ORACLE_Z_TO_R":((1,2),(3,4)),
    "ORACLE_R_TO_Z":((3,4),(1,2)),
    "SHAM_56_TO_R":((5,6),(3,4)),
}
COLS=("scenario","p","replicate","method")

def check(ok,msg):
    if not ok:raise SystemExit("STAGEH2_AUDIT_STOP: "+msg)

def load(p):
    check(p.is_file(),"archive missing "+str(p))
    with p.open(encoding="utf-8",newline="") as f: return list(csv.DictReader(f))

def f(x):
    z=float(x)
    check(math.isfinite(z),"nonfinite numeric")
    return z

def main():
    s=load(OUT/"study_summary.csv")
    t=load(OUT/"top6_candidates.csv")
    keys={(sc,p,r,m) for sc in SCENARIOS for p in PS
          for r in range(1,9) for m in ARMS}
    summaries={(x["scenario"],int(x["p"]),int(x["replicate"]),x["method"]):x
               for x in s}
    check(len(s)==288 and len(summaries)==288 and set(summaries)==keys,
          "unexpected 96-study / 288-scan inventory")
    blocks=defaultdict(list)
    for x in t:
        key=(x["scenario"],int(x["p"]),int(x["replicate"]),x["method"])
        check(key in keys,"unexpected candidate shortlist group")
        blocks[key].append(x)
    check(len(t)==1728 and set(blocks)==keys,
          "missing or duplicate top-six blocks")
    for sc,p,rep,method in sorted(keys):
        key=(sc,p,rep,method)
        x=summaries[key]
        pool=(p-2)*(p-3)//2
        anchor,target=ARMS[method]
        check(x["scientific_stage"]==
              "ORACLE_ANCHOR_D_ONLY_SYNTHETIC_DEVELOPMENT",
              "invalid scientific-stage label")
        check(x["scenario"]==sc and x["method"]==method and
              int(x["n_D"])==120 and
              int(x["seed"])==20270117+(SCENARIOS.index(sc)+1)*100000+
                                    p*311+rep*37,
              "misidentified D study or seed")
        check((int(x["anchor_gene_a"]),int(x["anchor_gene_b"]))==anchor and
              (int(x["target_gene_a"]),int(x["target_gene_b"]))==target and
              int(x["candidate_pair_count"])==pool,
              "oracle/sham anchor or pool cardinality mismatch")
        rank=int(x["true_target_rank"])
        check(1<=rank<=pool,"planted target candidate has invalid rank")
        check(int(x["target_in_top1"])==int(rank==1) and
              int(x["target_in_top6"])==int(rank<=6),
              "target rank / top-k flags inconsistent")
        check(f(x["screening_seconds"])>=0 and
              0<=f(x["true_target_score"])<=f(x["max_score"])+1e-7,
              "invalid screening run cost or target score")
        M=p*(p-1)//2
        expectation=320/M*319/(M-1)
        check(abs(f(x["random_two_pair_coverage_K320"])-expectation)<1e-12,
              "analytic random 2-pair inclusion chance wrong")
        group=sorted(blocks[key],key=lambda r:int(r["joint_rank"]))
        check(len(group)==6 and
              [int(v["joint_rank"]) for v in group]==list(range(1,7)),
              "top-six ranks incomplete")
        prior=float("inf")
        found=0
        for v in group:
            a=int(v["found_pair_a"]);b=int(v["found_pair_b"])
            score=f(v["score"])
            check(1<=a<b<=p and
                  a not in anchor and b not in anchor,
                  "candidate target pair overlaps oracle anchor")
            check(score<=prior+1e-7 and score>=0,
                  "not a descending nonnegative score")
            prior=score
            is_true=int((a,b)==target)
            check(int(v["true_target"])==is_true,
                  "saved truth flag disagrees with raw pair ID")
            found+=is_true
        check(found==int(rank<=6),
              "reported top-six target rank not in saved top-six")
        check(abs(f(group[0]["score"])-f(x["max_score"]))<1e-9,
              "maximum score differs from top-1 candidate")
        if rank<=6:
            z=[v for v in group if int(v["true_target"])==1]
            check(len(z)==1 and int(z[0]["joint_rank"])==rank and
                  abs(f(z[0]["score"])-f(x["true_target_score"]))<1e-9,
                  "true target rank/score mismatch")
    print("=== H2 TRUE PAIR COMPLETION UNDER ORACLE ANCHOR ===")
    for sc in SCENARIOS:
      for p in PS:
        for method in ARMS:
            rr=[summaries[sc,p,rep,method] for rep in range(1,9)]
            count1=sum(int(x["target_in_top1"]) for x in rr)
            count6=sum(int(x["target_in_top6"]) for x in rr)
            sec=mean(f(x["screening_seconds"]) for x in rr)
            print(f"{sc:4s} p={p:4d} {method:15s} "
                  f"top1={count1}/8 top6={count6}/8 scan_s={sec:.3f}")
    print("STAGEH2_288_PAIR_AND_METRIC_PROVENANCE=PASS")
    print("STAGEH2_ORACLE_VS_SHAM_TARGET_IDENTITIES_AUDITED=PASS")
    print("STAGEH2_NOT_A_FULL_RLR_CLASSIFIER")
    print("STAGEH2_SCIENTIFIC_STAGE=D_ONLY_SYNTHETIC_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
