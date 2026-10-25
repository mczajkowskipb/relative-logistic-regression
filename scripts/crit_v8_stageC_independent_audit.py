#!/usr/bin/env python3
"""Independent arithmetic, inventory and frozen-manifest check of RDS-derived Stage C."""
from __future__ import annotations
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_v8_stageC_selection_v1"
DATASETS=("chen-2002","gordon-2002")
FOLDS=("rep00_fold00","rep00_fold01","rep00_fold02")
ARMS=("GLOBAL_CAL","GATE_ONLY","ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")

def require(ok:bool,msg:str)->None:
    if not ok:raise SystemExit("STAGEC_INDEPENDENT_AUDIT_STOP: "+msg)

def data(path:Path)->list[dict[str,str]]:
    require(path.is_file(),"file missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        rows=list(csv.DictReader(f))
    require(bool(rows),"empty csv "+str(path))
    return rows

def number(x:str)->float:
    v=float(x)
    require(math.isfinite(v),"nonfinite numerical statistic")
    return v

def close(a:float,b:float,tol:float=1e-8)->bool:
    return abs(a-b)<=tol

def main():
    require(len(sys.argv)==3,"usage auditor.py CRIT7_REF CRIT8_REF")
    v7=Path(sys.argv[1]).resolve(strict=True)
    v8=Path(sys.argv[2]).resolve(strict=True)
    old7=v7/"results/crit_v7_original_rlr_2x3_v1"
    old8=v8/"results/crit_v8_one_shot_2x3_v1"
    folds=data(OUT/"selection_folds.csv")
    candidates=data(OUT/"v7_candidate_rankings.csv")
    arms=data(OUT/"v8_confirmation_arms.csv")
    m7=data(old7/"run_manifest.csv")
    m8=data(old8/"run_manifest.csv")
    require(len(folds)==6 and len(arms)==30 and
            len(m7)==24 and len(m8)==36,"unexpected file dimensions")
    keys={(ds,fold) for ds in DATASETS for fold in FOLDS}
    fold_map={(r["dataset"],r["split_id"]):r for r in folds}
    require(set(fold_map)==keys and len(fold_map)==6,
            "missing or duplicate fold keys")
    m7_map={(r["dataset"],r["split_id"],r["method"]):r for r in m7}
    m8_map={(r["dataset"],r["split_id"],r["method"]):r for r in m8}
    grouped=defaultdict(list)
    for r in candidates:
        k=(r["dataset"],r["split_id"])
        require(k in keys,"unexpected v7 fold")
        grouped[k].append(r)
    agroups=defaultdict(dict)
    for r in arms:
        k=(r["dataset"],r["split_id"])
        require(k in keys and r["method"] in ARMS,"unknown V arm")
        require(r["method"] not in agroups[k],"duplicate V arm")
        agroups[k][r["method"]]=r
    require(set(grouped)==keys and set(agroups)==keys,"missing fold data")
    for k in sorted(keys):
        f=fold_map[k]
        ds,sid=k
        cs=sorted(grouped[k],key=lambda r:int(r["D_rank"]))
        require(1<=len(cs)<=6 and
                [int(r["D_rank"]) for r in cs]==list(range(1,len(cs)+1)),
                "D rank sequence corrupted")
        lower=[number(r["confirmation_lower"]) for r in cs]
        meanv=[number(r["confirmation_mean"]) for r in cs]
        se=[number(r["confirmation_se"]) for r in cs]
        require(all(close(z,mu-s-0.002)
                    for z,mu,s in zip(lower,meanv,se)),
                "v7 confirmation arithmetic does not match 1SE+0.002")
        best=max(range(len(cs)),key=lambda i:lower[i])+1
        require(best==int(f["v7_best_D_rank"]) and
                close(max(lower),number(f["v7_best_confirm_lower"])) and
                close(lower[0],number(f["first_D_confirm_lower"])) and
                int(f["first_D_equals_v7_best"])==int(best==1),
                "v7 candidate ranking/selection inconsistent")
        require(sum(int(r["selected_by_v7"]) for r in cs)==1 and
                all(int(r["selected_by_v7"])==int(i+1==best)
                    for i,r in enumerate(cs)),
                "v7 best candidate tag incorrect")
        require(all(int(r["first_D_nominee"])==int(i==0)
                    for i,r in enumerate(cs)),
                "D first candidate tag incorrect")
        expected7=max(lower)>0.002
        expected8=lower[0]>0.002
        require(int(f["v7_accepted"])==int(expected7) and
                int(f["v8_accepted"])==int(expected8),
                "one-shot / candidate-maximum decisions incorrect")
        require(close(number(f["v8_confirm_lower"]),lower[0]) and
                number(f["v7_v8_first_candidate_parity_err"])<1e-8,
                "v8 nomination differs from v7 first D candidate")
        require(int(f["candidate_count"])==len(cs) and
                int(f["n_D"])+int(f["n_V"]) in (143,144,145),
                "candidate count or outer train support incorrect")
        old7row=m7_map[(ds,sid,"CRIT_TREE")]
        old8row=m8_map[(ds,sid,"CONDITIONAL_TREE")]
        require(int(old7row["tree_accepted"])==int(expected7) and
                int(old8row["accepted"])==int(expected8) and
                close(number(old7row["delta_vs_original"]),
                      number(f["v7_outer_deltaLL"])) and
                close(number(old8row["delta_vs_original"]),
                      number(f["v8_outer_deltaLL"])),
                "archived test scores or acceptance mismatched")

        method_map=agroups[k]
        require(set(method_map)==set(ARMS),"incomplete V arm inventory")
        for method, arm in method_map.items():
            available=int(arm["available"])
            accepted=int(arm["accepted"])
            require(available==1,"unexpected unavailable correction")
            mu=number(arm["validation_mean"])
            se=number(arm["validation_se"])
            low=number(arm["confirmation_lower"])
            require(close(low,mu-se-0.002) and
                    accepted==int(low>0.002),
                    "one-shot V arm decision arithmetic incorrect")
            historic=m8_map[(ds,sid,method)]
            require(accepted==int(historic["accepted"]) and
                    close(low,number(historic["confirmation_lower"])) and
                    close(number(arm["outer_deltaLL"]),
                          number(historic["delta_vs_original"])),
                    "v8 archived one-shot test/outer audit discrepancy")
        require(close(number(method_map["CONDITIONAL_TREE"]["confirmation_lower"]),
                      lower[0]),"candidate #1 conditional score mismatch")
        print(f"{ds}:{sid} candidates={len(cs)} V_best_D_rank={best} "
              f"V_best_lower={max(lower):+.7f} D_first_lower={lower[0]:+.7f} "
              f"v7_accept={int(expected7)} v8_accept={int(expected8)}")
    print("STAGEC_V7_V8_SELECTION_ARITHMETIC=PASS")
    print("STAGEC_ALL_30_VALIDATION_DECISIONS_AUDITED=PASS")
    print("STAGEC_FROZEN_MANIFEST_PARITY=PASS")
    print("STAGEC_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")

if __name__=="__main__":
    main()
