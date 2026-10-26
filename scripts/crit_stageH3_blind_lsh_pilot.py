#!/usr/bin/env python3
"""CRIT H3: blind XOR pair-of-rank-pairs discovery by binary Hamming LSH.

All training samples are D, independent synthetic; true pair IDs are used
ONLY for retrospective scoring, NEVER proposal/selection. This does not fit
RLR, use V, or evaluate an E set. Hashes depend on D outcomes and can overfit.
"""
from __future__ import annotations
import csv
import math
import resource
import time
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent.parent
OUT=ROOT/"results/crit_stageH3_blind_lsh_fix1_v1"
PS=(32,64,128,256,512)
REPS=4
N=120
HASH_BITS=18
HASH_TABLES=24
SCORE_BATCH=4096
TRUTH=((0,1),(2,3))
SCENARIOS=("null","xor")

def need(ok,why):
    if not ok:
        raise RuntimeError("STAGEH3_STOP: "+why)

def sample_D(p:int,seed:int,scenario:str):
    rng=np.random.default_rng(seed)
    X=rng.standard_normal((N,p)).astype(np.float32)
    z=np.where(X[:,0]>X[:,1],1.0,-1.0)
    r=np.where(X[:,2]>X[:,3],1.0,-1.0)
    logit=2.4*z*r if scenario=="xor" else np.zeros(N,dtype=float)
    prob=1/(1+np.exp(-logit))
    y=(rng.random(N)<prob).astype(np.uint8)
    need(np.unique(y).size==2,"binary label support missing")
    return X,y

def make_signatures(X:np.ndarray):
    p=X.shape[1]
    a,b=np.triu_indices(p,k=1)
    # Each signature is a fixed binary rank comparison across 120 patients.
    S=np.ascontiguousarray((X[:,a]>X[:,b]).T)
    need(S.shape==(p*(p-1)//2,N),"signature dimensions wrong")
    return S,a.astype(np.int32),b.astype(np.int32)

def join_bucket(code, ycode:int, opposite:bool):
    """Vectorized hash-bucket join for R_i XOR R_j ~= (NOT) y bits."""
    M=len(code)
    mask=(1<<HASH_BITS)-1
    q=(code ^ ycode ^ (mask if opposite else 0)).astype(np.int32)
    sorted_ix=np.argsort(code,kind="stable").astype(np.int32)
    ordered=code[sorted_ix]
    left=np.searchsorted(ordered,q,side="left")
    right=np.searchsorted(ordered,q,side="right")
    counts=right-left
    total=int(counts.sum())
    if not total:
        return np.empty(0,dtype=np.int32),np.empty(0,dtype=np.int32)
    src=np.repeat(np.arange(M,dtype=np.int32),counts)
    lower=np.repeat(np.cumsum(counts,dtype=np.int64)-counts,counts)
    offsets=np.arange(total,dtype=np.int64)-lower
    dest=sorted_ix[np.repeat(left,counts)+offsets]
    idx=src<dest
    return src[idx],dest[idx]

def blind_rank_pair_lsh(X,y,seed:int):
    S,pa,pb=make_signatures(X)
    M=len(pa)
    true_a=int(np.where((pa==0)&(pb==1))[0][0])
    true_b=int(np.where((pa==2)&(pb==3))[0][0])
    truth=(min(true_a,true_b),max(true_a,true_b))
    rng=np.random.default_rng(seed+93217)
    ybits=np.asarray(y,dtype=np.uint8)
    powers=(1<<np.arange(HASH_BITS,dtype=np.int64))
    seen_truth=False
    proposals=0
    scored=0
    best={}
    t0=time.perf_counter()
    for tab in range(HASH_TABLES):
        patient_ix=rng.choice(N,HASH_BITS,replace=False)
        code=(S[:,patient_ix].astype(np.int32) @ powers).astype(np.int32)
        ycode=int(ybits[patient_ix].astype(np.int32) @ powers)
        for flip in (False,True):
            src,dst=join_bucket(code,ycode,flip)
            if not len(src):
                continue
            allowed=(pa[src]!=pa[dst])&(pa[src]!=pb[dst])&(
                pb[src]!=pa[dst])&(pb[src]!=pb[dst])
            src=src[allowed]
            dst=dst[allowed]
            proposals+=len(src)
            if not len(src):
                continue
            seen_truth |= bool(np.any((src==truth[0])&(dst==truth[1])))
            scores=np.empty(len(src),dtype=np.float32)
            for off in range(0,len(src),SCORE_BATCH):
                end=min(off+SCORE_BATCH,len(src))
                agreement=np.count_nonzero(
                    ((S[src[off:end]] ^ S[dst[off:end]]) !=
                     ybits[None,:]),axis=1)
                scores[off:end]=abs((2.0*agreement/N)-1)
            scored+=len(src)
            take=min(6,len(src))
            candidates=np.argpartition(-scores,take-1)[:take]
            for j in candidates:
                key=(int(src[j]),int(dst[j]))
                val=float(scores[j])
                if key not in best or best[key]<val:
                    best[key]=val
    ranked=sorted(best.items(),
                  key=lambda kv:(-kv[1],kv[0][0],kv[0][1]))
    chosen=ranked[:6]
    exact_hit=any(pair==truth for pair,_ in chosen)
    rank=next((i+1 for i,(pair,_) in enumerate(chosen)
               if pair==truth),None)
    rows=[]
    for i,(pair,score) in enumerate(chosen):
        j,k=pair
        rows.append(dict(rank=i+1,pairA_a=int(pa[j])+1,
                         pairA_b=int(pb[j])+1,
                         pairB_a=int(pa[k])+1,
                         pairB_b=int(pb[k])+1,
                         score=score,
                         true_pair_product=int(pair==truth)))
    return dict(
        pair_feature_count=M,tables=HASH_TABLES,bits=HASH_BITS,
        candidate_pairs_proposed=proposals,
        candidate_pairs_scored=scored,
        true_joint_hash_collision=int(seen_truth),
        true_joint_in_top6=int(exact_hit),
        true_joint_top6_rank=rank or 0,
        best_score=float(chosen[0][1]) if chosen else float("nan"),
        seconds=time.perf_counter()-t0,
        peak_rss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        top=rows
    )

def main():
    need(OUT.is_dir() and (OUT/"SOURCE_SHA.txt").is_file() and
         not (OUT/"study_summary.csv").exists(),
         "output missing or already written")
    summary=[]
    top=[]
    for sc_index,sc in enumerate(SCENARIOS,1):
      for p in PS:
        for rep in range(1,REPS+1):
          seed=20270211+sc_index*100003+p*317+rep*43
          D,y=sample_D(p,seed,sc)
          v=blind_rank_pair_lsh(D,y,seed)
          need(v["candidate_pairs_scored"]>0 and len(v["top"])==6,
               "LSH produced fewer than six candidates")
          key={"scenario":sc,"p":p,"replicate":rep,"seed":seed,
               "n_D":N,"method":"BLIND_LSH_PAIR_PRODUCT",
               "scientific_stage":"D_ONLY_SYNTHETIC_DEVELOPMENT"}
          summary.append(dict(
            **key,
            pair_feature_count=v["pair_feature_count"],
            hash_tables=v["tables"],hash_bits=v["bits"],
            candidate_pairs_proposed=v["candidate_pairs_proposed"],
            candidate_pairs_scored=v["candidate_pairs_scored"],
            true_joint_hash_collision=v["true_joint_hash_collision"],
            true_joint_in_top6=v["true_joint_in_top6"],
            true_joint_top6_rank=v["true_joint_top6_rank"],
            best_score=v["best_score"],seconds=v["seconds"],
            peak_rss_kb=v["peak_rss_kb"]))
          for entry in v["top"]:
            top.append(dict(**{k:key[k] for k in
                         ("scenario","p","replicate","seed","method")},
                         **entry))
          print(
            f"H3_CASE {sc} p={p} rep={rep} M={v['pair_feature_count']} "
            f"proposals={v['candidate_pairs_scored']} "
            f"truth_collision={v['true_joint_hash_collision']} "
            f"truth_top6={v['true_joint_in_top6']} "
            f"best={v['best_score']:.4f} secs={v['seconds']:.2f}",
            flush=True)
    need(len(summary)==2*len(PS)*REPS and len(top)==6*len(summary),
         "H3 40 x 6 inventory missing")
    for filename,data in (("study_summary.csv",summary),
                          ("top6_candidates.csv",top)):
        with (OUT/filename).open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=data[0].keys())
            w.writeheader();w.writerows(data)
    print("STAGEH3_40_BLIND_D_ONLY_STUDIES=PASS")
    print("STAGEH3_24_TABLE_LSH_XOR_PAIR_DISCOVERY=PASS")
    print("STAGEH3_TOP6_PROVENANCE=PASS")
    print("STAGEH3_SCIENTIFIC_STAGE=D_ONLY_SYNTHETIC_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
