#!/usr/bin/env python3
"""Differential unit tests for H3 bucket joins and planted XOR retrieval."""
import itertools
import numpy as np
from crit_stageH3_blind_lsh_pilot import (
    HASH_BITS,N,blind_rank_pair_lsh,join_bucket
)

def test_bucket():
    codes=np.array([0,1,1,2,5,7,11,255],dtype=np.int32)
    ycode=3
    for flip in (False,True):
        a,b=join_bucket(codes,ycode,flip)
        got=set(zip(a.tolist(),b.tolist()))
        mask=(1<<HASH_BITS)-1
        goal={ (i,j) for i,j in itertools.combinations(range(len(codes)),2)
             if (int(codes[i])^int(codes[j])) ==
                (ycode ^ (mask if flip else 0)) }
        assert got==goal,(flip,got,goal)
    print("STAGEH3_EXHAUSTIVE_HASH_BUCKET_DIFFERENTIAL=PASS")

def test_noiseless_xor():
    rng=np.random.default_rng(202612345)
    X=rng.standard_normal((N,16)).astype(np.float32)
    z=(X[:,0]>X[:,1])
    r=(X[:,2]>X[:,3])
    y=(~(z^r)).astype(np.uint8)
    got=blind_rank_pair_lsh(X,y,seed=459)
    assert got["true_joint_hash_collision"]==1,got
    assert got["true_joint_in_top6"]==1,got
    assert got["true_joint_top6_rank"]==1,got
    assert len(got["top"])==6,got
    assert abs(got["top"][0]["score"]-1.0)<1e-6,got
    print("STAGEH3_NOISELESS_BLIND_XOR_EXACT_RECOVERY=PASS")

if __name__=="__main__":
    test_bucket()
    test_noiseless_xor()
    print("STAGEH3_UNIT_TESTS=PASS")
