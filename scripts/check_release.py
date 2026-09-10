#!/usr/bin/env python3
from pathlib import Path
import argparse, hashlib, json, re, sys

ROOT=Path(__file__).resolve().parents[1]

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()

def parse_hashes(p):
    out={}
    for line in p.read_text().splitlines():
        if not line.strip(): continue
        h,rel=line.split(None,1)
        out[rel.strip()]=h
    return out

ap=argparse.ArgumentParser()
ap.add_argument('--quick', action='store_true')
a=ap.parse_args()

cfgdir=ROOT/'configs/paper_frozen'
req=['compcancer_9x10.json','compcancer_2x2.json','compcancer_collect.json','ptb_fixed.json','ptb_loco.json','legacy8.json']
for n in req:
    p=cfgdir/n
    assert p.is_file(), f'missing {p}'
    c=json.loads(p.read_text())
    m=c['models']
    assert m['seed']==42
    assert m['top_pool_values']==[1000]
    assert m['alpha']==0.5
    assert m['lambda_choice']=='lambda.min'
    assert m['rlr_grid']['m']==[10] and m['rlr_grid']['d']==[2]

assert json.loads((cfgdir/'compcancer_9x10.json').read_text())['worlds'][0]['n_repeats']==10
assert json.loads((cfgdir/'compcancer_2x2.json').read_text())['worlds'][0]['n_repeats']==2
assert json.loads((cfgdir/'legacy8.json').read_text())['worlds'][0]['n_repeats']==10
ptb=json.loads((cfgdir/'ptb_fixed.json').read_text())['worlds'][0]
assert ptb['train_file']=='data/cohort_norm/01_train_data.arff'
assert ptb['test_file']=='data/cohort_norm/02_test_data.arff'

# No machine-specific paths in publication-default code/config/docs/scripts.
bad=[]
pat=[re.compile(r'/opt/2026/RLR'), re.compile(r'\b[A-Za-z]:[\\/]')]
for base in [ROOT/'R',ROOT/'configs',ROOT/'docs',ROOT/'scripts',ROOT/'README.md']:
    fs=[base] if base.is_file() else [p for p in base.rglob('*') if p.is_file()]
    for p in fs:
        if p.resolve()==Path(__file__).resolve(): continue
        txt=p.read_text(errors='ignore')
        if any(q.search(txt) for q in pat): bad.append(str(p.relative_to(ROOT)))
assert not bad, 'machine-specific paths: '+', '.join(sorted(set(bad)))

# Bundled GMT checksums.
for rel,expected in parse_hashes(ROOT/'checksums/MSIGDB_2026_1_HS_SHA256.txt').items():
    # Stage6 checksum file may contain original absolute path; match by basename.
    name=Path(rel).name
    p=ROOT/'reference_outputs/ptb_ora/gmt'/name
    assert p.is_file(), f'missing {p}'
    assert sha256(p)==expected, f'GMT checksum mismatch: {name}'

# Manifest invariants.
import csv
rows=list(csv.DictReader((ROOT/'paper_numbers_manifest_v3.csv').open(encoding='utf-8')))
assert len(rows)>=266
assert not any('external cohort split' in r.get('dataset_scope','').lower() for r in rows)
assert not any(r.get('status','').endswith('TODO') for r in rows)

if not a.quick:
    # Optional processed-data verification when files are present.
    missing=0; mismatch=0
    for rel,expected in parse_hashes(ROOT/'checksums/DATA_SHA256_EXPECTED.txt').items():
        p=ROOT/rel
        if not p.exists(): missing+=1; continue
        if sha256(p)!=expected:
            print('MISMATCH',rel); mismatch+=1
    assert mismatch==0, f'{mismatch} data checksum mismatches'
    print(f'Processed data not present (allowed): {missing}')

print('Release checks: PASS')
