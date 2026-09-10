#!/usr/bin/env python3
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
manifest=ROOT/'checksums/DATA_SHA256_EXPECTED.txt'
bad=[]; missing=[]
for line in manifest.read_text().splitlines():
    if not line.strip(): continue
    expected,rel=line.split(None,1); p=ROOT/rel.strip()
    if not p.exists(): missing.append(rel.strip()); continue
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h!=expected: bad.append(rel.strip())
print('present_ok=', len(manifest.read_text().splitlines())-len(missing)-len(bad))
print('missing=',len(missing)); print('mismatch=',len(bad))
if bad:
    print('\n'.join('MISMATCH '+x for x in bad)); raise SystemExit(1)
