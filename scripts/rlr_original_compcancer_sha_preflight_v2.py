#!/usr/bin/env python3
"""READ ONLY: verify all 11 original processed publication ARFF SHA256 hashes.
Never changes files, manifest, repository, or original sample data.
Run: python3 scripts/rlr_original_compcancer_sha_preflight_v2.py --data-root /PATH/TO/DATA
"""
import argparse,hashlib,pathlib,re,sys
p=argparse.ArgumentParser()
p.add_argument("--data-root",required=True,help="Parent of paper original ARFF folders; recursive .arff search only within this root")
p.add_argument("--expected",default="checksums/DATA_SHA256_EXPECTED.txt")
a=p.parse_args()
root=pathlib.Path(a.data_root).expanduser().resolve()
assert root.is_dir(),f"Data root is not a directory: {root}"
data=pathlib.Path(a.expected)
assert data.is_file(),f"SHA256 manifest missing. Run from RLR repo: {data}"
names=["alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
"chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
"shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1"]
expected={}
for line in data.read_text().splitlines():
 m=re.fullmatch(r"([0-9a-f]{64})\s+(.*)",line)
 if m:expected[m.group(2)]=m.group(1)
ok=missing=mismatch=duplicates=0
for name in names:
 key=f"data/CompCancer/{name}.arff"
 sha=expected.get(key)
 assert sha,"Missing original locked checksum: "+key
 hits=sorted(root.rglob(name+".arff"))
 if not hits:
  missing+=1
  print(f"{name}: MISSING")
  continue
 valid=[]
 for path in hits:
  digest=hashlib.sha256()
  with path.open("rb") as h:
   for block in iter(lambda:h.read(2**20),b""):digest.update(block)
  actual=digest.hexdigest()
  valid.append((path,actual))
 winners=[path for path,digest in valid if digest==sha]
 if winners:
  ok+=1
  if len(valid)>1:duplicates+=1
  print(f"{name}: PASS_SHA256  n_candidates={len(valid)}  path={winners[0]}")
 else:
  mismatch+=1
  print(f"{name}: MISMATCH  n_candidates={len(valid)}  expected={sha}  detected="+",".join(d for _,d in valid))
print(f"SHA_PREFLIGHT_PASS={ok}/11  MISSING={missing} MISMATCH={mismatch} DUPLICATE_NAMES={duplicates}")
print("ORIGINAL_PROCESSED_ARFF_PARITY="+("PASS" if ok==11 else "NOT_VERIFIED"))
sys.exit(0 if ok==11 else 2)
