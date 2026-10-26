#!/usr/bin/env python3
"""Strict 11-dataset *original processed ARFF* byte verification.

Run on machine with the actual locally processed published inputs.
Does NOT download or reconstruct them. No patient data transferred.
"""
import argparse,hashlib,json
from pathlib import Path
SETS=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
def sha256(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for block in iter(lambda:f.read(1<<20),b""):
   h.update(block)
 return h.hexdigest()
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--root",default=".")
 parser.add_argument("--output",default="results/rlr_full11_exact_arff_gate_v1.json")
 args=parser.parse_args()
 root=Path(args.root)
 manifest=root/"checksums/DATA_SHA256_EXPECTED.txt"
 if not manifest.is_file():raise SystemExit("ARFF_GATE_STOP: original manifest missing")
 expected={}
 for line in manifest.read_text().splitlines():
  parts=line.split(maxsplit=1)
  if len(parts)==2:expected[parts[1]]=parts[0]
 result=[]
 for ds in SETS:
  shard="CompCancer_truth_2x2" if ds in SETS[-2:] else "CompCancer_truth_9x10"
  variants=(f"data/{shard}/{ds}.arff",f"data/CompCancer/{ds}.arff")
  ref=f"data/CompCancer/{ds}.arff"
  want=expected.get(ref)
  if not want:raise SystemExit(f"ARFF_GATE_STOP: expected hash absent for {ds}")
  p=next((root/f for f in variants if (root/f).is_file()),None)
  got=sha256(p) if p else None
  matched=(got==want)
  result.append(dict(dataset=ds,path=str(p) if p else None,
                     expected_sha256=want,found_sha256=got,
                     matched_original_published_ARFF=matched))
  print(f"ARFF_GATE {ds} status={'PASS' if matched else 'BLOCKED'}",flush=True)
 out=root/args.output
 out.parent.mkdir(parents=True,exist_ok=True)
 out.write_text(json.dumps(result,indent=2)+"\n")
 if not all(v["matched_original_published_ARFF"] for v in result):
  raise SystemExit("ARFF_GATE_STOP: EXACT_ORIGINAL_11_ARFF_DATASET_PARITY_NOT_ESTABLISHED")
 print("ARFF_GATE_EXACT_ORIGINAL_11_ARFF_SHA256=PASS")
if __name__=="__main__":main()
