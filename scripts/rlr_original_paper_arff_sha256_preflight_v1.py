#!/usr/bin/env python3
"""Inventory ORIGINAL processed ARFF files referenced by frozen RLR paper.

Only source file paths/size/checksum are exported, never individual patients.
Use on author's Ubuntu project root (not the reconstructed public-source
GitHub Actions dataset directory). Exit 2 with --strict on missing inputs.
"""
import argparse
import csv
import hashlib
from pathlib import Path

P=argparse.ArgumentParser()
P.add_argument("--root",required=True,help="root containing original frozen data/ directory")
P.add_argument("--out",required=True,help="output directory for source manifest")
P.add_argument("--strict",action="store_true",help="exit nonzero if any original cohort world incomplete")
args=P.parse_args()
root=Path(args.root).expanduser().resolve()
out=Path(args.out).expanduser().resolve()
out.mkdir(parents=True,exist_ok=True)
worlds=[
 ("compcancer_original_9x10",root/"data"/"CompCancer_truth_9x10",9),
 ("compcancer_original_2x2",root/"data"/"CompCancer_truth_2x2",2),
 ("ptb_fixed",root/"data"/"cohort_norm",2),
 ("ptb_loco",root/"data"/"cohort_prelim",4),
 ("legacy8",root/"data"/"comparison",8)]
legacy={"adenocarcinoma","breast2","cns","colon","dlbcl","leukemia","ovarian","prostate"}
rows=[];coverage=[]
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  while True:
   data=f.read(1048576)
   if not data:break
   h.update(data)
 return h.hexdigest()
for name,folder,required in worlds:
 items=sorted([p for p in folder.rglob("*") if p.is_file() and p.suffix.lower()==".arff"]) if folder.is_dir() else []
 seen=set()
 for path in items:
  relative=path.relative_to(root).as_posix()
  normalized=path.stem.casefold().replace("_","-")
  seen.add(normalized)
  rows.append({"world":name,"file":relative,"size_bytes":path.stat().st_size,
               "sha256":digest(path)})
 if name=="ptb_fixed":
  complete=(folder/"01_train_data.arff").is_file() and (folder/"02_test_data.arff").is_file()
 elif name=="legacy8":
  complete=legacy.issubset({p.stem.casefold() for p in items})
 else:complete=len(items)>=required
 coverage.append({"world":name,"directory":str(folder),
                  "required_arff":required,"found_arff":len(items),
                  "complete":int(complete),
                  "note":"original processed ARFF; confirm original class/probe identity separately"})
with (out/"original_arff_sha256_manifest.csv").open("w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=["world","file","size_bytes","sha256"])
 w.writeheader();w.writerows(rows)
with (out/"original_paper_world_inventory.csv").open("w",encoding="utf-8",newline="") as f:
 w=csv.DictWriter(f,fieldnames=coverage[0].keys())
 w.writeheader();w.writerows(coverage)
for z in coverage:
 print("ORIGINAL_ARFF_WORLD",z["world"],"present",z["found_arff"],
       "required",z["required_arff"],"complete",bool(z["complete"]))
passed=all(z["complete"] for z in coverage)
print("ORIGINAL_RLR_ALL_PAPER_WORLD_ARFF_INVENTORY="+("PASS" if passed else "INCOMPLETE"))
print("ORIGINAL_RLR_RAW_PATIENT_DATA_NOT_EXPORTED=PASS")
if args.strict and not passed:raise SystemExit(2)
