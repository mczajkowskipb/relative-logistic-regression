#!/usr/bin/env python3
"""Fetch two publicly archived filtered CompCancer expression matrices.

Sources are public Schliep/de Souto supplementary matrices. This is
REUSED/exposed methodological development data, not fresh biomedical testing.
Downloader stores output only on CI runner; no secrets.
"""
import csv
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT=Path("results/rlr_sparse_leaf_real_pilot_v3")
DATA=ROOT/"source_data"
SETS={
 "chowdary-2006":("https://schlieplab.org/Static/Supplements/CompCancer/"
                  "Affymetrix/chowdary-2006/chowdary-2006_database.txt",
                  104,182),
 "singh-2002":("https://schlieplab.org/Static/Supplements/CompCancer/"
               "Affymetrix/singh-2002/singh-2002_database.txt",
               102,339)
}
def need(x,msg):
    if not x: raise SystemExit("SPARSE_DATA_STOP: "+msg)
def main():
    DATA.mkdir(parents=True,exist_ok=False)
    manifests=[]
    for name,(url,exp_samples,exp_genes) in SETS.items():
      r=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (research reproducibility)"})
      with urllib.request.urlopen(r,timeout=45) as f:
        raw=f.read()
      need(1024<len(raw)<25_000_000,"size invalid "+name)
      sha=hashlib.sha256(raw).hexdigest()
      lines=list(csv.reader(raw.decode("utf-8-sig").splitlines(),delimiter="\t"))
      need(len(lines)==exp_genes+2,
           f"unexpected source gene lines {name}: {len(lines)}")
      ids=[x.strip().strip('"') for x in lines[0] if x.strip()]
      labels=[x.strip().strip('"') for x in lines[1] if x.strip()]
      need(len(ids)==exp_samples and len(labels)==exp_samples,
           f"sample header/class inventory {name}: ids={len(ids)} classes={len(labels)}")
      need(len(set(ids))==exp_samples and len(set(labels))==2,
           "duplicate patient IDs or nonbinary phenotype")
      features=[]
      values=[]
      for row in lines[2:]:
        need(len(row)==exp_samples+1,
             f"not gene-name and {exp_samples} sample values in {name}")
        features.append(row[0].strip().strip('"'))
        values.append([float(v) for v in row[1:]])
      need(len(set(features))==exp_genes,"duplicate probe IDs")
      need(all(v==v and abs(v)<1e30 for col in values for v in col),
           "nonfinite or implausible expression")
      path=DATA/(name+".csv")
      with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["sample_id","class"]+features)
        for i in range(exp_samples):
          w.writerow([ids[i],labels[i]]+[values[g][i] for g in range(exp_genes)])
      manifests.append({"name":name,"original_url":url,
                        "download_sha256":sha,"samples":exp_samples,
                        "features":exp_genes,"classes":sorted(set(labels)),
                        "processed_csv":str(path)})
      print(f"SPARSE_SOURCE={name} samples={exp_samples} features={exp_genes} sha256={sha}",flush=True)
    (ROOT/"input_manifest.json").write_text(
        json.dumps(manifests,indent=2)+"\n",encoding="utf-8")
    print("SPARSE_PUBLIC_REAL_COMP_CANCER_DATA_2_OF_2=PASS")
    print("SPARSE_SOURCE_SCIENTIFIC_STAGE=REUSED_BIOMEDICAL_DEVELOPMENT_ONLY")
if __name__=="__main__":main()
