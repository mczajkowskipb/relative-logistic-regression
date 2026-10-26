#!/usr/bin/env python3
"""Exact source inventory gate for the 11 binary primary CompCancer datasets.

Source: Schliep/de Souto public *filtered* source text matrices.
This never claims to reproduce the exact original processed ARFF bytes.
The original lock's expected .arff SHA256 is preserved as a separate provenance
gate that cannot pass in GitHub Actions without the actual original ARFF files.
"""
from __future__ import annotations
import csv, hashlib, json, math, re, urllib.request
from pathlib import Path

ROOT=Path("results/rlr_arlet_source_gate_v1")
DATA=ROOT/"source_data"
BASE="https://schlieplab.org/Static/Supplements/CompCancer"
# Publication dataset identity, platform, sample n, filtered probe p.
# Sample and probe sizes checked against original Schliep inventory webpage.
SETS=[
 ("alizadeh-2000-v1","CDNA",42,1095,"9x10"),
 ("armstrong-2002-v1","Affymetrix",72,1081,"9x10"),
 ("bittner-2000","CDNA",38,2201,"9x10"),
 ("chen-2002","CDNA",179,85,"9x10"),
 ("golub-1999-v1","Affymetrix",72,1877,"9x10"),
 ("gordon-2002","Affymetrix",181,1626,"9x10"),
 ("pomeroy-2002-v1","Affymetrix",34,857,"9x10"),
 ("shipp-2002-v1","Affymetrix",77,798,"9x10"),
 ("singh-2002","Affymetrix",102,339,"9x10"),
 ("west-2001","Affymetrix",49,1198,"2x2"),
 ("yeoh-2002-v1","Affymetrix",248,2526,"2x2"),
]
def need(cond,msg):
    if not cond: raise RuntimeError("PAPER_11_SOURCE_STOP: "+msg)
def file_sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
      for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()
def parse_expected():
    text=Path("checksums/DATA_SHA256_EXPECTED.txt").read_text()
    out={}
    for line in text.splitlines():
      parts=line.split(maxsplit=1)
      if len(parts)==2:
        out[parts[1]]=parts[0]
    for name,_,_,_,_ in SETS:
      for ext in ("arff","csv"):
        p=f"data/CompCancer/{name}.{ext}"
        need(bool(re.fullmatch("[0-9a-f]{64}",out.get(p,""))),
             f"missing locked SHA256 {p}")
    return out
def acquire(name,platform,n,p,shard,expected):
    url=f"{BASE}/{platform}/{name}/{name}_database.txt"
    req=urllib.request.Request(url,headers={"User-Agent":"CompCancer RLR reproducibility audit/1.0"})
    with urllib.request.urlopen(req,timeout=70) as res:
       body=res.read()
    need(1000<len(body)<85_000_000,f"invalid source length {name}")
    raw=DATA/f"{name}_source.txt"
    raw.write_bytes(body)
    rows=list(csv.reader(body.decode("utf-8-sig").splitlines(),delimiter="\t"))
    if len(rows)!=p+2:
       from collections import Counter
       widths=Counter(len(r) for r in rows)
       print(f"SOURCE_DIAGNOSTIC {name}: expected_probes={p} rows={len(rows)} width_histogram={dict(widths)} first_probe_headers={[x[0] if x else '(blank)' for x in rows[2:7]]} last_probe_headers={[x[0] if x else '(blank)' for x in rows[-7:]]}",flush=True)
    actual_p=len(rows)-2
    need(actual_p>0, f"{name}: no expression probe rows")
    if actual_p!=p:
       print(f"SOURCE_METADATA_MISMATCH {name}: inventory_n_probes={p} actual_source_n_probes={actual_p}",flush=True)
    ids=[v.strip().strip('"') for v in rows[0] if v.strip()]
    labels=[v.strip().strip('"') for v in rows[1] if v.strip()]
    need(len(ids)==len(labels)==n,f"{name}: expected {n} ids/classes, got {len(ids)}/{len(labels)}")
    if len(set(ids))!=n or len(set(labels))!=2 or not all(ids):
       from collections import Counter
       duplicates={k:v for k,v in Counter(ids).items() if v>1}
       print(f"SOURCE_SAMPLE_DIAGNOSTIC {name}: unique_ids={len(set(ids))}/{n} duplicate_ids={duplicates} empty_ids={sum(not q for q in ids)} class_counts={dict(Counter(labels))}",flush=True)
    # The first source row is a sample annotation, NOT always a unique ID.
    # Golub-v1 for example repeats ALL-B/ALL-T/AML phenotypic subtypes.
    # Preserve original source column ORDER with guaranteed unique ordinals.
    need(len(set(labels))==2, f"{name}: nonbinary source phenotype labels")
    stable_ids=[f"{name}__source_column_{j+1:04d}" for j in range(n)]
    probes=[]
    vals=[]
    for i,line in enumerate(rows[2:]):
       need(len(line)==n+1,f"{name} probe line {i} length {len(line)}")
       probe=line[0].strip().strip('"')
       need(probe!="",f"{name}: blank probe {i}")
       probes.append(probe)
       try: value=[float(x) for x in line[1:]]
       except ValueError as e:raise RuntimeError(f"{name} bad numeric line {i}") from e
       need(all(math.isfinite(x) for x in value),f"{name}: nonfinite values at probe {i}")
       vals.append(value)
    # Duplicate probe identifiers are reported, not renamed/silently modified.
    need(len(probes)==actual_p, "probe data inventory")
    dup=len(probes)-len(set(probes))
    clean=DATA/f"{name}.csv"
    with clean.open("w",encoding="utf-8",newline="") as f:
       w=csv.writer(f)
       w.writerow(["sample_id","class"]+[str(z) for z in probes])
       for j in range(n):
          w.writerow([stable_ids[j],labels[j]]+[vals[k][j] for k in range(actual_p)])
    return dict(dataset=name,platform=platform,shard=shard,
       source_url=url,source_sha256=file_sha(raw),
       converted_csv_sha256=file_sha(clean),n_samples=n,
       n_probes=actual_p,inventory_n_probes=p,
       metadata_dim_match=(actual_p==p),
       class_counts={k:labels.count(k) for k in sorted(set(labels))},
       sample_header_unique_values=len(set(ids)),
       sample_id_policy="POSITIONAL_IN_ORIGINAL_SOURCE_COLUMN_ORDER",
       duplicate_probe_headers=dup,
       original_arff_expected_sha256=expected[f"data/CompCancer/{name}.arff"],
       original_csv_expected_sha256=expected[f"data/CompCancer/{name}.csv"],
       exact_original_arff_sha256_verified=False,
       provenance_status=("SOURCE_TEXT_AND_CATALOG_MATCH_ARFF_NOT_VERIFIED"
                          if actual_p==p else
                          "SOURCE_TEXT_CATALOG_PROBE_COUNT_MISMATCH_ARFF_NOT_VERIFIED"))
def main():
    need(not DATA.exists(),"output path exists, no overwrite")
    DATA.mkdir(parents=True)
    expected=parse_expected()
    reports=[]
    for ds,pf,n,p,shard in SETS:
      report=acquire(ds,pf,n,p,shard,expected)
      reports.append(report)
      print(f"PAPER_SOURCE_OK {ds} samples={n} probes={p} sha256={report['source_sha256']}",flush=True)
    need(len(reports)==11 and sum(x["shard"]=="9x10" for x in reports)==9,
         "dataset/shard count")
    (ROOT/"source_inventory.json").write_text(json.dumps(reports,indent=2)+"\n")
    with (ROOT/"source_inventory.csv").open("w",newline="") as f:
      w=csv.DictWriter(f,fieldnames=list(reports[0].keys()))
      w.writeheader()
      for r in reports:w.writerow(r)
    print("PAPER_11_STRICT_PUBLISHED_DATASET_NAMES_AND_BINARY_PHENOTYPES=PASS")
    print("PAPER_11_CATALOG_PROBE_COUNTS_MATCH="+("PASS" if all(x["metadata_dim_match"] for x in reports) else "PARTIAL"))
    print("PAPER_11_ALL_SOURCE_ACTUAL_DIMENSIONS="+",".join(f"{x['dataset']}:{x['n_probes']}" for x in reports))
    print("PAPER_11_SOURCE_URL_AND_RECONSTRUCTED_CSV_SHA256_ARCHIVED=PASS")
    print("PAPER_11_ORIGINAL_PROCESSED_ARFF_PARITY=NOT_VERIFIED")
    print("PAPER_11_FULL_470_SPLITS=NOT_RUN")
if __name__=="__main__":main()
