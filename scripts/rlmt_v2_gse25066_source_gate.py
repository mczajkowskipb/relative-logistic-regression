#!/usr/bin/env python3
"""Metadata/source-integrity ONLY for the externally nominated GSE25066 cohort.

Download original submitter GEO series matrices GSE25055 discovery and GSE25065
validation, create SHA256 + basic accession/metadata inventory. NO LABEL MAPPING,
NO MODEL FITTING, NO OUTCOME INSPECTION, NO use of previously viewed 19 cohorts.

Requires internet access to NCBI on operator's Ubuntu. Fail closed if incomplete.
NEVER overwrite original downloaded files, output dirs or old RLR results.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

GSES={"GSE25055":310,"GSE25065":198}
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn"
LIMIT=350_000_000

def check(ok,msg):
    if not ok:
        raise SystemExit("GEO_PROVENANCE_STOP: "+msg)

def download(url,path):
    request=urllib.request.Request(url,headers={
        "User-Agent":"rlmt-v2-independent-cohort-research/0.1"})
    with urllib.request.urlopen(request,timeout=40) as response:
        check(response.status==200,"HTTP status not 200")
        length=0
        with path.open("xb") as f:
            while True:
                chunk=response.read(1024*1024)
                if not chunk:break
                length+=len(chunk)
                check(length<=LIMIT,"GEO matrix download too large")
                f.write(chunk)
    check(length>1000,"GEO source suspiciously small")
    return length

def hashfile(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()

def summarize(gse,path,expected):
    fields={}
    marker=False
    rowcount=0
    with gzip.open(path,"rt",encoding="utf-8",errors="strict") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                marker=True
                break
            if line.startswith("!"):
                key,*values=next(csv.reader([line.rstrip("\n")],delimiter="\t"))
                fields.setdefault(key,[]).append(values)
        check(marker,"no GEO series matrix expression-table marker")
        header=next(csv.reader([next(f).rstrip("\n")],delimiter="\t"))
        check(header and header[0].strip('"')=="ID_REF",
              "expression matrix first column is not ID_REF")
        gsm=[v.strip('"') for v in header[1:]]
        check(len(gsm)==expected and len(set(gsm))==expected and
              all(x.startswith("GSM") and x[3:].isdigit() for x in gsm),
              "unexpected sample count or duplicated GSM")
        for line in f:
            if line.startswith("!series_matrix_table_end"):
                break
            rowcount+=1
        else:
            check(False,"no end of GEO expression matrix")
    sample_geo=fields.get("!Sample_geo_accession",[])
    check(len(sample_geo)==1 and
          [x.strip('"') for x in sample_geo[0]]==gsm,
          "series sample IDs do not match matrix columns")
    return {
      "accession":gse,"expected_samples":expected,
      "observed_samples":len(gsm),
      "unique_gsms":len(set(gsm)),
      "expression_probe_rows":rowcount,
      "source_format":"GEO submitter processed series matrix gzip",
      "sample_metadata_field_names":sorted(
        k for k in fields if k.startswith("!Sample_")),
      "characteristics_records":len(
        fields.get("!Sample_characteristics_ch1",[])),
      "label_definition_verified":False,
      "patient_group_identifiers_verified":False,
      "cohort_disjoint_patient_ids_verified":False,
      "preprocessing_scale_verified":False,
      "actual_license_terms_verified":False,
    },set(gsm)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output-dir",required=True)
    args=ap.parse_args()
    out=Path(args.output_dir).absolute()
    check(not out.exists(),"metadata output already exists")
    check(out.parent.is_dir(),"output parent must already exist")
    check(not out.is_symlink(),"symlink output is not allowed")
    out.mkdir()
    summary={}
    unique={}
    try:
        for gse,n in GSES.items():
            url=f"{BASE}/{gse}/matrix/{gse}_series_matrix.txt.gz"
            p=out/(gse+"_series_matrix.txt.gz")
            print("GEO_SOURCE="+url,flush=True)
            length=download(url,p)
            details,gsms=summarize(gse,p,n)
            details.update(source_url=url,compressed_bytes=length,
                           sha256=hashfile(p))
            unique[gse]=gsms
            summary[gse]=details
            print(f"GEO_SOURCE_VERIFIED={gse}:{n}:{details['sha256']}",
                  flush=True)
        overlap=unique["GSE25055"]&unique["GSE25065"]
        check(not overlap,
              "same GSM accession present in both discovery and validation")
        output={
          "scientific_stage":"SOURCE_AND_METADATA_ONLY_NO_OUTCOMES",
          "parent_series":"GSE25066",
          "discovery":"GSE25055",
          "validation":"GSE25065",
          "source":"NCBI GEO",
          "sample_id_overlap_n":len(overlap),
          "subject_id_overlap_NOT_verified":True,
          "response_definition_NOT_verified":True,
          "endpoints_NOT_LOADED":True,
          "no_training_or_prediction":True,
          "subseries":summary}
        (out/"PROVENANCE_METADATA_ONLY.json").write_text(
            json.dumps(output,indent=2,ensure_ascii=False)+"\n",
            encoding="utf-8")
        print("GEO_TRAINING_SUBSERIES=310")
        print("GEO_VALIDATION_SUBSERIES=198")
        print("GEO_GSM_ID_DUPLICATION=0")
        print("GEO_NEW_COHORT_METADATA_GATE=PASS")
        print("GEO_CLINICAL_LABEL_AND_GROUP_GATE=PENDING_MANUAL_AUDIT")
        print("NO_MODEL_FITTING=TRUE")
    except Exception:
        (out/"INCOMPLETE_DO_NOT_USE").write_text(
          "DO NOT USE: provenance download or integrity gate failed.\n")
        raise

if __name__=="__main__":
    main()
