#!/usr/bin/env python3
"""Clinical LABEL-SCHEMA AND ID audit for already-downloaded GSE25055/25065.

NO model fitting. Reads the original GEO gzip series-matrix HEADER ONLY, not
gene expression values. On DISCOVERY GSE25055, inspects observed pCR/RD
labels and reports only AGGREGATED COUNTS. On GSE25065 VALIDATION,
does NOT interpret/print/count/export response label VALUES -- only checks field PRESENCE.
Checks GSM uniqueness and submitting site's sample-id uniqueness across the
two subseries, but cannot prove clinical SUBJECT identity from GEO alone.
Writes only summary JSON, NO sample-level labels or patient IDs are exported.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

SHA = {
    "GSE25055": "9f8a94a9226f38d16380d776ab007c0453c3d31c9915d69562d5854baf3d6777",
    "GSE25065": "adac6710a8452be19e0f73070668b186a282ebeb34f9aaec09c884906034f6e9",
}
EXPECTED = {"GSE25055": 310, "GSE25065": 198}
SCHEMA_RESPONSE = "pathologic_response_pcr_rd"
SCHEMA_SAMPLE_ID = "sample id"
NA_VALUES = {"", "NA", "N/A", "NULL", "NOT AVAILABLE", "MISSING", "UNKNOWN"}

def need(ok, msg):
    if not ok:
        raise SystemExit("RLMT_V2_CLINICAL_AUDIT_STOP: " + msg)

def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1048576),b""):
            h.update(block)
    return h.hexdigest()

def parse_discovery_and_id_fields(gse,path):
    # Keep only minimal metadata needed for identity and endpoint audits.
    headers={}
    gsm=None
    by_sample=defaultdict(dict)
    clinical=defaultdict(dict)
    count_matrix=0
    with gzip.open(path,"rt",encoding="utf-8",errors="strict") as stream:
        for line in stream:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"): continue
            cells=next(csv.reader([line.rstrip("\n\r")],delimiter="\t"))
            field=cells[0]
            values=[v.strip() for v in cells[1:]]
            headers.setdefault(field,[]).append(len(values))
            if len(values)!=EXPECTED[gse]:
                need(False, f"metadata header dimension mismatch: {gse} {field}")
            if field=="!Sample_geo_accession":
                need(gsm is None,"multiple GEO accession header lines")
                gsm=values
            if field!="!Sample_characteristics_ch1": continue
            for i,entry in enumerate(values):
                if ":" not in entry: continue
                key, val=entry.split(":",1)
                key=key.strip().lower()
                val=val.strip()
                if key==SCHEMA_SAMPLE_ID:
                    need(i not in by_sample[gse],
                         "duplicate sample id field at one GSM")
                    by_sample[gse][i]=val
                elif key==SCHEMA_RESPONSE:
                    # Never access/count validation outcomes.
                    if gse=="GSE25055":
                        need(i not in clinical[gse],
                             "duplicate observed pCR/RD field")
                        clinical[gse][i]=val
                    else:
                        clinical[gse][i]="PRESENT_NOT_EVALUATED"
        else:
            need(False,f"no matrix begin delimiter {gse}")
        # Deliberately DO NOT read the expression matrix.
    need(gsm is not None and len(gsm)==EXPECTED[gse] and
         len(set(gsm))==EXPECTED[gse] and
         all(v.startswith("GSM") for v in gsm),
         "GSM count or uniqueness differs from source")
    need(len(by_sample[gse])==EXPECTED[gse],
         f"sample id field incomplete {gse}: {len(by_sample[gse])}")
    need(len(clinical[gse])==EXPECTED[gse],
         f"response field incomplete {gse}: {len(clinical[gse])}")
    ids=[by_sample[gse][i].strip() for i in range(EXPECTED[gse])]
    need(all(x and x.upper() not in NA_VALUES for x in ids),
         "missing submitter sample id")
    dup=len(ids)-len(set(ids))
    if gse=="GSE25055":
        observations=[clinical[gse][i].strip() for i in range(EXPECTED[gse])]
        normalized=[(v.upper() if v.upper() in NA_VALUES
                     else v) for v in observations]
        allowed={"pCR","RD"}.union(NA_VALUES)
        unexpected=sorted(set(normalized)-allowed)
        need(not unexpected,
             f"unexpected DISCOVERY observed endpoint categories: {unexpected}")
        counts=Counter("NA" if x in NA_VALUES else x for x in normalized)
    else:
        counts=None
    return {
        "gse":gse,
        "n_series_gsm":len(gsm),
        "within_series_submitter_id_duplicates":dup,
        "characteristic_rows":len(headers.get("!Sample_characteristics_ch1",[])),
        "observed_endpoint_field_count":len(clinical[gse]),
        "response_value_counts_discovery_only":dict(counts) if counts is not None else None,
        "heldout_response_values_not_evaluated":gse=="GSE25065",
    },set(ids),set(gsm)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir",required=True)
    ap.add_argument("--output-dir",required=True)
    a=ap.parse_args()
    source=Path(a.source_dir).resolve(strict=True)
    dest=Path(a.output_dir).absolute()
    need(dest.parent.is_dir() and not dest.exists(),
         "output must be a new, nonexisting directory under a known parent")
    manifests={}
    ids={}
    gsms={}
    for gse in ("GSE25055","GSE25065"):
        file=source/(gse+"_series_matrix.txt.gz")
        need(file.is_file(),"missing complete original source "+gse)
        need(digest(file)==SHA[gse],"GEO source SHA256 changed "+gse)
        meta, samples, accessions=parse_discovery_and_id_fields(gse,file)
        manifests[gse]=meta
        ids[gse]=samples
        gsms[gse]=accessions
    between=len(ids["GSE25055"].intersection(ids["GSE25065"]))
    gsm_overlap=len(gsms["GSE25055"].intersection(gsms["GSE25065"]))
    need(gsm_overlap==0,"GSM accessions overlap across GEO subseries")
    result={
        "scientific_stage":"DISCOVERY_LABEL_COUNTS_AND_VALIDATION_ID_ONLY",
        "data_integrity":"SHA256_MATCHES_STAGE1_ARCHIVE",
        "study":"GSE25066",
        "endpoint_candidate":"observed pathologic_response_pcr_rd",
        "binary_observed_classes":["pCR","RD"],
        "not_acceptable_truth":["dlda30_prediction","rcb_0_i_prediction",
                                "chemosensitivity_prediction"],
        "no_model_fit":True,
        "validation_response_values_not_evaluated_or_exported":True,
        "clinical_subject_independence_not_proven":True,
        "within_series_GSM_overlaps":gsm_overlap,
        "between_submitter_sample_id_overlap":between,
        "discovery":manifests["GSE25055"],
        "validation":manifests["GSE25065"],
        "final_clinical_group_gate":"BLOCKED_PENDING_SOURCE_VERIFICATION",
        "validation_label_values_not_reported_or_used":True,
    }
    dest.mkdir()
    (dest/"CLINICAL_SCHEMA_DISCOVERY_ONLY.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",
        encoding="utf-8")
    print("GSE25066_GSM_DISJOINT=PASS")
    print("GSE25066_DISCOVERY_OBSERVED_PCR_RD="+
          json.dumps(manifests["GSE25055"]["response_value_counts_discovery_only"],
                     sort_keys=True))
    print(f"GSE25066_SUBMITTER_SAMPLE_ID_OVERLAP={between}")
    print("GSE25066_VALIDATION_RESPONSE_VALUES_NOT_EVALUATED=TRUE")
    print("GSE25066_CLINICAL_SCHEMA_AUDIT=PASS")
    print("GSE25066_SUBJECT_INDEPENDENCE=NOT_PROVEN")
    print("GSE25066_CLINICAL_MODEL_FITS=0")
    print("REPORT="+str(dest/"CLINICAL_SCHEMA_DISCOVERY_ONLY.json"))

if __name__=="__main__":
    main()
