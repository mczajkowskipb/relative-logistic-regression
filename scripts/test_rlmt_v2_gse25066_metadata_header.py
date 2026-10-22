#!/usr/bin/env python3
"""Artificial GEO header tests; NO real clinical response values inspected."""
from __future__ import annotations
import gzip
import importlib.util
import json
import tempfile
from pathlib import Path

here=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location(
    "geo_clinical_audit",here/"rlmt_v2_gse25066_clinical_metadata_audit.py")
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def make_head(where,gse,n,one_id=None,missing_clinical=False):
    samples=[f"GSM{700000+i}" for i in range(n)]
    ids=[f"{gse}_P{i}" for i in range(n)]
    if one_id is not None:
        ids[-1]=one_id
    labels=(["RD","pCR","NA"]*(n//3+1))[:n]
    lines=["!Sample_geo_accession\t"+'\t'.join(samples),
           "!Sample_characteristics_ch1\t"+
           '\t'.join('"sample id: '+v+'"' for v in ids)]
    if not missing_clinical:
        lines.append("!Sample_characteristics_ch1\t"+
            '\t'.join('"pathologic_response_pcr_rd: '+v+'"' for v in labels))
    lines.append("!series_matrix_table_begin")
    lines.append("ID_REF\t"+'\t'.join(samples))
    with gzip.open(where,"wt",encoding="utf-8") as f:
        f.write('\n'.join(lines)+'\n')
    return ids,labels

with tempfile.TemporaryDirectory(prefix="v2_fake_geo_test_") as p:
    root=Path(p)
    a=root/"discovery.gz"
    b=root/"validation.gz"
    ids_a,labels_a=make_head(a,"GSE25055",310)
    ids_b,labels_b=make_head(b,"GSE25065",198)
    meta_a,keys_a,gsms_a=mod.parse_discovery_and_id_fields("GSE25055",a)
    meta_b,keys_b,gsms_b=mod.parse_discovery_and_id_fields("GSE25065",b)
    assert meta_a["response_value_counts_discovery_only"]["pCR"]>0
    assert meta_a["response_value_counts_discovery_only"]["RD"]>0
    assert meta_a["response_value_counts_discovery_only"]["NA"]>0
    assert meta_b["response_value_counts_discovery_only"] is None
    assert meta_b["heldout_response_values_not_evaluated"] is True
    assert keys_a.isdisjoint(keys_b) and gsms_a==gsms_b.__class__(gsms_a)
    assert len(gsms_a)==310 and len(gsms_b)==198
    # Check no validation label values can appear in returned metadata.
    assert all(str(x) not in json.dumps(meta_b)
               for x in ("pCR","RD","NA"))
    bad=root/"bad.gz"
    make_head(bad,"GSE25055",310,missing_clinical=True)
    try:
        mod.parse_discovery_and_id_fields("GSE25055",bad)
    except SystemExit as exc:
        assert "response field incomplete" in str(exc)
    else:
        raise AssertionError("missing observed label field not detected")
    # Patient-proxy duplicate reports must be surfaced, not silently ignored.
    repeated=root/"repeated.gz"
    make_head(repeated,"GSE25055",310,one_id=ids_a[0])
    meta,*_=mod.parse_discovery_and_id_fields("GSE25055",repeated)
    assert meta["within_series_submitter_id_duplicates"]==1
print("RLMT_V2_GEO_SYNTHETIC_HEADER_UNIT=PASS")
print("VALIDATION_RESPONSE_VALUES_NOT_EVALUATED=TRUE")
print("NO_BIOMEDICAL_DATA=TRUE")
