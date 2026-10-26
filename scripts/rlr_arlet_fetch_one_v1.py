#!/usr/bin/env python3
"""Download ONE of precisely 11 manuscript CompCancer public source matrices."""
import argparse, importlib.util, json
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument("--dataset",required=True)
args=parser.parse_args()
spec=importlib.util.spec_from_file_location("source_gate",
    "scripts/rlr_arlet_source_gate_v1.py")
import sys
m=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=m
spec.loader.exec_module(m)
rows=[s for s in m.SETS if s[0]==args.dataset]
m.need(len(rows)==1,"dataset not in canonical 11")
name,platform,n,p,shard=rows[0]
m.ROOT=Path("results/rlr_arlet_v1")/name
m.DATA=m.ROOT/"source_data"
m.need(not m.DATA.exists(),"refusing overwrite of existing source matrix")
m.DATA.mkdir(parents=True)
e=m.parse_expected()
inventory=m.acquire(name,platform,n,p,shard,e)
(m.ROOT/"input_source_manifest.json").write_text(json.dumps(inventory,indent=2)+"\n")
print("RLR_11_ONE_SOURCE_DATASET_VERIFIED="+name)
print("RLR_11_ONE_SOURCE_SAMPLE_COUNT="+str(n))
print("RLR_11_ONE_SOURCE_ACTUAL_FEATURES="+str(inventory["n_probes"]))
print("RLR_11_ONE_SOURCE_DUP_PROBE_HEADERS="+str(inventory["duplicate_probe_headers"]))
print("RLR_11_ONE_PROCESSED_ARFF_SHA_PARITY=NOT_VERIFIED")
