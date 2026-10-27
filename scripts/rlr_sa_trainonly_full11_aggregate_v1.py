#!/usr/bin/env python3
"""Independently aggregate 11 SA-RLR training-only feasibility outcomes."""
from __future__ import annotations
import argparse,csv
from pathlib import Path
P=argparse.ArgumentParser()
P.add_argument("--input",required=True)
P.add_argument("--out",required=True)
opt=P.parse_args()
SETS=["alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
      "chen-2002","golub-1999-v1","gordon-2002",
      "pomeroy-2002-v1","shipp-2002-v1","singh-2002",
      "west-2001","yeoh-2002-v1"]
data=[]
for ds in SETS:
    files=list(Path(opt.input).rglob(f"{ds}/train_only_selection.csv"))
    if len(files)!=1:
        raise SystemExit(f"SA_TRAIN_AUDIT_STOP missing/duplicated dataset {ds} count {len(files)}")
    with files[0].open(encoding="utf-8",newline="") as h:
        rows=list(csv.DictReader(h))
    if len(rows)!=1 or rows[0]["dataset"]!=ds:
        raise SystemExit("SA_TRAIN_AUDIT_STOP dataset row unexpected")
    x=rows[0]
    for a in ("n_outer_T","n_outer_E_UNREAD","n_inner_D","n_inner_V","n_genes",
              "n_unique_proposed_ordinal_roots","n_roots_with_one_viable_branch",
              "n_roots_with_two_viable_branches"):
        if not x[a].isdigit() or int(x[a])<0:
            raise SystemExit("SA_TRAIN_AUDIT_STOP invalid count "+ds+":"+a)
    if int(x["n_inner_D"])+int(x["n_inner_V"])!=int(x["n_outer_T"]):
        raise SystemExit("SA_TRAIN_AUDIT_STOP D/V overlap or incomplete")
    if int(x["n_roots_with_one_viable_branch"])+int(x["n_roots_with_two_viable_branches"])>int(x["n_unique_proposed_ordinal_roots"]):
        raise SystemExit("SA_TRAIN_AUDIT_STOP root counts")
    if x["original_protocol_split"]!="rep00_fold00":
        raise SystemExit("SA_TRAIN_AUDIT_STOP wrong initial canonical split")
    if x["evidence_level"]!="INNER_TRAINING_ONLY_NOT_EXTERNAL_VALIDATION":
        raise SystemExit("SA_TRAIN_AUDIT_STOP false evidence label")
    state=x["selection_state"]
    if state not in ("ERROR","FALLBACK_NO_VALID_TREE","ACCEPTED_INNER_V","REJECTED_INNER_V"):
        raise SystemExit("SA_TRAIN_AUDIT_STOP unknown selection state")
    if x["accepted_V"] not in ("0","1") or (x["accepted_V"]=="1")!=(state=="ACCEPTED_INNER_V"):
        raise SystemExit("SA_TRAIN_AUDIT_STOP V gate inconsistency")
    if state=="ACCEPTED_INNER_V":
        if not (float(x["V_logloss_gain_vs_RLR50"])>=.005 and
                float(x["V_logloss_gain_vs_samepair_nontree"])>=.002 and
                float(x["V_BA_change"])>=-.02):
            raise SystemExit("SA_TRAIN_AUDIT_STOP accepted V gate violates thresholds")
    data.append(x)
output=Path(opt.out);output.mkdir(parents=True,exist_ok=True)
with (output/"all11_training_only_feasibility.csv").open("w",encoding="utf-8",newline="") as h:
    w=csv.DictWriter(h,fieldnames=data[0].keys())
    w.writeheader();w.writerows(data)
states={s:sum(z["selection_state"]==s for z in data) for s in
        ("ERROR","FALLBACK_NO_VALID_TREE","ACCEPTED_INNER_V","REJECTED_INNER_V")}
print("SA_TRAIN_11_DATASETS_INDEPENDENT_INVENTORY=PASS")
print("SA_TRAIN_CANONICAL_FIRST_SPLIT_AND_DV_DISJOINT_COUNTS=PASS")
print("SA_TRAIN_V_ACCEPTANCE_CRITERIA_REPLAYED=PASS")
print("SA_TRAIN_STATES "+str(states))
print("SA_TRAIN_TOTAL_PROPOSED_ROOTS="+str(sum(int(x["n_unique_proposed_ordinal_roots"]) for x in data)))
print("SA_TRAIN_TOTAL_ONE_VIABLE_BRANCH_ROOTS="+str(sum(int(x["n_roots_with_one_viable_branch"]) for x in data)))
print("SA_TRAIN_TOTAL_TWO_VIABLE_BRANCH_ROOTS="+str(sum(int(x["n_roots_with_two_viable_branches"]) for x in data)))
print("SA_TRAIN_NO_EXTERNAL_E_RESULTS_ORIGINAL_ARFF_PARITY_UNVERIFIED")
