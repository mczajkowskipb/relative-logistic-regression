#!/usr/bin/env python3
"""Read-only gate-versus-test diagnosis for all RLMT-DC-v0 accepted trees.

All 95 original outer folds are included; accepted models are printed in full.
This uses preserved development/confirmation metadata and already-audited
OUTER-TEST losses, never re-fits and NEVER proposes post-hoc new thresholds.
Does not write/modify result files.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OLD="results/rlmt_dc_published_data_5fold_v1"
PAIR="results/tspdt_recon_19x5_v1/four_method_comparison_v1/matched_95_folds.csv"


def require(ok,why):
    if not ok:
        raise SystemExit("DC_CONFIRMATION_AUDIT_STOP: "+why)


def csvrows(path: Path):
    require(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as fd:
        return list(csv.DictReader(fd))


def number(value,who):
    require(value is not None,"missing "+who)
    z=float(value)
    require(math.isfinite(z),"invalid numeric "+who)
    return z


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-repo",required=True)
    args=parser.parse_args()
    require(Path.cwd().resolve()==ROOT,"run in separate TSPDT worktree")
    ref=Path(args.reference_repo).resolve(strict=True)
    rroot=ref/OLD
    require((rroot/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
            "original 285-fit run incomplete")
    require((rroot/"published_5fold_audit_v1/independent_fold_metrics.csv").is_file()
            and (rroot/"published_5fold_accepted_replay_v1/accepted_models.csv").is_file()
            and (rroot/"published_5fold_summary_v1/descriptive_aggregates.json").is_file(),
            "original 285-fit independent audits missing")
    manifest=csvrows(rroot/"run_manifest.csv")
    folds=csvrows(ROOT/PAIR)
    require(len(manifest)==285 and len(folds)==95,
            "original benchmark/reporter records incomplete")
    lookup={}
    for row in folds:
        key=(row["world"],row["dataset"],row["split_id"])
        require(key not in lookup,"duplicate fold "+str(key))
        lookup[key]=row
    all_dc=[m for m in manifest if m["method"]=="RLMT_DC_v0"]
    require(len(all_dc)==95 and all(m["status"]=="OK" for m in all_dc),
            "95 successful DC model records required")
    accepted=[]
    rejected=[]
    for r in all_dc:
        key=(r["world"],r["dataset"],r["split_id"])
        require(key in lookup,"missing four-way matched comparator fold")
        folder=(ref/r["model_dir"]).resolve()
        require(folder.is_relative_to(rroot.resolve()),
                "model directory escapes frozen reference root")
        decision=json.loads((folder/"dc"/"decision.json").read_text())
        is_accept=r["accepted_split"].strip().lower()=="true"
        if not is_accept:
            require(r["accepted_split"].strip().lower()=="false",
                    "DC split missing Boolean value")
            rejected.append(key)
            continue
        require(decision["status"]=="ACCEPTED",
                "DC manifest accepts while metadata rejects "+str(key))
        confirm=decision["confirmation"]
        gain=number(confirm["gain"],"confirmation mean gain")
        se=number(confirm["se_gain"],"confirmation gain SE")
        lower=number(confirm["lower_gain"],"confirmation lower gain")
        n=int(confirm["n"])
        require(confirm["accepted"] is True and
                n>=20 and gain>0 and lower>0 and se>=0 and
                abs((gain-se)-lower)<1e-9,
                "original confirmation 1-SE criterion violated "+str(key))
        tree=json.loads((folder/"tree.json").read_text())
        require(tree["accepted_split"] is True and int(tree["size"])==2,
                "tree and confirmation metadata disagree "+str(key))
        nodes=csvrows(folder/"nodes.csv")
        require(len(nodes)==3 and nodes[1]["is_leaf"]=="TRUE" and
                nodes[2]["is_leaf"]=="TRUE",
                "unexpected two-leaf serialized node inventory "+str(key))
        sizes=(int(nodes[1]["n"]),int(nodes[2]["n"]))
        counts=(int(nodes[1]["class_0"]),int(nodes[1]["class_1"]),
                int(nodes[2]["class_0"]),int(nodes[2]["class_1"]))
        require(min(sizes)>=16 and min(counts)>=6,
                "final tree violates frozen child-support constraint "+str(key))
        row=lookup[key]
        delta=number(row["RLMT_DC_v0__delta_logloss_vs_RLR"],
                     "paired heldout logloss")
        brier=number(row["RLMT_DC_v0__delta_brier_vs_RLR"],
                     "paired heldout brier")
        accuracy=number(row["RLMT_DC_v0__delta_accuracy_vs_RLR"],
                        "paired heldout accuracy")
        require(int(r["n_test"])==int(row["n_test"]),
                "outer test count mismatch "+str(key))
        accepted.append(dict(
            world=key[0],dataset=key[1],fold=key[2],
            discovery_n=int(decision["discovery_n"]),
            confirmation_n=n,
            confirmation_mean_gain=gain,
            confirmation_se=se,
            confirmation_lower=gain-se,
            parent_validation_ll=number(confirm["parent_logloss"],
                                        "parent confirm loss"),
            split_validation_ll=number(confirm["split_logloss"],
                                       "split confirm loss"),
            final_left_n=sizes[0],final_right_n=sizes[1],
            heldout_n=int(row["n_test"]),
            delta_heldout_logloss=delta,
            delta_heldout_brier=brier,
            delta_heldout_accuracy=accuracy))
    require(len(accepted)==9 and len(rejected)==86,
            f"unexpected DC gate count: {len(accepted)}/95")
    print("=== NINE ACCEPTED RLMT-DC MODELS: TRAIN-ONLY CONFIRMATION vs OUTER TEST ===")
    print(f"{'dataset':22} {'fold':12} {'nV':>3} {'gainV':>10} "
          f"{'SEV':>10} {'lowerV':>10} {'dLL_outer':>12} {'dAcc':>8}")
    for row in sorted(accepted,key=lambda r:-r["delta_heldout_logloss"]):
        print(f"{row['dataset']:22} {row['fold']:12} "
              f"{row['confirmation_n']:3d} "
              f"{row['confirmation_mean_gain']:+10.6f} "
              f"{row['confirmation_se']:10.6f} "
              f"{row['confirmation_lower']:+10.6f} "
              f"{row['delta_heldout_logloss']:+12.6f} "
              f"{row['delta_heldout_accuracy']:+8.4f}")
    discordance=sum(row["delta_heldout_logloss"]>1e-10 for row in accepted)
    concordance=sum(row["delta_heldout_logloss"]< -1e-10 for row in accepted)
    print(f"DC_CONFIRMATION_ACCEPTED=9/95")
    print(f"DC_ACCEPTED_BUT_HARMFUL_ON_OUTER_TEST={discordance}/9")
    print(f"DC_ACCEPTED_AND_BENEFICIAL_ON_OUTER_TEST={concordance}/9")
    print("DC_CONFIRMATION_VS_TEST_AUDIT=PASS")
    print("NO_RETRAINING_NO_THRESHOLD_CHANGE=TRUE")


if __name__=="__main__":
    main()
