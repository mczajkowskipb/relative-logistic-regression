#!/usr/bin/env python3
"""Read-only 5-fold exploratory summary, dataset as unit (NOT independent cohorts)."""
from __future__ import annotations

import collections
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/rlmt_dc_published_data_5fold_v1"
AUDIT = OUT / "published_5fold_audit_v1/independent_fold_metrics.csv"
REPLAY = OUT / "published_5fold_accepted_replay_v1/accepted_models.csv"
MANIFEST = OUT / "run_manifest.csv"
SUMMARY = OUT / "published_5fold_summary_v1"
METHODS = ("RLR","RLMT_v1","RLMT_DC_v0")
FOLD_IDS = tuple(f"rep00_fold{i:02d}" for i in range(5))


def need(test: bool, reason: str) -> None:
    if not test:
        raise SystemExit("FIVEFOLD_SUMMARY_STOP: " + reason)


def readcsv(path: Path) -> list[dict[str, str]]:
    need(path.is_file(), "missing prerequisite: " + str(path))
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def numeric(row: dict[str, str], field: str) -> float:
    v = float(row[field])
    need(math.isfinite(v), f"nonfinite {field}: {row}")
    return v


def mean(a: list[float]) -> float:
    need(len(a)>0,"empty mean")
    return sum(a)/len(a)


def median(a: list[float]) -> float:
    b = sorted(a)
    n = len(b)
    return b[n//2] if n%2 else (b[n//2-1]+b[n//2])/2


def main() -> None:
    need(Path.cwd().resolve()==ROOT,"run from repository root")
    need(not SUMMARY.exists(),"refusing overwrite of previous summary")
    need((OUT/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "original full runner not successful")
    paper_sha = (OUT/"PUBLICATION_MAIN_SHA.txt").read_text().strip()
    need(paper_sha=="4079a5ac400f1d2cfded523c7000f655470925ce",
         "original published RLR SHA differs")
    cfg = json.loads(
        (ROOT/"configs/rlmt_dc/published_data_5fold_v1.json").read_text())
    ordered = [(shard["world"],name)
               for shard in cfg["cohorts"] for name in shard["datasets"]]
    need(len(ordered)==19 and len(set(ordered))==19,
         "original 19 datasets registration mismatch")
    raw = readcsv(MANIFEST)
    audited = readcsv(AUDIT)
    accepted = readcsv(REPLAY)
    need(len(raw)==len(audited)==285,
         "need 285 completed and 285 independently audited model rows")
    need(all(r["status"]=="OK" for r in raw) and
         all(r["status"]=="OK" for r in audited),
         "a fit/audit error is present")
    fold_map = {}
    for r in audited:
        key = r["world"], r["dataset"], r["split_id"]
        need((key[:2] in ordered) and
             (key[2] in FOLD_IDS),
             "unexpected world/dataset/fold in independent audit")
        fold_map.setdefault(key,{})
        need(r["method"] in METHODS and
             r["method"] not in fold_map[key],
             "missing/duplicate comparator in a fold")
        fold_map[key][r["method"]] = r
    need(len(fold_map)==95 and
         all(set(m)==set(METHODS) for m in fold_map.values()),
         "not 95 complete three-method fold groups")

    manifest_map={}
    for r in raw:
        key=(r["world"],r["dataset"],r["split_id"],r["method"])
        need(key not in manifest_map,"duplicate model in manifest")
        manifest_map[key]=r
    replayed={(r["dataset"],r["split_id"],r["method"]) for r in accepted}
    actual_accepted={
        (r["dataset"],r["split_id"],r["method"])
        for r in raw if r["accepted_split"]=="TRUE"
    }
    need(len(actual_accepted)==len(replayed) and
         replayed==actual_accepted,
         "not every accepted two-leaf tree passed serialized replay")

    dataset_records = []
    fold_records = []
    for world,dataset in ordered:
        fold_triplets=[]
        for f in FOLD_IDS:
            key=(world,dataset,f)
            need(key in fold_map,"missing five-fold dataset:"+str(key))
            group=fold_map[key]
            row={}
            for method in METHODS:
                entry=group[method]
                original=manifest_map[(world,dataset,f,method)]
                need(original["status"]=="OK", "method failed")
                row[method]={
                    "logloss":numeric(entry,"logloss_recomputed"),
                    "brier":numeric(entry,"brier_recomputed"),
                    "leaves":int(original["realized_leaves"]),
                    "accepted":original["accepted_split"]=="TRUE",
                }
            fold_triplets.append(row)
            fold_records.append({
                "world":world, "dataset":dataset, "split_id":f,
                "rlr_logloss":row["RLR"]["logloss"],
                "v1_logloss":row["RLMT_v1"]["logloss"],
                "dc_logloss":row["RLMT_DC_v0"]["logloss"],
                "delta_v1_logloss":row["RLMT_v1"]["logloss"] -
                                   row["RLR"]["logloss"],
                "delta_dc_logloss":row["RLMT_DC_v0"]["logloss"] -
                                   row["RLR"]["logloss"],
                "delta_dc_brier":row["RLMT_DC_v0"]["brier"] -
                                 row["RLR"]["brier"],
                "accepted_v1":row["RLMT_v1"]["accepted"],
                "accepted_dc":row["RLMT_DC_v0"]["accepted"],
            })

        differences=[
            r["RLMT_DC_v0"]["logloss"]-r["RLR"]["logloss"]
            for r in fold_triplets
        ]
        v1_differences=[
            r["RLMT_v1"]["logloss"]-r["RLR"]["logloss"]
            for r in fold_triplets
        ]
        brier_differences=[
            r["RLMT_DC_v0"]["brier"]-r["RLR"]["brier"]
            for r in fold_triplets
        ]
        dataset_records.append({
            "world":world, "dataset":dataset,
            "n_folds":5,
            "rlr_mean_logloss":mean([
                r["RLR"]["logloss"] for r in fold_triplets]),
            "v1_mean_logloss":mean([
                r["RLMT_v1"]["logloss"] for r in fold_triplets]),
            "dc_mean_logloss":mean([
                r["RLMT_DC_v0"]["logloss"] for r in fold_triplets]),
            "mean_delta_v1_minus_RLR":mean(v1_differences),
            "mean_delta_dc_minus_RLR":mean(differences),
            "mean_delta_dc_brier_minus_RLR":mean(brier_differences),
            "dc_winning_folds":sum(z < -1e-10 for z in differences),
            "dc_losing_folds":sum(z > 1e-10 for z in differences),
            "dc_tied_folds":sum(abs(z)<=1e-10 for z in differences),
            "accepted_v1_folds":sum(
                z["RLMT_v1"]["accepted"] for z in fold_triplets),
            "accepted_dc_folds":sum(
                z["RLMT_DC_v0"]["accepted"] for z in fold_triplets),
        })

    need(len(dataset_records)==19 and len(fold_records)==95,
         "unexpected dataset/fold summary count")
    delta_dc=[d["mean_delta_dc_minus_RLR"] for d in dataset_records]
    delta_v1=[d["mean_delta_v1_minus_RLR"] for d in dataset_records]
    desc={
        "analysis":"RLMT-DC 19 publication datasets x 5 folds exploratory",
        "independent_validation":False,
        "n_datasets":19,"n_fold_triplets":95,"n_model_fits":285,
        "n_fits_with_error":0,
        "n_accepted_v1":sum(d["accepted_v1_folds"] for d in dataset_records),
        "n_accepted_dc":sum(d["accepted_dc_folds"] for d in dataset_records),
        "dataset_level_dc_vs_rlr":{
            "mean_delta":mean(delta_dc),"median_delta":median(delta_dc),
            "win":sum(z < -1e-10 for z in delta_dc),
            "loss":sum(z > 1e-10 for z in delta_dc),
            "tie":sum(abs(z)<=1e-10 for z in delta_dc)},
        "dataset_level_v1_vs_rlr":{
            "mean_delta":mean(delta_v1),"median_delta":median(delta_v1),
            "win":sum(z < -1e-10 for z in delta_v1),
            "loss":sum(z > 1e-10 for z in delta_v1),
            "tie":sum(abs(z)<=1e-10 for z in delta_v1)},
        "code_sha":(OUT/"CODE_SHA.txt").read_text().strip(),
        "cohort_unit":"historical dataset, not independent prospective cohort",
        "folds_are_not_independent_units":True,
        "test_based_hyperparameter_selection_permitted":False,
    }
    SUMMARY.mkdir(parents=True)
    for name,data in (
        ("per_dataset.csv",dataset_records),
        ("per_fold.csv",fold_records),
    ):
        with (SUMMARY/name).open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (SUMMARY/"descriptive_aggregates.json").write_text(
        json.dumps(desc,indent=2)+"\n")
    (SUMMARY/"README.txt").write_text(
        "Aggregates computed AFTER 285-fit full audit and replay.\n"
        "Five outer folds per historically exposed publication dataset.\n"
        "NO independent validation; NO inferential p-values; NO retuning.\n")
    print("=== PER-DATASET FIVE-FOLD LOGLOSS, DELTA_DC_MINUS_RLR ===")
    for d in dataset_records:
        print(f"{d['dataset']:24s} {d['mean_delta_dc_minus_RLR']:+.7f}"
              f"  accepted_dc={d['accepted_dc_folds']}/5")
    print("DATASET_LEVEL_DC_WIN_LOSS_TIE="
          f"{desc['dataset_level_dc_vs_rlr']['win']}/"
          f"{desc['dataset_level_dc_vs_rlr']['loss']}/"
          f"{desc['dataset_level_dc_vs_rlr']['tie']}")
    print(f"DATASET_LEVEL_MEAN_DC_DELTA={mean(delta_dc):+.8f}")
    print("TOTAL_ACCEPTED_DC="+str(desc["n_accepted_dc"]))
    print("TOTAL_ACCEPTED_V1="+str(desc["n_accepted_v1"]))
    print("PUBLISHED_5FOLD_SUMMARY=PASS")
    print("REPORT="+str(SUMMARY))


if __name__ == "__main__":
    main()
