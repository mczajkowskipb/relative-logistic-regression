#!/usr/bin/env python3
"""Independent read-only verification of new Stage B against frozen original and CRIT7."""
from __future__ import annotations
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/crit_v7_stageB_mechanistic_v1"
FIX3_BASE = "results/crit_v7_original_rlr_2x3_v1"
METHODS = (
    "ORIGINAL_RLR", "GLOBAL_INTERACTION", "CRIT_TREE",
    "GATE_ONLY", "ADDITIVE_PAIRS", "UNCONDITIONAL_TREE",
)
DATASETS = ("chen-2002", "gordon-2002")
FOLDS = ("rep00_fold00", "rep00_fold01", "rep00_fold02")


def check(ok: bool, msg: str) -> None:
    if not ok:
        raise SystemExit("STAGEB_AUDIT_STOP: " + msg)


def rows(path: Path) -> list[dict]:
    check(path.is_file(), "missing file: " + str(path))
    with path.open(newline="", encoding="utf-8") as f:
        data = list(csv.DictReader(f))
    check(bool(data), "empty CSV: " + str(path))
    return data


def num(value: str) -> float:
    v = float(value)
    check(math.isfinite(v), "non-finite numeric value")
    return v


def logloss(y: int, p: float) -> float:
    p = max(1e-15, min(1 - 1e-15, p))
    return -(y * math.log(p) + (1 - y) * math.log1p(-p))


def main() -> None:
    check(len(sys.argv) == 3, "usage auditor.py ORIGINAL_REF FIX3_REF")
    ref = Path(sys.argv[1]).resolve(strict=True)
    fix3 = Path(sys.argv[2]).resolve(strict=True)
    archived = ref / "results/rlmt_dc_published_data_5fold_v1"
    prior = fix3 / FIX3_BASE
    check((archived / "status.txt").read_text().strip() == "RUNNER_EXIT=0",
          "archived original incomplete")
    check((prior / "status.txt").read_text().strip() ==
          "CRIT7_PILOT_BATCH_EXIT=0", "fix3 incomplete")
    manifest = rows(OUT / "comparison_manifest.csv")
    predictions = rows(OUT / "patient_predictions.csv")
    prior_manifest = rows(prior / "run_manifest.csv")
    original_manifest = rows(archived / "run_manifest.csv")
    check(len(manifest) == 36 and len(prior_manifest) == 24 and
          len(original_manifest) == 285, "incorrect frozen row count")
    keyset = {
        (dataset, fold, method)
        for dataset in DATASETS for fold in FOLDS for method in METHODS
    }
    keyrows = {(x["dataset"], x["split_id"], x["method"]): x
               for x in manifest}
    check(len(keyrows) == 36 and set(keyrows) == keyset,
          "incorrect six-fold by six-method matrix")
    groups = defaultdict(list)
    for r in predictions:
        key = (r["dataset"], r["split_id"], r["method"])
        check(key in keyrows, "extraneous patient prediction")
        groups[key].append(r)
    check(set(groups) == keyset and len(predictions) == sum(
        int(m["n_test"]) for m in manifest), "incomplete patient rows")

    prior_by_key = {(r["dataset"], r["split_id"], r["method"]):r
                    for r in prior_manifest}
    orig_by_key = {(r["dataset"], r["split_id"], r["method"]):r
                   for r in original_manifest}
    for dataset in DATASETS:
        for fold in FOLDS:
            old = orig_by_key[(dataset, fold, "RLR")]
            olddir = (ref / old["model_dir"]).resolve(strict=True)
            check(olddir.is_relative_to(archived.resolve(strict=True)),
                  "original archived model_dir escaped")
            original_preds = rows(olddir / "predictions_test.csv")
            patient_ids = [(r["sample_id"], r["y_true"]) for r in original_preds]
            original_probs = [num(r["prob"]) for r in original_preds]
            for method in METHODS:
                key = (dataset, fold, method)
                met = keyrows[key]
                group = groups[key]
                check(met["status"] == "OK" and met["stage"] ==
                      "EXPOSED_RETROSPECTIVE_DEV_ONLY",
                      "method status or scientific stage missing")
                check(len(group) == len(patient_ids) == int(met["n_test"]),
                      "unexpected outer test sample count")
                check([(r["sample_id"], r["y_true"]) for r in group] ==
                      patient_ids, "patient IDs or true labels disagree")
                check(len({r["sample_id"] for r in group}) == len(group),
                      "duplicate test IDs")
                check(max(abs(num(r["original_prob"]) - v)
                          for r,v in zip(group, original_probs)) < 1e-9,
                      "original patient probabilities differ")
                ps = [num(r["prob"]) for r in group]
                yy = [int(r["y_positive"]) for r in group]
                check(all(y in (0,1) for y in yy) and set(yy) == {0,1},
                      "invalid positive-class encoding")
                check(all(0<=p<=1 for p in ps), "probability outside [0,1]")
                losses = [logloss(y,p) for y,p in zip(yy,ps)]
                base_losses = [logloss(y,p) for y,p in zip(yy,original_probs)]
                check(abs(mean(losses)-num(met["logloss"]))<1e-9 and
                      abs(mean(base_losses)-num(met["original_logloss"]))<1e-9 and
                      abs(mean(losses)-mean(base_losses)-num(met["delta_vs_original"]))<1e-9 and
                      abs(mean([(p-y)**2 for y,p in zip(yy,ps)]) -
                          num(met["brier"]))<1e-9,
                      "recomputed metrics disagree")
                check(num(met["matched_glm_parity"])<1e-12,
                      "GLM parity gate failed")
                check(all(int(row["leaf_id"]) in (1,2) for row in group),
                      "expected nominated gate in every fold")
                if method in ("ORIGINAL_RLR","GLOBAL_INTERACTION","CRIT_TREE"):
                    prior_row = prior_by_key[key]
                    preddir = (fix3 / prior_row["model_dir"]).resolve(strict=True)
                    check(preddir.is_relative_to(prior.resolve(strict=True)),
                          "fix3 model_dir escaped")
                    saved = rows(preddir / "predictions_test.csv")
                    check([(r["sample_id"],r["y_true"]) for r in saved] ==
                          patient_ids, "prior patient ID/label mismatch")
                    check(max(abs(num(r["prob"])-p) for r,p in zip(saved,ps)) < 1e-9,
                          "frozen prior probabilities do not replay")
                if method == "CRIT_TREE" and met["tree_accepted"] == "0":
                    check(max(abs(p-q) for p,q in zip(ps,original_probs))<1e-12,
                          "rejected tree failed exact RLR fallback")
            # Report robustness of accepted trees against deletion of ONE patient.
            row = keyrows[(dataset,fold,"CRIT_TREE")]
            if row["tree_accepted"] == "1":
                group = groups[(dataset,fold,"CRIT_TREE")]
                yy = [int(p["y_positive"]) for p in group]
                original_probs = [num(r["original_prob"]) for r in group]
                ds = [logloss(y,num(r["prob"]))-logloss(y,p)
                      for y,r,p in zip(yy,group,original_probs)]
                leave_one = [(sum(ds)-v)/(len(ds)-1) for v in ds]
                wins = sum(v < -1e-10 for v in ds)
                losses = sum(v > 1e-10 for v in ds)
                leaves = defaultdict(list)
                for rowp,d in zip(group,ds):
                    leaves[int(rowp["leaf_id"])].append((int(rowp["y_positive"]),d))
                print(f"TREE_{dataset}_{fold}: n={len(ds)} improves={wins} "
                      f"worsens={losses} mean_delta={mean(ds):+.9f} "
                      f"leave_one_out_best={min(leave_one):+.9f} "
                      f"leave_one_out_worst={max(leave_one):+.9f}")
                for leaf, obs in sorted(leaves.items()):
                    print(f"  leaf_{leaf}: n={len(obs)} positive={sum(x[0] for x in obs)} "
                          f"delta={mean(x[1] for x in obs):+.9f}")

    print("=== STAGEB METHOD MEANS ACROSS SIX EXPOSED FOLDS ===")
    for method in METHODS:
        data = [keyrows[(d,f,method)] for d in DATASETS for f in FOLDS]
        ds = [num(x["delta_vs_original"]) for x in data]
        print(f"{method}: mean_delta={mean(ds):+.9f} "
              f"wins={sum(d < -1e-10 for d in ds)} "
              f"losses={sum(d > 1e-10 for d in ds)}")
    print("STAGEB_36_VECTORS_INDEPENDENTLY_AUDITED=PASS")
    print("STAGEB_EXISTING_METHOD_PATIENT_PARITY=PASS")
    print("STAGEB_SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")


if __name__ == "__main__":
    main()
