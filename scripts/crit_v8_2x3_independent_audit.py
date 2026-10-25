#!/usr/bin/env python3
"""CRIT-v8 independent archived patient-level replay (read-only)."""
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/crit_v8_one_shot_2x3_v1"
OLD = "results/rlmt_dc_published_data_5fold_v1"
DS = ("chen-2002", "gordon-2002")
FOLDS = ("rep00_fold00", "rep00_fold01", "rep00_fold02")
METHODS = ("ORIGINAL_RLR", "GLOBAL_CAL", "GATE_ONLY",
           "ADDITIVE_PAIRS", "GLOBAL_INTERACTION", "CONDITIONAL_TREE")


def need(ok: bool, why: str) -> None:
    if not ok:
        raise SystemExit("CRIT8_AUDIT_STOP: " + why)


def csvrows(p: Path) -> list[dict[str, str]]:
    need(p.is_file(), "missing " + str(p))
    with p.open(encoding="utf-8", newline="") as f:
        value = list(csv.DictReader(f))
    need(bool(value), "empty " + str(p))
    return value


def num(x: str) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError) as e:
        raise SystemExit("CRIT8_AUDIT_STOP: nonnumeric value") from e
    need(math.isfinite(v), "nonfinite value")
    return v


def metric(y: list[int], p: list[float]) -> tuple[float, float]:
    need(len(y) == len(p) > 0 and all(0 <= x <= 1 for x in p),
         "invalid probabilities")
    clipped = [max(1e-15, min(1 - 1e-15, v)) for v in p]
    loss = mean(-(t * math.log(v) + (1-t) * math.log1p(-v))
                for t, v in zip(y, clipped))
    brier = mean((t-v)**2 for t, v in zip(y, p))
    return loss, brier


def main() -> None:
    need(len(sys.argv) == 2, "usage auditor.py ORIGINAL_RLR_CHECKOUT")
    ref = Path(sys.argv[1]).resolve(strict=True)
    old = (ref / OLD).resolve(strict=True)
    need(Path.cwd().resolve() == ROOT, "must run within pinned v8 worktree")
    need((old / "status.txt").read_text().strip() == "RUNNER_EXIT=0",
         "incomplete frozen RLR reference")
    current = csvrows(OUT / "run_manifest.csv")
    original = csvrows(old / "run_manifest.csv")
    need(len(current) == 36 and len(original) == 285,
         "missing or duplicate run_manifest records")
    key = {(r["dataset"], r["split_id"], r["method"]): r for r in current}
    legacy = {(r["dataset"], r["split_id"], r["method"]): r
              for r in original}
    required = {(d,s,m) for d in DS for s in FOLDS for m in METHODS}
    need(len(key) == 36 and set(key) == required,
         "unexpected method-fold inventory")
    scores: dict[str,list[tuple[str,float,int]]] = {m:[] for m in METHODS}
    for dataset in DS:
        for fold in FOLDS:
            historic = legacy[(dataset,fold,"RLR")]
            model_dir = (ref / historic["model_dir"]).resolve(strict=True)
            need(model_dir.is_relative_to(old), "archive path traversal")
            oldrows = csvrows(model_dir / "predictions_test.csv")
            class_support = csvrows(model_dir.parent / "class_support.csv")
            need(len(class_support)==2, "class orientation unavailable")
            pos = class_support[1]["class"]
            sample_labels = [(r["sample_id"],r["y_true"]) for r in oldrows]
            need(len(set(sample_labels)) == len(oldrows),
                 "original test sample duplicates")
            y = [int(r["y_true"]==pos) for r in oldrows]
            oldprobs = [num(r["prob"]) for r in oldrows]
            oldll,oldbrier = metric(y,oldprobs)
            for method in METHODS:
                row = key[(dataset,fold,method)]
                need(row["status"]=="OK" and
                     row["scientific_stage"]=="EXPOSED_RETROSPECTIVE_DEV_ONLY",
                     "status or honest scientific designation missing")
                sub = (ROOT / row["model_dir"]).resolve(strict=True)
                need(sub.is_relative_to(OUT.resolve(strict=True)),
                     "pilot model path traversal")
                pr = csvrows(sub / "predictions_test.csv")
                need([(r["sample_id"],r["y_true"]) for r in pr] ==
                     sample_labels, "heldout patient IDs/labels mismatch")
                need(int(row["n_test"]) == len(pr) and
                     int(row["nomination_count"]) in (0,1),
                     "patient count or single nomination violated")
                p = [num(r["prob"]) for r in pr]
                root = [num(r["original_RLR_prob"]) for r in pr]
                need(max(abs(a-b) for a,b in zip(root,oldprobs))<1e-9,
                     "original paper baseline predictions differ")
                ll,br = metric(y,p)
                need(abs(ll-num(row["logloss"]))<1e-9 and
                     abs(br-num(row["brier"]))<1e-9 and
                     abs(oldll-num(row["original_logloss"]))<1e-9 and
                     abs(oldbrier-num(row["original_brier"]))<1e-9 and
                     abs(ll-oldll-num(row["delta_vs_original"]))<1e-9,
                     "saved metrics do not agree with patient replay")
                accepted = int(row["accepted"])
                available = int(row["validation_available"])
                need(accepted in (0,1) and available in (0,1),
                     "invalid acceptance flags")
                if method=="ORIGINAL_RLR":
                    need(accepted==0 and available==0,
                         "original model must not be V-gated")
                else:
                    if available:
                        lower = num(row["confirmation_lower"])
                        need((lower > 0.002) == bool(accepted),
                             "V decision not equal to predeclared one-shot gate")
                    else:
                        need(not accepted,
                             "unavailable correction was accepted")
                if method not in ("ORIGINAL_RLR","GLOBAL_CAL"):
                    need(available <= int(row["nomination_count"]),
                         "model available without D-nominated gate")
                if not accepted:
                    need(max(abs(a-b) for a,b in zip(p,root))<1e-12,
                         "unaccepted correction differs from exact RLR")
                leaves=[int(r["leaf_id"]) for r in pr]
                gates=[int(r["gate_side"]) for r in pr]
                need(all(z in (0,1,2) for z in gates),
                     "invalid physical gate IDs")
                need(all(z in (1,2) for z in leaves) if
                     method=="CONDITIONAL_TREE" and accepted
                     else all(z==0 for z in leaves),
                     "leaf-routing / status mismatch")
                if accepted and method!="GLOBAL_CAL":
                    need(row["gate_gene_a"] and row["gate_gene_b"] and
                         row["false_gene_a"] and row["false_gene_b"] and
                         row["true_gene_a"] and row["true_gene_b"],
                         "accepted relational correction without named features")
                scores[method].append((dataset,ll-oldll,accepted))
    for method in METHODS:
        parts = scores[method]
        need(len(parts)==6,"incomplete method scores")
        gain = [x[1] for x in parts]
        wins=sum(v < -1e-10 for v in gain)
        losses=sum(v > 1e-10 for v in gain)
        print(f"{method}: accepted={sum(x[2] for x in parts)}/6 "
              f"W/L/T={wins}/{losses}/{6-wins-losses} "
              f"mean_deltaLL={mean(gain):+.9f}")
        for ds in DS:
            print(f"  {ds}: mean_deltaLL="
                  f"{mean(x[1] for x in parts if x[0]==ds):+.9f}")
    print("CRIT8_ORIGINAL_RLR_PATIENT_PARITY_AUDIT=PASS")
    print("CRIT8_ONE_SHOT_36_VECTORS_AUDITED=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEV_ONLY")


if __name__ == "__main__":
    main()
