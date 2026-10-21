#!/usr/bin/env python3
"""Read-only matched four-method report across 19 publication ARFF x five folds.

Use saved predictions only. No fitting, no hyperparameter search and no
claim of independent external validation. Requires already PASSed 285-fit
RLR/RLMT audits and the separately audited 95 TSPDT-inspired trees.

Usage: python3 scripts/tspdt_four_method_published_report.py \
    --reference-repo /opt/2026/RLR/relative-logistic-regression
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MAIN_SHA = "4079a5ac400f1d2cfded523c7000f655470925ce"
TSPDT_OUTPUT = ROOT / "results/tspdt_recon_19x5_v1"
REF_OUTPUT = "results/rlmt_dc_published_data_5fold_v1"
METHODS = ("RLR", "RLMT_v1", "RLMT_DC_v0", "TSPDT2011_INSPIRED_RECON_V1")
FIELDS = ("accuracy", "balanced_accuracy", "macro_f1", "logloss", "brier")
COMPARISON_DIR = TSPDT_OUTPUT / "four_method_comparison_v1"
FOLDS = tuple(f"rep00_fold{i:02d}" for i in range(5))


def require(ok: bool, why: str) -> None:
    if not ok:
        raise SystemExit("FOUR_METHOD_REPORT_STOP: " + why)


def csvrows(path: Path) -> list[dict[str, str]]:
    require(path.is_file(), "missing " + str(path))
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def asfloat(x: str, description: str) -> float:
    try:
        z = float(x)
    except (ValueError, TypeError) as e:
        raise SystemExit("FOUR_METHOD_REPORT_STOP: bad " + description) from e
    require(math.isfinite(z), "nonfinite " + description)
    return z


def score(predictions: list[dict[str, str]],
          expected_ids: list[str], expected_truth: list[str],
          levels: tuple[str, str], label: str) -> dict[str, float]:
    require(len(predictions) == len(expected_ids),
            "wrong prediction length for " + label)
    require([r["sample_id"] for r in predictions] == expected_ids,
            "test ID/order mismatch " + label)
    require([r["y_true"] for r in predictions] == expected_truth,
            "test truth mismatch " + label)
    predicted = [r["y_pred"] for r in predictions]
    require(all(z in levels for z in predicted),
            "out-of-vocabulary prediction " + label)
    probability = [asfloat(r["prob"], label + " probability")
                   for r in predictions]
    require(all(0.0 <= z <= 1.0 for z in probability),
            "invalid class probability " + label)
    n = len(expected_truth)
    accuracy = sum(t == p for t, p in zip(expected_truth,predicted)) / n
    recalls = []
    f1s = []
    for k in levels:
        n_true = sum(t == k for t in expected_truth)
        require(n_true > 0, "undefined test recall " + label + "/" + k)
        tp = sum(t == k and p == k
                 for t,p in zip(expected_truth,predicted))
        fp = sum(t != k and p == k
                 for t,p in zip(expected_truth,predicted))
        fn = n_true - tp
        recalls.append(tp/n_true)
        f1s.append(2*tp/(2*tp+fp+fn))
    binary = [int(t == levels[1]) for t in expected_truth]
    clip = [max(1e-15, min(1.0-1e-15, p)) for p in probability]
    logloss = -sum(y*math.log(p)+(1-y)*math.log1p(-p)
                   for y,p in zip(binary,clip))/n
    brier = sum((y-p)**2 for y,p in zip(binary,probability))/n
    return {
        "accuracy": accuracy,
        "balanced_accuracy": statistics.mean(recalls),
        "macro_f1": statistics.mean(f1s),
        "logloss": logloss, "brier": brier
    }



def infer_reference_class_order(reference: list[dict[str, str]],
                                observed_labels: list[str],
                                key: tuple[str, str, str]) -> tuple[str, str]:
    """Recover R's ARFF factor order from the audited RLR >=.5 rule.

    Sorting class labels alphabetically is WRONG: foreign::read.arff()
    retains the nominal ARFF level order. The independently audited
    original RLR branch uses class1 iff p >= .5; exactly one binary
    orientation must reproduce every reference hard-class prediction.
    """
    labels = set(observed_labels)
    require(len(labels) == 2, "not a binary test fold " + str(key))
    orders = []
    for pos in labels:
        neg = next(k for k in labels if k != pos)
        valid = all(
            row["y_pred"] == (
                pos if asfloat(row["prob"], "RLR probability "+str(key)) >= 0.5
                else neg
            )
            for row in reference
        )
        if valid:
            orders.append((neg, pos))
    require(len(orders) == 1,
            "cannot uniquely recover original R factor class order "
            "from independently audited RLR predictions "+str(key))
    return orders[0]

def ordered_map(rows: list[dict[str,str]],
                model_field: str) -> dict[tuple[str,str,str,str],dict[str,str]]:
    result = {}
    for row in rows:
        key = (row["world"],row["dataset"],row["split_id"],
               row[model_field])
        require(key not in result, "duplicate model " + str(key))
        result[key] = row
    return result


def safe_path(root: Path, relative: str, why: str) -> Path:
    path = (root / relative).resolve()
    require(path.is_relative_to(root.resolve()),
            "model artifact escapes result root: " + why)
    return path


def writecsv(path: Path, rows: list[dict]) -> None:
    require(bool(rows), "empty report " + str(path))
    with path.open("w",encoding="utf-8",newline="") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-repo",required=True)
    args = ap.parse_args()
    require(Path.cwd().resolve() == ROOT, "run from separate TSPDT worktree")
    require(not COMPARISON_DIR.exists(), "existing report: refuse overwrite")

    ref = Path(args.reference_repo).resolve(strict=True)
    refout = ref / REF_OUTPUT
    require(ref != ROOT, "reference and TSPDT worktrees cannot be the same")
    for out,name in ((TSPDT_OUTPUT,"TSPDT"),(refout,"RLR/RLMT")):
        require((out/"status.txt").is_file() and
                (out/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
                "experiment not completed successfully: " + name)
    for path in (
        refout/"published_5fold_audit_v1/independent_fold_metrics.csv",
        refout/"published_5fold_accepted_replay_v1/accepted_models.csv",
        refout/"published_5fold_firstfold_parity_v1/firstfold_prediction_parity.csv",
        refout/"published_5fold_summary_v1/per_dataset.csv",
        TSPDT_OUTPUT/"audit_v1/per_fold.csv",
        TSPDT_OUTPUT/"audit_v1/per_dataset.csv",
    ):
        require(path.is_file(), "missing passed prior audit " + str(path))
    require((refout/"PUBLICATION_MAIN_SHA.txt").read_text().strip()==MAIN_SHA
            and (TSPDT_OUTPUT/"PUBLICATION_MAIN_SHA.txt").read_text().strip()
            == MAIN_SHA, "not publication-frozen main")

    refinputs = csvrows(refout/"input_manifest.csv")
    tspinputs = csvrows(TSPDT_OUTPUT/"input_manifest.csv")
    require(len(refinputs)==len(tspinputs)==19,
            "expected identical 19-source inventories")
    for a,b in zip(refinputs,tspinputs):
        require((a["world"],a["dataset"],a["sha256"]) ==
                (b["world"],b["dataset"],b["sha256"]),
                "source differs between methods")

    refrows = csvrows(refout/"run_manifest.csv")
    tsp_rows = csvrows(TSPDT_OUTPUT/"run_manifest.csv")
    require(len(refrows)==285 and len(tsp_rows)==95,
            "need 285 ref and 95 reconstructed TSPDT model rows")
    require(all(r["status"]=="OK" for r in refrows+tsp_rows),
            "some original model fit failed")
    refmap=ordered_map(refrows,"method")
    tspmap=ordered_map(tsp_rows,"method")
    require(len(refmap)==285 and len(tspmap)==95,
            "duplicate recorded attempts")

    independent_ref = ordered_map(
        csvrows(refout/"published_5fold_audit_v1/independent_fold_metrics.csv"),
        "method")
    require(len(independent_ref)==285,
            "not 285 independently checked RLR/RLMT metrics")
    independent_tsp = {}
    for r in csvrows(TSPDT_OUTPUT/"audit_v1/per_fold.csv"):
        key=(r["world"],r["dataset"],r["split_id"])
        require(key not in independent_tsp, "duplicate TSPDT audit fold")
        independent_tsp[key] = r
    require(len(independent_tsp)==95,"not 95 independently checked TSPDT fits")

    all_models = []
    fold_comparisons = []
    dataset_comparisons = []
    for inputrow in refinputs:
        world=inputrow["world"]
        dataset=inputrow["dataset"]
        per_dataset=[]
        for fold in FOLDS:
            key=(world,dataset,fold)
            modelmap={}
            for method in METHODS:
                fullkey=(*key,method)
                source = tspmap if method==METHODS[-1] else refmap
                require(fullkey in source,
                        "missing model in paired fold "+str(fullkey))
                entry=source[fullkey]
                source_root=ROOT if method==METHODS[-1] else ref
                artifact=safe_path(source_root,entry["model_dir"],
                                   "/".join(fullkey))
                preds=csvrows(artifact/"predictions_test.csv")
                modelmap[method]=(entry,preds)
            reference=modelmap["RLR"][1]
            ids=[r["sample_id"] for r in reference]
            truths=[r["y_true"] for r in reference]
            require(len(ids)==len(set(ids))==len(truths),
                    "test ID duplication or dimension mismatch "+str(key))
            require(len(ids)==int(modelmap["RLR"][0]["n_test"]),
                    "wrong reference heldout size "+str(key))
            levels=infer_reference_class_order(reference,truths,key)

            foldresult={"world":world,"dataset":dataset,
                        "split_id":fold,"n_test":len(ids),
                        "class0_reference_ARFF_order":levels[0],
                        "class1_reference_ARFF_order":levels[1]}
            scores={}
            for method in METHODS:
                entry,pred=modelmap[method]
                require(int(entry["n_test"])==len(ids),
                        "different test counts: "+"/".join((*key,method)))
                sc=score(pred,ids,truths,levels,"/".join((*key,method)))
                scores[method]=sc
                if method=="TSPDT2011_INSPIRED_RECON_V1":
                    independent=independent_tsp[key]
                    require(abs(sc["logloss"] -
                                asfloat(independent["tspdt_logloss_laplace"],
                                        "audited TSPDT logloss"))<1e-10,
                            "TSPDT audit logloss not reproduced "+str(key))
                    for metric, col in (
                        ("accuracy","tspdt_accuracy"),
                        ("balanced_accuracy","tspdt_balanced_accuracy"),
                        ("macro_f1","tspdt_macro_f1"),
                        ("brier","tspdt_brier_laplace"),
                    ):
                        require(abs(sc[metric] -
                                    asfloat(independent[col],
                                            "audited TSPDT "+col))<1e-10,
                                "TSPDT audited metric mismatch "+
                                str((key,metric)))
                else:
                    independent=independent_ref[(*key,method)]
                    for metric, col in (
                        ("logloss","logloss_recomputed"),
                        ("brier","brier_recomputed"),
                    ):
                        require(abs(sc[metric] -
                                    asfloat(independent[col],
                                            "audited ref "+col))<1e-10,
                                "original ref audited metric mismatch "+
                                str((*key,method,metric))+
                                " ref_ARFF_class_order="+str(levels)+
                                " calculated="+format(sc[metric],".15g")+
                                " audited="+format(
                                    asfloat(independent[col],col),".15g"))
                all_models.append({
                    "world":world,"dataset":dataset,"split_id":fold,
                    "method":method,"n_test":len(ids),
                    "class0_reference_ARFF_order":levels[0],
                    "class1_reference_ARFF_order":levels[1],
                    **sc,
                    "leaf_count":entry.get("leaves","") if method==METHODS[-1]
                        else entry.get("realized_leaves",""),
                })
                for met,val in sc.items():
                    foldresult[f"{method}__{met}"]=val
                    foldresult[f"{method}__delta_{met}_vs_RLR"]=(
                        val-scores["RLR"][met]
                    )
            fold_comparisons.append(foldresult)
            per_dataset.append(foldresult)
        summary={"world":world,"dataset":dataset,"fold_count":5}
        for method in METHODS:
            for met in FIELDS:
                scorekey=f"{method}__{met}"
                vals=[v[scorekey] for v in per_dataset]
                avg=statistics.mean(vals)
                summary[f"{method}__mean_{met}"]=avg
                summary[f"{method}__delta_mean_{met}_vs_RLR"]=(
                    avg-statistics.mean([
                        v[f"RLR__{met}"] for v in per_dataset]))
        dataset_comparisons.append(summary)
    require(len(all_models)==380 and len(fold_comparisons)==95 and
            len(dataset_comparisons)==19,
            "four-method result dimensions invalid")
    rankings={}
    for method in METHODS[1:]:
        rankings[method]={}
        for met in FIELDS:
            # Lower is better for losses, higher for classification metrics.
            orientation= -1 if met in ("logloss","brier") else +1
            deltas=[r[f"{method}__delta_mean_{met}_vs_RLR"]
                    for r in dataset_comparisons]
            rankings[method][met]={
                "dataset_wins":sum(orientation*d>1e-10 for d in deltas),
                "dataset_losses":sum(orientation*d< -1e-10 for d in deltas),
                "dataset_ties":sum(abs(d)<=1e-10 for d in deltas),
                "mean_delta":statistics.mean(deltas),
                "median_delta":statistics.median(deltas),
            }
    COMPARISON_DIR.mkdir(parents=True)
    writecsv(COMPARISON_DIR/"all_380_model_fold_metrics.csv",all_models)
    writecsv(COMPARISON_DIR/"matched_95_folds.csv",fold_comparisons)
    writecsv(COMPARISON_DIR/"matched_19_dataset_means.csv",
             dataset_comparisons)
    (COMPARISON_DIR/"dataset_level_contrasts.json").write_text(
        json.dumps({
            "model_identity":"TSPDT-2011-inspired, not original TSPDT",
            "is_prospective_independent_validation":False,
            "unit_of_summary":"dataset average over 5 correlated outer folds",
            "n_datasets":19,"n_folds":95,"n_models":380,
            "classes_from":"saved original RLR heldout truth",
            "positive_class_recovered_from":"audited RLR p>=0.5 and saved y_pred; never alphabetically sorted",
            "probability_definition_TSPDT":"training leaf Laplace alpha=1",
            "comparators":rankings
        },indent=2)+"\n",encoding="utf-8")
    (COMPARISON_DIR/"README.txt").write_text(
        "Read-only common-metric comparison of 285 frozen RLR/RLMT "
        "fits and 95 independently replayed TSPDT-inspired fits.\n"
        "All metrics recomputed from prediction CSVs on exactly matched "
        "outer test IDs and true labels; no model fitting performed.\n"
        "R's original ARFF nominal factor order recovered from audited "
        "RLR predictions using the p>=0.5 hard-class rule; sorted "
        "string labels are never assumed to be the class order.\n"
        "TSPDT is reconstruction with train-only MAD f100, depth3, "
        "Wilson pruning, Laplace leaf scores, not original 2011 TSPDT.\n"
        "All 19 sets were already used in method development. "
        "Fold observations are not independent patients/datasets. "
        "No retrospective hyperparameter optimization or p-values.\n",
        encoding="utf-8")
    print("=== DATASET-LEVEL MEAN DIFFERENCE vs RLR (negative logloss better) ===")
    for method in METHODS[1:]:
        a=rankings[method]
        print(f"{method}: accuracy W/L/T="
              f"{a['accuracy']['dataset_wins']}/"
              f"{a['accuracy']['dataset_losses']}/"
              f"{a['accuracy']['dataset_ties']}; "
              f"balanced W/L/T={a['balanced_accuracy']['dataset_wins']}/"
              f"{a['balanced_accuracy']['dataset_losses']}/"
              f"{a['balanced_accuracy']['dataset_ties']}; "
              f"mean Δlogloss={a['logloss']['mean_delta']:+.6f}")
    print("COMPARISON_MODELS=380")
    print("COMPARISON_PAIRED_OUTER_FOLDS=95")
    print("COMPARISON_DATASETS=19")
    print("RLR_RLMT_TSPDT_FOUR_METHOD_REPORT=PASS")
    print("REPORT="+str(COMPARISON_DIR))


if __name__ == "__main__":
    main()
