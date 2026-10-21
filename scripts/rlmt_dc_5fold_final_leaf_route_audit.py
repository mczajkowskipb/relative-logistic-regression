#!/usr/bin/env python3
"""Read-only routed-leaf failure audit for all nine frozen RLMT-DC splits.

Reconstructs paired per-leaf logloss/Brier and test class mix from immutable
prediction CSVs and routing CSV; compares to final leaf training-class mix.
No model fitting, no selected-fold filtering, no output writes, no thresholds
for making future training decisions. This is post-hoc diagnosis only.
"""
from __future__ import annotations
import argparse
import collections
import csv
import json
import math
from pathlib import Path

REPO=Path(__file__).resolve().parent.parent
REF_RESULTS="results/rlmt_dc_published_data_5fold_v1"
FOLDS=REPO/"results/tspdt_recon_19x5_v1/four_method_comparison_v1/matched_95_folds.csv"


def need(condition,why):
    if not condition:
        raise SystemExit("DC_LEAF_DIAG_STOP: "+why)


def csvrows(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def num(x,why):
    try:
        v=float(x)
    except (ValueError,TypeError) as e:
        raise SystemExit("DC_LEAF_DIAG_STOP: not numeric "+why) from e
    need(math.isfinite(v),"nonfinite "+why)
    return v


def loss(y,p):
    p=min(1-1e-15,max(1e-15,p))
    return -(y*math.log(p)+(1-y)*math.log1p(-p))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference-repo",required=True)
    args=ap.parse_args()
    need(Path.cwd().resolve()==REPO,"run in separate TSPDT worktree")
    ref=Path(args.reference_repo).resolve(strict=True)
    results=(ref/REF_RESULTS).resolve(strict=True)
    need((results/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
         "full original 285 fits incomplete")
    need((results/"published_5fold_accepted_replay_v1/accepted_models.csv").is_file(),
         "independent two-leaf replay absent")
    folds=csvrows(FOLDS)
    need(len(folds)==95,"matched 95-fold report absent")
    compare={(r["world"],r["dataset"],r["split_id"]):r for r in folds}
    need(len(compare)==95,"duplicate fold in matched comparison")

    all_runs=csvrows(results/"run_manifest.csv")
    dc=[r for r in all_runs if r["method"]=="RLMT_DC_v0"]
    rlr={(r["world"],r["dataset"],r["split_id"]):r for r in all_runs
         if r["method"]=="RLR"}
    need(len(dc)==95 and len(rlr)==95 and
         all(r["status"]=="OK" for r in dc),
         "expected successful 95 DC and RLR folds")
    accepted=[r for r in dc if r["accepted_split"].lower()=="true"]
    need(len(accepted)==9,"expected exactly nine accepted DC splits")

    lines=[]
    for run in accepted:
        key=(run["world"],run["dataset"],run["split_id"])
        need(key in compare and key in rlr,"missing pair comparator "+str(key))
        folder=(ref/run["model_dir"]).resolve(strict=True)
        root=(ref/rlr[key]["model_dir"]).resolve(strict=True)
        need(folder.is_relative_to(results) and root.is_relative_to(results),
             "model artifact escapes immutable results")
        node=csvrows(folder/"nodes.csv")
        route=csvrows(folder/"test_routes.csv")
        dcp=csvrows(folder/"predictions_test.csv")
        rlp=csvrows(root/"predictions_test.csv")
        need(len(node)==3 and
             [x["node_id"] for x in node]==["1","2","3"],
             "wrong tree structure "+str(key))
        need(len(route)==len(dcp)==len(rlp)==int(run["n_test"]),
             "unpaired DC/RLR predictions "+str(key))
        need([x["sample_id"] for x in route] ==
             [x["sample_id"] for x in dcp] ==
             [x["sample_id"] for x in rlp] and
             [x["y_true"] for x in dcp] ==
             [x["y_true"] for x in rlp],
             "test ID/order/labels differ "+str(key))
        meta=json.loads((folder/"dc/decision.json").read_text())
        need(meta["status"]=="ACCEPTED","missing original positive confirmation gate")
        level=compare[key]["class1_reference_ARFF_order"]
        label0=compare[key]["class0_reference_ARFF_order"]
        for leaf in (2,3):
            training=node[leaf-1]
            subset=[i for i,x in enumerate(route)
                    if int(x["leaf_id"])==leaf]
            saved_dir=folder/"leaf_models"/f"node_{leaf}"
            active=csvrows(saved_dir/"active_pairs.csv")
            coeff=csvrows(saved_dir/"coefficients.csv")
            need(len(coeff)>=1,"empty saved final leaf coefficients "+str(key))
            count=int(training["n"])
            n0=int(training["class_0"])
            n1=int(training["class_1"])
            need(count==n0+n1 and n0>=6 and n1>=6 and count>=16,
                 "invalid final training class support "+str(key))
            scores=[]
            labels=[]
            dcwrong=0
            rlrwrong=0
            confident_wrong=0
            very_confident_wrong=0
            for ix in subset:
                rlrp=rlp[ix]
                dcpred=dcp[ix]
                pos = int(dcpred["y_true"]==level)
                need(dcpred["y_true"] in (level,label0),
                     "test class outside original ARFF levels "+str(key))
                pd=num(dcpred["prob"],"DC probability")
                pr=num(rlrp["prob"],"RLR probability")
                need(0<=pd<=1 and 0<=pr<=1,
                     "invalid probability range "+str(key))
                truth=dcpred["y_true"]
                labels.append(pos)
                dcwrong+=dcpred["y_pred"]!=truth
                rlrwrong+=rlrp["y_pred"]!=truth
                p_true=pd if pos else 1-pd
                confident_wrong+=p_true<0.05
                very_confident_wrong+=p_true<0.01
                scores.append((loss(pos,pd)-loss(pos,pr),
                               (pos-pd)**2-(pos-pr)**2))
            ntest=len(subset)
            sumloss=sum(t[0] for t in scores)
            sumbrier=sum(t[1] for t in scores)
            lines.append(dict(
                key=key,leaf=leaf,
                training_n=count,training_class0=n0,training_class1=n1,
                training_class1_frac=n1/count,
                test_n=ntest,test_class1=sum(labels),
                test_class1_frac=sum(labels)/ntest if ntest else math.nan,
                test_wrong_dc=dcwrong,test_wrong_rlr=rlrwrong,
                p_true_below_0_05=confident_wrong,
                p_true_below_0_01=very_confident_wrong,
                active_pairs=len(active),
                logloss_delta_sum=sumloss,
                logloss_delta_per_test= sumloss/ntest if ntest else math.nan,
                brier_delta_sum=sumbrier,
            ))
        group=[x for x in lines if x["key"]==key]
        need(len(group)==2,"missing one routed leaf "+str(key))
        observed=sum(x["logloss_delta_sum"] for x in group)/int(run["n_test"])
        expected=num(compare[key]["RLMT_DC_v0__delta_logloss_vs_RLR"],
                     "matched reference fold Δlogloss")
        need(abs(observed-expected)<1e-9,
             "leaf-wise test losses do not sum to frozen fold delta "+
             str(key))
        need(sum(x["test_n"] for x in group)==int(run["n_test"]),
             "test routing leaves do not partition fold "+str(key))

    need(len(lines)==18,"expected 18 terminal leaves from nine accepted trees")
    print("=== ALL ACCEPTED DC TREES: FINAL FIT TRAIN LABELS vs ROUTED HOLDOUT LABELS ===")
    print("Dataset/fold                Leaf   Train(n,1)  Test(n,1) WrongDC/RLR  "
          "Ptrue<.01  Active  ΣΔLL")
    for row in sorted(lines,key=lambda x:(-sum(
                 z["logloss_delta_sum"] for z in lines if z["key"]==x["key"]),
                 x["key"],x["leaf"])):
        k=row["key"]
        print(f"{k[1]:17} {k[2][-2:]:>2}       {row['leaf']:>2} "
              f" {row['training_n']:>3},{row['training_class1']:>3}"
              f"     {row['test_n']:>2},{row['test_class1']:>2}  "
              f" {row['test_wrong_dc']:>2}/{row['test_wrong_rlr']:<2}"
              f"         {row['p_true_below_0_01']:>2}"
              f"      {row['active_pairs']:>2}"
              f"    {row['logloss_delta_sum']:+9.4f}")
    print("DC_ACCEPTED_LEAVES=18")
    print("DC_ACCEPTED_TREES=9")
    print("DC_LEAF_ROUTE_OUTCOME_AUDIT=PASS")
    print("INFERENCE=POSTHOC_DESCRIPTIVE_ONLY")


if __name__=="__main__":
    main()
