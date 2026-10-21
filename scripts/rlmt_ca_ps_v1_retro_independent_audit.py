#!/usr/bin/env python3
"""Read-only-to-models independent four-method retrospective 95-fold audit.

Methods: original RLR, original RLMT-v1, original RLMT-DC-v0, NEW RLMT-CA-PS-v1.
Writes ONLY a fresh summary dir under NEW CA-PS result root.
This is EXPOSED retrospective method-development evidence, not external CV.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
OLD="results/rlmt_dc_published_data_5fold_v1"
NEW="results/rlmt_ca_ps_v1_retrospective_19x5"
METHODS=("RLR","RLMT_v1","RLMT_DC_v0","RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2")
METRICS=("accuracy","balanced_accuracy","macro_f1","logloss","brier")


def require(p,msg):
    if not p:
        raise SystemExit("RETRO_CA_AUDIT_STOP: "+msg)


def csvrows(path):
    require(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def numeric(v,desc):
    try: z=float(v)
    except (ValueError,TypeError) as e:
        raise SystemExit("RETRO_CA_AUDIT_STOP: invalid number "+desc) from e
    require(math.isfinite(z),"nonfinite "+desc)
    return z


def within(root,rel):
    p=(root/rel).resolve()
    require(p.is_relative_to(root.resolve()),
            "artifact escaped expected result directory")
    return p


def score(predict,class0,class1):
    require(len(predict)>0,"zero test predictions")
    labels=[r["y_true"] for r in predict]
    require(set(labels)=={class0,class1},
            "held-out fold missing a class or labels inconsistent")
    probs=[numeric(r["prob"],"test probability") for r in predict]
    require(all(0<=p<=1 for p in probs),"probability outside [0,1]")
    hard=[class1 if p>=0.5 else class0 for p in probs]
    require(all(r["y_pred"]==y for r,y in zip(predict,hard)),
            "recorded hard label differs from 0.5 decision")
    recalls=[]
    f1s=[]
    for c in (class0,class1):
        tp=sum(y==c and z==c for y,z in zip(labels,hard))
        fp=sum(y!=c and z==c for y,z in zip(labels,hard))
        fn=sum(y==c and z!=c for y,z in zip(labels,hard))
        recalls.append(tp/(tp+fn))
        f1s.append(2*tp/(2*tp+fp+fn) if (2*tp+fp+fn)>0 else 0)
    ys=[1 if y==class1 else 0 for y in labels]
    losses=[]
    for y,p in zip(ys,probs):
        q=min(1-1e-15,max(1e-15,p))
        losses.append(-(y*math.log(q)+(1-y)*math.log1p(-q)))
    n=len(ys)
    return dict(accuracy=sum(y==z for y,z in zip(labels,hard))/n,
                balanced_accuracy=sum(recalls)/2,
                macro_f1=sum(f1s)/2,
                logloss=statistics.mean(losses),
                brier=sum((y-p)**2 for y,p in zip(ys,probs))/n)


def identity(row):
    return row["world"],row["dataset"],row["split_id"]


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--reference-repo",required=True)
    args=a.parse_args()
    require(ROOT==Path.cwd().resolve(),
            "audit must run in separate new-method worktree")
    ref=Path(args.reference_repo).resolve(strict=True)
    old=(ref/OLD).resolve(strict=True)
    fresh=(ROOT/NEW).resolve(strict=True)
    for dest,expected in ((old,"RUNNER_EXIT=0"),
                          (fresh,"RLMT_CA_PS_RETRO_RUN_EXIT=0")):
        require((dest/"status.txt").read_text().strip()==expected,
                "incomplete experiment "+str(dest))
    oldrows=csvrows(old/"run_manifest.csv")
    newrows=csvrows(fresh/"run_manifest.csv")
    require(len(oldrows)==285 and len(newrows)==95,
            "old/new manifest dimensions not 285/95")
    oldlookup={}
    for r in oldrows:
        k=(*identity(r),r["method"])
        require(k not in oldlookup and r["status"]=="OK",
                "duplicate or invalid original fit "+str(k))
        oldlookup[k]=r
    newlookup={}
    for r in newrows:
        k=identity(r)
        require(k not in newlookup and r["status"]=="OK" and
                r["method"]==METHODS[-1],
                "duplicate or invalid new fit "+str(k))
        newlookup[k]=r
    require(len(newlookup)==95,"new fits not complete")
    cohorts=defaultdict(set)
    folds=[]
    all_models=[]
    pair_rows=[]
    for key in sorted(newlookup):
        world,dataset,fold=key
        cohorts[(world,dataset)].add(fold)
        row=newlookup[key]
        oldbaseline=oldlookup.get((*key,"RLR"))
        require(oldbaseline is not None,"missing frozen baseline "+str(key))
        original_base=within(ref,oldbaseline["model_dir"])
        support=csvrows(original_base.parent/"class_support.csv")
        require(len(support)==2,"original factor class support missing")
        class0,class1=(r["class"] for r in support)
        require(class0!=class1,"two distinct class levels required")
        paths={}
        for method in METHODS:
            record=(newlookup[key] if method==METHODS[-1] else
                    oldlookup.get((*key,method)))
            require(record is not None,"missing method "+str((*key,method)))
            root=ref if method!=METHODS[-1] else ROOT
            dst=old if method!=METHODS[-1] else fresh
            directory=within(root,record["model_dir"])
            require(directory.is_relative_to(dst) and directory.is_dir(),
                    "model path wrong "+str((*key,method)))
            preds=csvrows(directory/"predictions_test.csv")
            require(len(preds)==int(row["n_test"]),
                    "test dimensions differ "+str((*key,method)))
            paths[method]=(record,directory,preds)
        base_preds=paths["RLR"][2]
        ids=[r["sample_id"] for r in base_preds]
        truth=[r["y_true"] for r in base_preds]
        require(len(ids)==len(set(ids)),"test sample IDs not unique")
        scores={}
        for method in METHODS:
            record,directory,preds=paths[method]
            require([p["sample_id"] for p in preds]==ids and
                    [p["y_true"] for p in preds]==truth,
                    "unpaired test identities "+str((*key,method)))
            metrics=score(preds,class0,class1)
            scores[method]=metrics
            if method==METHODS[-1]:
                for m in METRICS:
                    require(abs(metrics[m]-numeric(record[m],m))<1e-9,
                            "new manifest score differs "+str((key,m)))
                decision=json.loads((directory/"ca_ps/decision.json").read_text())
                accepted=row["accepted_split"].lower()=="true"
                require(row["accepted_split"].lower() in ("true","false") and
                        decision["status"]==("ACCEPTED" if accepted else "NO_SPLIT"),
                        "confirmation/manifest status mismatch "+str(key))
                components=csvrows(directory/"prediction_components_test.csv")
                parent=csvrows(directory/"rlr_reference/predictions_test.csv")
                require(len(components)==len(parent)==len(preds),
                        "missing new components/parent "+str(key))
                for p,base,cmp,original in zip(preds,base_preds,components,parent):
                    require(p["sample_id"]==cmp["sample_id"]==
                            original["sample_id"]==base["sample_id"],
                            "parent/child sample order mismatch")
                    po=numeric(base["prob"],"archived RLR")
                    pn=numeric(original["prob"],"new original RLR")
                    pm=numeric(p["prob"],"CA-PS")
                    a1=numeric(cmp["parent_prob"],"mixture parent")
                    a2=numeric(cmp["child_prob"],"mixture child")
                    g=numeric(cmp["shrink_weight"],"mixture weight")
                    calc=a1+g*(a2-a1)
                    require(abs(po-pn)<1e-12 and abs(a1-po)<1e-12 and
                            0<=g<1 and abs(calc-pm)<1e-10 and
                            abs(numeric(cmp["mixed_prob"],"mixed")-pm)<1e-12,
                            "inconsistent parent RLR / convex mixture "+str(key))
                    if not accepted:
                        require(g==0 and pn==pm,
                                "no-split must equal RLR bitwise "+str(key))
                if accepted:
                    require((directory/"leaf_models/node_2/model_complete.rds").is_file() and
                            (directory/"leaf_models/node_3/model_complete.rds").is_file(),
                            "accepted split must preserve both fitted leaf RDS")
                pair_rows.append(dict(world=world,dataset=dataset,
                    split_id=fold,n_test=len(preds),accepted_split=accepted))
            all_models.append(dict(world=world,dataset=dataset,
                split_id=fold,method=method,n_test=len(preds),**metrics))
        foldrow=dict(world=world,dataset=dataset,split_id=fold,
                     accepted_split=row["accepted_split"].lower()=="true")
        for method in METHODS:
            for m in METRICS:
                foldrow[f"{method}__{m}"]=scores[method][m]
                foldrow[f"{method}__delta_{m}_vs_RLR"]=(
                    scores[method][m]-scores["RLR"][m])
        folds.append(foldrow)
    require(len(cohorts)==19 and all(v=={
        f"rep00_fold{i:02d}" for i in range(5)
    } for v in cohorts.values()),"not all 19 matched datasets x five folds")
    summaries=[]
    for world,dataset in sorted(cohorts):
        group=[r for r in folds
               if r["world"]==world and r["dataset"]==dataset]
        summary=dict(world=world,dataset=dataset,
                     n_folds=5,accepted_CA=sum(r["accepted_split"] for r in group))
        for method in METHODS:
            for metric in METRICS:
                v=statistics.mean(r[f"{method}__{metric}"] for r in group)
                baseline=statistics.mean(r[f"RLR__{metric}"] for r in group)
                summary[f"{method}__mean_{metric}"]=v
                summary[f"{method}__delta_mean_{metric}_vs_RLR"]=v-baseline
        summaries.append(summary)
    require(len(all_models)==380 and len(folds)==95 and
            len(summaries)==19,
            "final four-method analysis dimensions fail")
    ranking={}
    for method in METHODS[1:]:
        rank={}
        for m in METRICS:
            values=[r[f"{method}__delta_mean_{m}_vs_RLR"] for r in summaries]
            orient=-1 if m in ("logloss","brier") else 1
            rank[m]=dict(wins=sum(orient*v>1e-10 for v in values),
                         losses=sum(orient*v< -1e-10 for v in values),
                         ties=sum(abs(v)<=1e-10 for v in values),
                         mean_delta=statistics.mean(values))
        ranking[method]=rank
    out=fresh/"four_method_retro_audit_v1"
    require(not out.exists(),"read-only audit output already exists")
    out.mkdir()
    for fn,records in (("all_380_models.csv",all_models),
                       ("matched_95_folds.csv",folds),
                       ("dataset_19_means.csv",summaries)):
        with (out/fn).open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    (out/"descriptive_rankings.json").write_text(json.dumps(
        {"scientific_stage":"EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY",
         "not_prospective_or_independent_validation":True,
         "original_RLR_RLMT_fits_retrained":False,
         "new_CA_fits":95,"datasets":19,
         "four_methods":list(METHODS),
         "ranking_dataset_means":ranking},indent=2)+"\n")
    print("RETRO_CA_MODELS=95")
    print(f"RETRO_CA_ACCEPTED={sum(r['accepted_split'] for r in pair_rows)}/95")
    print("RETRO_FOUR_METHOD_FOLD_MODELS=380")
    print("RETRO_FOUR_METHOD_DATASETS=19")
    for method in METHODS[1:]:
        rank=ranking[method]
        print(f"{method}: ACC {rank['accuracy']['wins']}/"
              f"{rank['accuracy']['losses']}/"
              f"{rank['accuracy']['ties']} "
              f"mean_delta_logloss={rank['logloss']['mean_delta']:+.9f}")
    print("RLMT_CA_PS_RETRO_INDEPENDENT_19X5_AUDIT=PASS")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_DEVELOPMENT_ONLY")
    print("REPORT="+str(out))


if __name__=="__main__":
    main()
