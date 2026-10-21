#!/usr/bin/env python3
"""Immutable-results read-only follow-up: all 95 folds / 13 accepted CA-PS trees.

No refitting, no new splits, no rank tuning, no patient-level export.
Creates one NEW report directory, leaving old manifests/models/metrics intact.
This is posthoc EXPOSED cohort diagnostics, not method validation.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

HERE=Path(__file__).resolve().parent.parent
REL="results/rlmt_ca_ps_v1_retrospective_19x5"
ORIG="results/rlmt_dc_published_data_5fold_v1"
METHOD="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2"
FROZEN_CODE="e049fe20530ae7c5378db73f4e9186b5f83f789e"
REPORT="accepted_13_posthoc_diagnostic_v1"

def need(check,msg):
    if not check:
        raise SystemExit("CA_PS_POSTHOC_STOP: "+msg)

def rows(path):
    need(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        data=list(csv.DictReader(f))
    need(bool(data),"empty "+str(path))
    return data

def num(v,label):
    try:
        x=float(v)
    except (TypeError,ValueError) as e:
        raise SystemExit("CA_PS_POSTHOC_STOP: numeric "+label) from e
    need(math.isfinite(x),"nonfinite "+label)
    return x

def status_dir(path,expected):
    need((path/"status.txt").is_file(),"missing status "+str(path))
    need((path/"status.txt").read_text().strip()==expected,
         "incomplete source "+str(path))

def csvout(path,items):
    need(items,"no report data "+str(path))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]))
        w.writeheader()
        w.writerows(items)

def loss(p,y):
    q=min(1-1e-15,max(1e-15,p))
    return -(y*math.log(q)+(1-y)*math.log1p(-q))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-repo",required=True)
    opts=parser.parse_args()
    need(Path.cwd().resolve()==HERE,"run from new method worktree")
    reference=Path(opts.reference_repo).resolve(strict=True)
    new=(HERE/REL).resolve(strict=True)
    old=(reference/ORIG).resolve(strict=True)
    status_dir(new,"RLMT_CA_PS_RETRO_RUN_EXIT=0")
    status_dir(old,"RUNNER_EXIT=0")
    need((new/"CODE_SHA.txt").read_text().strip()==FROZEN_CODE,
         "experiment source is not the frozen v1 experiment")
    prior=(new/"independent_audit.log").read_text()
    need("RLMT_CA_PS_RETRO_INDEPENDENT_19X5_AUDIT=PASS" in prior,
         "original independent 380-model audit did not pass")
    indexed=rows(new/"four_method_retro_audit_v1"/"matched_95_folds.csv")
    original=rows(old/"run_manifest.csv")
    now=rows(new/"run_manifest.csv")
    need(len(indexed)==95 and len(now)==95 and len(original)==285,
         "model/fold inventory mismatch")
    original_bykey={(r["world"],r["dataset"],r["split_id"],r["method"]):r
                    for r in original}
    bykey={(r["world"],r["dataset"],r["split_id"]):r for r in now}
    need(len(original_bykey)==285 and len(bykey)==95,
         "duplicate source manifests")
    sample_totals=[]
    accepted=[]
    leaves=[]
    no_split=0
    all_folds=[]
    reason=Counter()
    for row in indexed:
        key=row["world"],row["dataset"],row["split_id"]
        need(key in bykey,"unpaired source fold "+str(key))
        nr=bykey[key]
        is_accepted=nr["accepted_split"].lower()=="true"
        need(nr["accepted_split"].lower() in ("true","false"),
             "bad decision manifest "+str(key))
        need(is_accepted==(row["accepted_split"].lower()=="true"),
             "manifest-vs-independent-audit acceptance mismatch")
        ca_dir=(HERE/nr["model_dir"]).resolve(strict=True)
        need(ca_dir.is_relative_to(new),"CA model path escapes result root")
        base_record=original_bykey[(*key,"RLR")]
        old_dc=original_bykey[(*key,"RLMT_DC_v0")]
        base_dir=(reference/base_record["model_dir"]).resolve(strict=True)
        dc_dir=(reference/old_dc["model_dir"]).resolve(strict=True)
        need(base_dir.is_relative_to(old) and dc_dir.is_relative_to(old),
             "archived model path escaped reference")
        class_support=rows(base_dir.parent/"class_support.csv")
        need(len(class_support)==2,"class levels unavailable")
        class0,class1=(r["class"] for r in class_support)
        model=rows(ca_dir/"predictions_test.csv")
        rlr=rows(base_dir/"predictions_test.csv")
        dc=rows(dc_dir/"predictions_test.csv")
        components=rows(ca_dir/"prediction_components_test.csv")
        n=len(model)
        need(n==int(nr["n_test"])==len(rlr)==len(dc)==len(components),
             "paired prediction length differs")
        d=json.loads((ca_dir/"ca_ps"/"decision.json").read_text())
        need(d["status"]==("ACCEPTED" if is_accepted else "NO_SPLIT"),
             "decision JSON diverges")
        reason[d.get("reason","unknown")]+=1
        total_delta=0.
        counts=defaultdict(lambda: dict(n=0,n1=0,ca_wrong=0,
             rlr_wrong=0,dc_wrong=0,true_prob_lt_001=0,
             true_prob_lt_005=0,raw_child_true_prob_lt_001=0,
             delta_ll_sum=0.,delta_brier_sum=0.,
             min_true_prob=1.,max_abs_prob_change=0.))
        per_sample=[]
        for a,b,c,mix in zip(model,rlr,dc,components):
            need(a["sample_id"]==b["sample_id"]==c["sample_id"]==mix["sample_id"] and
                 a["y_true"]==b["y_true"]==c["y_true"],
                 "test IDs or labels differ across methods")
            need(a["y_true"] in (class0,class1),"unknown true class")
            positive=int(a["y_true"]==class1)
            pa=num(a["prob"],"new prob")
            pb=num(b["prob"],"original RLR prob")
            pc=num(c["prob"],"old DC prob")
            parent=num(mix["parent_prob"],"component parent")
            child=num(mix["child_prob"],"component child")
            g=num(mix["shrink_weight"],"component g")
            l=int(a["leaf_id"])
            need(l in ((2,3) if is_accepted else (1,)),
                 "unexpected leaf ID")
            need(abs(parent-pb)<1e-12 and
                 abs(pa-(parent+g*(child-parent)))<1e-10,
                 "incorrect probability replay")
            local=counts[l]
            local["n"]+=1
            local["n1"]+=positive
            local["ca_wrong"]+=int((pa>=.5)!=bool(positive))
            local["rlr_wrong"]+=int((pb>=.5)!=bool(positive))
            local["dc_wrong"]+=int((pc>=.5)!=bool(positive))
            ptrue=pa if positive else (1-pa)
            child_true=child if positive else (1-child)
            local["true_prob_lt_001"]+=int(ptrue<.01)
            local["true_prob_lt_005"]+=int(ptrue<.05)
            local["raw_child_true_prob_lt_001"]+=int(child_true<.01)
            local["min_true_prob"]=min(local["min_true_prob"],ptrue)
            local["max_abs_prob_change"]=max(local["max_abs_prob_change"],
                                             abs(pa-pb))
            dl=loss(pa,positive)-loss(pb,positive)
            local["delta_ll_sum"]+=dl
            local["delta_brier_sum"]+=((pa-positive)**2-(pb-positive)**2)
            total_delta+=dl
            per_sample.append((l,g))
        delta=total_delta/n
        stored=num(row[METHOD+"__delta_logloss_vs_RLR"],
                   "independent fold delta")
        need(abs(delta-stored)<1e-9,"fold Δlogloss mismatch")
        dc_delta=num(row["RLMT_DC_v0__delta_logloss_vs_RLR"],
                     "old DC fold delta")
        if not is_accepted:
            no_split+=1
            need(abs(delta)<1e-12 and
                 len(counts)==1 and 1 in counts and
                 all(g==0 for _,g in per_sample),
                 "nonaccepted tree is not exact RLR fallback")
        else:
            need(len(counts)==2 and all(l in counts for l in (2,3)),
                 "accepted tree requires both observed test leaves")
            support=d.get("final_leaf_support")
            need(isinstance(support,dict),"missing final leaf class support")
            check_keys=("counts_true","counts_false","weight_true","weight_false")
            need(all(k in support for k in check_keys),
                 "incomplete leaf support/weights")
            expected_weights={}
            for leaf,side in ((2,"true"),(3,"false")):
                values=support["counts_"+side]
                need(len(values)==2,"invalid class support")
                n0,n1=map(int,values)
                total=n0+n1
                need(n0>=6 and n1>=6 and total>=16,
                     "original final per-class support invariant")
                balance=4*n0*n1/(total**2)
                expected=balance*total/(total+20)
                weight=num(support["weight_"+side],"trained leaf weight")
                need(abs(weight-expected)<1e-12,"support-based weight mismatch")
                expected_weights[leaf]=weight
                for routed,g in per_sample:
                    if routed==leaf:
                        need(abs(g-weight)<1e-12,"per-sample weight mismatch")
                z=counts[leaf]
                leaves.append(dict(
                    world=key[0],dataset=key[1],split_id=key[2],
                    leaf_id=leaf,train_n=total,train_class0=n0,train_class1=n1,
                    train_balance=balance,shrink_weight=weight,
                    test_n=z["n"],test_class1=z["n1"],
                    ca_wrong=z["ca_wrong"],rlr_wrong=z["rlr_wrong"],
                    dc_wrong=z["dc_wrong"],
                    ca_true_prob_lt_001=z["true_prob_lt_001"],
                    ca_true_prob_lt_005=z["true_prob_lt_005"],
                    raw_child_true_prob_lt_001=z["raw_child_true_prob_lt_001"],
                    min_ca_true_prob=z["min_true_prob"],
                    max_abs_ca_vs_rlr_prob_change=z["max_abs_prob_change"],
                    sum_delta_logloss_vs_rlr=z["delta_ll_sum"],
                    mean_delta_logloss_contribution=z["delta_ll_sum"]/n,
                    sum_delta_brier_vs_rlr=z["delta_brier_sum"]))
            need(abs(sum(counts[l]["delta_ll_sum"] for l in (2,3))/n-delta)<1e-12,
                 "leaf logloss decomposition doesn't sum to fold")
            confirmation=d.get("confirmation") or {}
            accepted.append(dict(
                world=key[0],dataset=key[1],split_id=key[2],
                n_test=n,confirmation_n=d.get("confirmation_n"),
                confirmation_gain=confirmation.get("gain"),
                confirmation_se=confirmation.get("se_gain"),
                confirmation_lower=confirmation.get("lower_gain"),
                gene_a=(d.get("nominated_pair") or {}).get("gene_a"),
                gene_b=(d.get("nominated_pair") or {}).get("gene_b"),
                ca_delta_logloss_vs_rlr=delta,
                old_dc_delta_logloss_vs_rlr=dc_delta,
                ca_delta_logloss_vs_old_dc=delta-dc_delta,
                ca_delta_accuracy_vs_rlr=num(
                    row[METHOD+"__delta_accuracy_vs_RLR"],"CA accuracy delta"),
                ca_delta_balanced_accuracy_vs_rlr=num(
                    row[METHOD+"__delta_balanced_accuracy_vs_RLR"],"CA BA delta"),
                ca_delta_brier_vs_rlr=num(
                    row[METHOD+"__delta_brier_vs_RLR"],"CA Brier delta"),
                ca_very_confident_wrong=sum(
                    counts[l]["true_prob_lt_001"] for l in (2,3)),
                rlr_wrong=sum(counts[l]["rlr_wrong"] for l in (2,3)),
                ca_wrong=sum(counts[l]["ca_wrong"] for l in (2,3))))
        all_folds.append((key,delta,is_accepted))
    need(len(accepted)==13 and len(leaves)==26 and no_split==82,
         "accepted/fallback/leaf model totals differ from validated 13/95")
    datasets=defaultdict(list)
    for key,dl,accepted_flag in all_folds:
        datasets[key[:2]].append((dl,accepted_flag))
    need(len(datasets)==19 and all(len(v)==5 for v in datasets.values()),
         "dataset x folds incomplete")
    posthoc_summary=[]
    for (world,dataset),cases in sorted(datasets.items()):
        posthoc_summary.append(dict(world=world,dataset=dataset,
            accepted=sum(flag for _,flag in cases),
            mean_delta_logloss_vs_rlr=statistics.mean(dl for dl,_ in cases),
            max_fold_delta_logloss=max(dl for dl,_ in cases),
            min_fold_delta_logloss=min(dl for dl,_ in cases)))
    globalmean=statistics.mean(row["mean_delta_logloss_vs_rlr"]
                               for row in posthoc_summary)
    oldmean=statistics.mean(
        num(row["RLMT_DC_v0__delta_mean_logloss_vs_RLR"],"old DC dataset mean")
        for row in rows(new/"four_method_retro_audit_v1"/"dataset_19_means.csv"))
    need(abs(globalmean-0.006177806)<1e-8 and
         abs(oldmean-0.017217934)<1e-8,
         "overall independent dataset-mean checkpoint changed")
    harms=sum(r["ca_delta_logloss_vs_rlr"]>1e-10 for r in accepted)
    benefits=sum(r["ca_delta_logloss_vs_rlr"]< -1e-10 for r in accepted)
    ties=13-harms-benefits
    accepted.sort(key=lambda r:-r["ca_delta_logloss_vs_rlr"])
    worst=accepted[:5]
    best=sorted(accepted,key=lambda r:r["ca_delta_logloss_vs_rlr"])[:5]
    report=new/REPORT
    need(not report.exists(),"diagnostic report already exists")
    report.mkdir()
    csvout(report/"accepted_13_folds.csv",accepted)
    csvout(report/"accepted_26_leaves.csv",leaves)
    csvout(report/"all_19_datasets.csv",posthoc_summary)
    summary=dict(
      scientific_stage="EXPOSED_RETROSPECTIVE_POSTHOC_DESCRIPTIVE_ONLY",
      analyzed_95_paired_folds=95,accepted_trees=13,accepted_leaves=26,
      no_split_parent_identical=82,
      accepted_outer_logloss_harm=harms,
      accepted_outer_logloss_benefit=benefits,
      accepted_outer_logloss_tie=ties,
      mean_dataset_delta_logloss_vs_rlr=globalmean,
      old_dc_mean_dataset_delta_logloss_vs_rlr=oldmean,
      accepted_ca_extreme_wrong_ptrue_under_001=sum(
        r["ca_true_prob_lt_001"] for r in leaves),
      old_model_files_touched=False,new_model_files_touched=False,
      training_or_tuning_performed=False,
      reason_counts=dict(reason))
    (report/"summary.json").write_text(json.dumps(
        summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("CA_PS_POSTHOC_REPLAY=PASS")
    print("CA_PS_ACCEPTED_TREES=13/95")
    print("CA_PS_ACCEPTED_LEAVES=26")
    print(f"ACCEPTED_LOGLOSS_BENEFIT/HARM/TIE={benefits}/{harms}/{ties}")
    print(f"EXCESS_WRONG_PTRUE_BELOW_001={summary['accepted_ca_extreme_wrong_ptrue_under_001']}")
    print(f"CA_DATASET_MEAN_DELTA_LOSS={globalmean:+.9f}")
    print(f"DC_DATASET_MEAN_DELTA_LOSS={oldmean:+.9f}")
    print("=== FIVE LARGEST HARMFUL ACCEPTED FOLDS (CA-vs-RLR) ===")
    for r in worst:
        print(f"{r['dataset']:20} {r['split_id']} "
              f"CA={r['ca_delta_logloss_vs_rlr']:+.6f} "
              f"DC={r['old_dc_delta_logloss_vs_rlr']:+.6f} "
              f"wrongCA/RLR={r['ca_wrong']}/{r['rlr_wrong']}")
    print("=== FIVE MOST FAVORABLE ACCEPTED FOLDS (CA-vs-RLR) ===")
    for r in best:
        print(f"{r['dataset']:20} {r['split_id']} "
              f"CA={r['ca_delta_logloss_vs_rlr']:+.6f} "
              f"DC={r['old_dc_delta_logloss_vs_rlr']:+.6f} "
              f"wrongCA/RLR={r['ca_wrong']}/{r['rlr_wrong']}")
    print("POSTHOC_NO_RETRAINING_NO_RETUNING=TRUE")
    print("SCIENTIFIC_STAGE=EXPOSED_RETROSPECTIVE_POSTHOC_ONLY")
    print("REPORT="+str(report))


if __name__=="__main__":
    main()
