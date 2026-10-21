#!/usr/bin/env python3
"""Independent read-only replay audit for ALL six synthetic CA-PS-v1 fits.

Recomputes standard classification and proper-score metrics from saved paired
heldout predictions. No R, no model fitting, no prior biomedical holdouts.
"""
import argparse
import csv
import json
import math
from pathlib import Path

SCENARIOS=("balanced_heterogeneous","minority_heterogeneous","balanced_null_tree")
SEEDS=(1139,2267)
NTEST=96
FIELDS=("accuracy","balanced_accuracy","logloss","brier")


def check(ok,msg):
    if not ok:
        raise SystemExit("RLMT_CA_PS_SYNTH_AUDIT_STOP: "+msg)


def readcsv(path):
    check(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def floatnum(v,label):
    try:
        x=float(v)
    except (ValueError,TypeError) as e:
        raise SystemExit("RLMT_CA_PS_SYNTH_AUDIT_STOP: invalid "+label) from e
    check(math.isfinite(x),"nonfinite "+label)
    return x


def score(rows):
    n=len(rows)
    check(n==NTEST,"wrong holdout length")
    y=[r["y_true"] for r in rows]
    check(set(y)=={"N","P"},"synthetic truth classes must be N/P")
    p=[floatnum(r["prob"],"probability") for r in rows]
    check(all(0<=v<=1 for v in p),"probability outside 0..1")
    yp=["P" if v>=0.5 else "N" for v in p]
    check(all(r["y_pred"]==t for r,t in zip(rows,yp)),
          "hard-class decision and probability mismatch")
    brier=sum((v-(int(t=="P")))**2 for t,v in zip(y,p))/n
    logloss=sum(-(int(t=="P")*math.log(max(1e-15,min(1-1e-15,v)))+
                  int(t=="N")*math.log1p(-max(1e-15,min(1-1e-15,v))))
                for t,v in zip(y,p))/n
    recall=lambda cls:sum(t==cls and z==cls for t,z in zip(y,yp))/sum(t==cls for t in y)
    return dict(accuracy=sum(t==z for t,z in zip(y,yp))/n,
                balanced_accuracy=(recall("N")+recall("P"))/2,
                logloss=logloss,brier=brier)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",required=True)
    args=ap.parse_args()
    root=Path(args.root).resolve(strict=True)
    project=Path(__file__).resolve().parent.parent
    check(root.is_relative_to((project/"results").resolve()),
          "outside isolated experimental result tree")
    p=json.loads((root/"suite_protocol.json").read_text(encoding="utf-8"))
    check(p["method"]=="RLMT_CA_PS_V1_K10_m10_d2_pool1000_L2",
          "wrong method identity")
    check(set(p["scenarios"])==set(SCENARIOS) and
          set(p["seeds"])==set(SEEDS),
          "simulation protocol/seed divergence")
    man=readcsv(root/"paired_six_scenarios.csv")
    aggr=readcsv(root/"scenario_means.csv")
    check(len(man)==6 and len(aggr)==3,
          "incomplete paired result grid")
    seen=set()
    saved={}
    for row in man:
        scenario=row["scenario"]
        seed=int(row["seed"])
        key=(scenario,seed)
        check(scenario in SCENARIOS and seed in SEEDS and key not in seen,
              "unexpected/duplicate synthetic result "+str(key))
        seen.add(key)
        check(int(row["n_train"])==240 and int(row["n_test"])==NTEST,
              "fixed sample sizes changed")
        folder=(root/"models"/scenario/f"seed_{seed}").resolve()
        check(folder.is_relative_to(root) and folder.is_dir(),
              "model path outside protocol")
        ca=readcsv(folder/"predictions_test.csv")
        rlr=readcsv(folder/"rlr_reference"/"predictions_test.csv")
        comps=readcsv(folder/"prediction_components_test.csv")
        check(len(ca)==len(rlr)==len(comps)==NTEST,
              "missing paired samples")
        for a,b,c in zip(ca,rlr,comps):
            check(a["sample_id"]==b["sample_id"]==c["sample_id"] and
                  a["y_true"]==b["y_true"],
                  "unpaired sample IDs or truth")
            pm=floatnum(a["prob"],"mixed")
            pr=floatnum(b["prob"],"original RLR")
            pc=floatnum(c["child_prob"],"child")
            parent=floatnum(c["parent_prob"],"parent")
            mixed=floatnum(c["mixed_prob"],"mixed saved")
            g=floatnum(c["shrink_weight"],"shrink weight")
            check(0<=g<1 and abs(parent-pr)<1e-9 and
                  abs(mixed-pm)<1e-9 and
                  abs(parent+g*(pc-parent)-pm)<1e-9,
                  "bad convex saved probability decomposition")
        d=json.loads((folder/"ca_ps"/"decision.json").read_text(encoding="utf-8"))
        accept=row["accepted"].lower()=="true"
        check(row["accepted"].lower() in ("true","false") and
              d["status"]==("ACCEPTED" if accept else "NO_SPLIT") and
              d["realized_leaves"]==(2 if accept else 1),
              "acceptance status mismatch")
        if not accept:
            check(all(floatnum(a["prob"],"ca")==
                      floatnum(b["prob"],"rlr")
                      for a,b in zip(ca,rlr)) and
                  all(floatnum(c["shrink_weight"],"g")==0 for c in comps),
                  "unsplit parent probabilities must be bit-exact")
        else:
            check((folder/"leaf_models/node_2/model_complete.rds").is_file() and
                  (folder/"leaf_models/node_3/model_complete.rds").is_file(),
                  "accepted model's leaf RDS absent")
        ref=score(rlr)
        new=score(ca)
        for met in FIELDS:
            for method,value in (("rlr",ref[met]),("ca",new[met])):
                check(abs(value-floatnum(row[f"{method}_{met}"],met))<1e-8,
                      "saved vs recomputed "+method+"/"+met+" "+str(key))
            check(abs((new[met]-ref[met])-
                      floatnum(row["delta_"+met],"delta "+met))<1e-8,
                  "paired metric difference mismatch "+str(key))
        saved[key]=row
    check(seen==set((s,t) for s in SCENARIOS for t in SEEDS),
          "missing a preregistered simulation case")
    for agg in aggr:
        scenario=agg["scenario"]
        check(scenario in SCENARIOS and int(agg["n"])==2,
              "bad scenario summary")
        group=[saved[(scenario,s)] for s in SEEDS]
        for metric in ("logloss","brier","balanced_accuracy"):
            expected=sum(floatnum(r["delta_"+metric],metric) for r in group)/2
            check(abs(expected-floatnum(agg["mean_delta_"+metric],metric))<1e-8,
                  "scenario means do not match six complete saved results")
        check(int(agg["accepted"])==
              sum(r["accepted"].lower()=="true" for r in group),
              "scenario acceptance count mismatch")
    print("RLMT_CA_PS_SYNTHETIC_AUDITED_FITS=6")
    print("RLMT_CA_PS_SYNTHETIC_DATASETS=3")
    print("RLMT_CA_PS_SYNTHETIC_INDEPENDENT_METRIC_REPLAY=PASS")
    print("NO_BIOMEDICAL_DATA=TRUE")


if __name__=="__main__":
    main()
