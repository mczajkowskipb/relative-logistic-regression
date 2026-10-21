#!/usr/bin/env python3
"""Compare all first folds in new five-fold run against frozen 57-fit pilot.

Both sides must have been independently audited; no model fitting.
"""
from __future__ import annotations
import csv
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
A=ROOT/"results/rlmt_dc_published_data_smoke_v1"
B=ROOT/"results/rlmt_dc_published_data_5fold_v1"
EXPECTED_OLD_SHA="814fee71fe6306db1e62b2ca47b1038df82d77a8"
RESULT=B/"published_5fold_firstfold_parity_v1"


def require(ok,msg):
    if not ok: raise SystemExit("FIVEFOLD_PILOT_PARITY_STOP: "+msg)


def read(path):
    require(path.is_file(),"missing "+str(path))
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    require(Path.cwd().resolve()==ROOT,"run from repo root")
    require(not RESULT.exists(),"will not overwrite prior parity audit")
    require((A/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
            "prior pilot did not complete")
    require((B/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
            "new five-fold run did not complete")
    require((A/"CODE_SHA.txt").read_text().strip()==EXPECTED_OLD_SHA,
            "unexpected first-fold pilot model SHA")
    require((B/"published_5fold_audit_v1/independent_fold_metrics.csv").is_file(),
            "independent five-fold audit absent")
    require((A/"published_audit_v1/independent_fold_metrics.csv").is_file(),
            "independent first-fold audit absent")
    originals=read(A/"run_manifest.csv")
    new=read(B/"run_manifest.csv")
    require(len(originals)==57 and len(new)==285,
            "unexpected manifest sizes")
    original_map={
      (x["world"],x["dataset"],x["split_id"],x["method"]):x
      for x in originals}
    require(len(original_map)==57,"duplicate original entries")
    new_map={
      (x["world"],x["dataset"],x["split_id"],x["method"]):x
      for x in new}
    require(len(new_map)==285,"duplicate five-fold entries")
    report=[]
    for key,old in original_map.items():
        require(key[2]=="rep00_fold00","unexpected old fold identifier")
        require(key in new_map,"missing first fold in 5-fold run")
        fresh=new_map[key]
        require(old["status"]==fresh["status"]=="OK",
                f"fit failure at {key}")
        op=read(ROOT/old["model_dir"]/"predictions_test.csv")
        np=read(ROOT/fresh["model_dir"]/"predictions_test.csv")
        require(len(op)==len(np),"prediction count mismatch "+str(key))
        maxabs=0.0
        for i,(a,b) in enumerate(zip(op,np)):
            require(a["sample_id"]==b["sample_id"] and
                    a["y_true"]==b["y_true"] and
                    a["y_pred"]==b["y_pred"],
                    f"first-fold sample/label/pred mismatch {key} row={i}")
            u,v=float(a["prob"]),float(b["prob"])
            require(math.isfinite(u) and math.isfinite(v),
                    "nonfinite probabilities")
            maxabs=max(maxabs,abs(u-v))
        require(maxabs<=1e-12,
                f"model output changed across pilot and 5-fold: {key}"
                f" maximum delta {maxabs}")
        require(old["realized_leaves"]==fresh["realized_leaves"] and
                old["accepted_split"]==fresh["accepted_split"],
                f"tree selection changed for {key}")
        report.append(dict(
            world=key[0],dataset=key[1],method=key[3],
            fold=key[2],n_test=len(op),max_abs_probability_diff=maxabs))
    RESULT.mkdir(parents=True)
    with (RESULT/"firstfold_prediction_parity.csv").open(
        "w",newline=""
    ) as stream:
        writer=csv.DictWriter(stream,fieldnames=list(report[0]))
        writer.writeheader()
        writer.writerows(report)
    (RESULT/"README.txt").write_text(
        "Read-only 57/57 identical first-fold predictions and tree decisions "
        "between original frozen 57-fit pilot and full five-fold run.\n"
        "This checks technical reproducibility, not independent validation.\n")
    print("FIRST_FOLD_MODEL_PARITY=",len(report),"/57",sep="")
    print("MAX_ABS_PREDICTION_DIFFERENCE=",
          max(r["max_abs_probability_diff"] for r in report),sep="")
    print("PUBLISHED_5FOLD_FIRSTFOLD_PARITY=PASS")
    print("REPORT="+str(RESULT))


if __name__=="__main__":
    main()
