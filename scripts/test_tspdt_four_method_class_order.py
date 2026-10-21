#!/usr/bin/env python3
"""Synthetic regression: original R/ARFF class order may invert lexical order.

No actual cohort files; no model fits, no result directory writes.
"""
from __future__ import annotations
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/"scripts"))
from tspdt_four_method_published_report import (
    infer_reference_class_order, score
)


def main() -> None:
    key=("compcancer_2x2","synthetic_reversed_levels","rep00_fold00")
    original=[
        {"sample_id":"s1","y_true":"0","y_pred":"0","prob":"0.9"},
        {"sample_id":"s2","y_true":"1","y_pred":"1","prob":"0.1"},
        {"sample_id":"s3","y_true":"0","y_pred":"0","prob":"0.8"},
        {"sample_id":"s4","y_true":"1","y_pred":"1","prob":"0.2"},
    ]
    ids=[r["sample_id"] for r in original]
    true=[r["y_true"] for r in original]
    levels=infer_reference_class_order(original,true,key)
    assert levels==("1","0"),levels
    actual=score(original,ids,true,levels,"synthetic-reversed")
    expected=-(math.log(0.9)+math.log(0.9)+
               math.log(0.8)+math.log(0.8))/4
    assert abs(actual["logloss"]-expected)<1e-14,actual
    assert abs(actual["brier"]-0.025)<1e-14,actual
    assert actual["accuracy"]==1.0
    assert actual["balanced_accuracy"]==1.0
    assert actual["macro_f1"]==1.0
    # A numeric label can be lexically sorted incorrectly while
    # remaining valid in the original categorical R factor order.
    wrong=score(original,ids,true,("0","1"),"incorrect-lexical")
    assert wrong["logloss"]>1.5
    assert wrong["logloss"]>actual["logloss"]+1.0

    # The exact R baseline threshold >=0.5 determines positive class.
    equal=[
        {"sample_id":"a","y_true":"1","y_pred":"0","prob":"0.5"},
        {"sample_id":"b","y_true":"0","y_pred":"1","prob":"0.1"},
    ]
    assert infer_reference_class_order(
        equal,[r["y_true"] for r in equal],key)==("1","0")

    # Label/probability contradictions must fail closed.
    bad=[dict(row) for row in original]
    bad[0]["y_pred"]="1"
    try:
        infer_reference_class_order(bad,true,key)
        raise AssertionError("inconsistent RLR predictions unexpectedly accepted")
    except SystemExit as e:
        assert str(e).startswith("FOUR_METHOD_REPORT_STOP:"),str(e)

    print("FOUR_METHOD_CLASS_ORDER_REGRESSION=PASS")
    print("REVERSED_ARFF_CLASS_ORDER=1,0")
    print(f"REFERENCE_LOGLOSS={expected:.12f}")


if __name__=="__main__":
    main()
