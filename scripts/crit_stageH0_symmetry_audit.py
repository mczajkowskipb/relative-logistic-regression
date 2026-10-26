#!/usr/bin/env python3
"""Stage H0: read-only, versioned, orientation-invariant planted XOR audit.

Inputs: frozen Stage-G archives. Creates *one* new H0 result directory;
never touches Stage-G. H0 is retrospective scientific annotation, not
an additional independent benchmark or an algorithm update.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict, Counter
from pathlib import Path
from statistics import mean

FROZEN_SHA="6ca44c30a439dbdb97b3a45f7e8f45432da3bd59"
EXPECTED_EXIT="STAGEG_BATCH_EXIT=0"
Z="gene01>gene02"
R="gene03>gene04"
SCENARIOS=("null","global_pair","gate_only","xor_interaction")
GROUPS=("REUSED_DEV","FRESH_DEV")
METHODS=("ORIGINAL_RLR","GLOBAL_CAL","GATE_ONLY",
         "ADDITIVE_PAIRS","GLOBAL_INTERACTION","CONDITIONAL_TREE")

def must(ok: bool,msg: str)->None:
    if not ok:
        raise SystemExit("H0_STOP: "+msg)

def read(p: Path):
    must(p.is_file(),"missing archive "+str(p))
    with p.open("r",encoding="utf-8",newline="") as f:
        a=list(csv.DictReader(f))
    must(len(a)>0,"empty archive "+str(p))
    return a

def val(s):
    n=float(s)
    must(math.isfinite(n),"nonfinite archived metric")
    return n

def classify(c):
    gate=c["gate"]
    interaction=c["interaction_pair"]
    left=c["false_pair"]
    right=c["true_pair"]
    direct=gate==Z and left==R and right==R
    swapped=gate==R and left==Z and right==Z
    # Check both leaf realization AND interaction provenance.
    if direct:
        must(interaction==R,
             "direct tree uses wrong interaction provenance")
        return "DIRECT_TREE"
    if swapped:
        must(interaction==Z,
             "swapped tree uses wrong interaction provenance")
        return "SWAPPED_TREE"
    return "OTHER"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--stageg",type=Path,required=True,
                    help="frozen Stage-G results directory")
    ap.add_argument("--out",type=Path,required=True,
                    help="new H0 results directory, must not exist")
    args=ap.parse_args()
    old=args.stageg.resolve(strict=True)
    out=args.out
    must(not out.exists() and not out.is_symlink(),
         "output exists or symlink")
    must((old/"status.txt").read_text().strip()==EXPECTED_EXIT,
         "Stage-G batch did not finish cleanly")
    must((old/"SOURCE_SHA.txt").read_text().strip()==FROZEN_SHA,
         "Stage-G source SHA mismatch")
    a=(old/"independent_audit.log").read_text(encoding="utf-8")
    for gate in (
      "STAGEG_34560_PATIENT_METRICS_INDEPENDENTLY_AUDITED=PASS",
      "STAGEG_48_STUDY_NOMINATION_AND_FALLBACK_AUDITED=PASS",
      "STAGEG_24_REUSED_STUDIES_EXACT_ORIGINAL_RLR_PARITY=PASS"):
        must(gate in a,"frozen Stage-G independent metric audit missing "+gate)

    candidates=read(old/"candidate_provenance.csv")
    study=read(old/"study_summary.csv")
    must(len(study)==288,"unexpected Stage-G study inventory")
    expected={(s,n,m) for s in SCENARIOS
              for n in range(1,13) for m in METHODS}
    results={(x["scenario"],int(x["replicate"]),x["method"]):x
              for x in study}
    must(set(results)==expected and len(results)==288,
         "Study method IDs duplicate/missing")
    grouped=defaultdict(list)
    for x in candidates:
        key=(x["scenario"],int(x["replicate"]))
        must(x["scenario"] in SCENARIOS and int(x["replicate"]) in range(1,13),
             "unexpected candidate group")
        grouped[key].append(x)
    must(len(grouped)==48,"missing candidate studies")

    xor_rows=[]
    table=[]
    for s in SCENARIOS:
      for n in range(1,13):
        key=(s,n)
        cs=sorted(grouped[key],key=lambda x:int(x["D_rank"]))
        must(1<=len(cs)<=6 and
             [int(x["D_rank"]) for x in cs]==list(range(1,len(cs)+1)),
             "candidate D ranks missing/duplicated")
        chosen=[x for x in cs if int(x["chosen_by_D"])==1]
        must(len(chosen)==1,"multiple D-selected candidates")
        c=chosen[0]
        must(int(c["D_rank"]) ==
             1+max(range(len(cs)),key=lambda i:val(cs[i]["D_cv_lower"])),
             "chosen candidate is not maximum within-D crossfit gain")
        group="REUSED_DEV" if n<=6 else "FRESH_DEV"
        base=results[s,n,"CONDITIONAL_TREE"]
        glob=results[s,n,"GLOBAL_INTERACTION"]
        must(base["group"]==group and glob["group"]==group and
             int(base["shortlisted_n"])==len(cs) and
             int(base["D_selected_rank"])==int(c["D_rank"]),
             "study/candidate provenance mismatch")
        if s=="xor_interaction":
            cls=classify(c)
            rec=dict(group=group,replicate=n,class_=cls,
                     D_rank=int(c["D_rank"]),gate=c["gate"],
                     false_pair=c["false_pair"],
                     true_pair=c["true_pair"],
                     interaction_pair=c["interaction_pair"],
                     tree_accepted=int(base["accepted"]),
                     tree_deltaLL=val(base["delta_vs_original"]),
                     global_accepted=int(glob["accepted"]),
                     global_deltaLL=val(glob["delta_vs_original"]))
            xor_rows.append(rec)
            must(int(base["shortlisted_planted_joint"])==1,
                 "planted pair absent from Stage-G top-six on XOR")
        if s=="null" and results[s,n,"CONDITIONAL_TREE"]["accepted"]!="0":
            raise SystemExit("H0_STOP: null tree accepted; update evidence")
    must(len(xor_rows)==12,"expected 12 XOR development studies")
    counts=Counter(r["class_"] for r in xor_rows)
    must(counts["DIRECT_TREE"]==5 and counts["SWAPPED_TREE"]==7 and
         counts["OTHER"]==0,"observed 5 direct + 7 swapped not reproduced")
    must(all(x["tree_accepted"]==1 and x["tree_deltaLL"]<0
             and x["global_accepted"]==1 and x["global_deltaLL"]<0
             for x in xor_rows),
         "model acceptance/heldout gain mismatch")
    for group in GROUPS:
        xs=[x for x in xor_rows if x["group"]==group]
        must(len(xs)==6,"XOR split not 6+6")
        table.append({
            "group":group,"studies":len(xs),
            "direct":sum(x["class_"]=="DIRECT_TREE" for x in xs),
            "swapped":sum(x["class_"]=="SWAPPED_TREE" for x in xs),
            "symmetric_exact_recovery":sum(x["class_"]!="OTHER" for x in xs),
            "tree_accepted":sum(x["tree_accepted"] for x in xs),
            "tree_E_wins":sum(x["tree_deltaLL"]<0 for x in xs),
            "tree_mean_deltaLL":mean(x["tree_deltaLL"] for x in xs),
            "global_accepted":sum(x["global_accepted"] for x in xs),
            "global_mean_deltaLL":mean(x["global_deltaLL"] for x in xs)
        })
    doc={
      "scientific_stage":"REUSED_AND_FRESH_SYNTHETIC_DEVELOPMENT_ONLY",
      "stageg_frozen_sha":FROZEN_SHA,
      "evidence":"Stage G exact candidate provenance and audited patient metrics",
      "planted_pair1":Z,"planted_pair2":R,
      "orientation_invariant_explanation":
        "z*r=r*z; design spaces overlap but regularized fitted predictions need not be identical",
      "classification_counts":dict(counts),
      "orientation_invariant_exact_recovery":12,
      "xor_studies":12,
      "null_conditional_tree_acceptance":"0/12",
      "families":table,
      "individual_studies":xor_rows,
      "limitations":[
        "Not new independent performance evidence: reused and inspected synthetic development experiments.",
        "No type-I error calibration; 0/12 null acceptances is insufficient.",
        "Full C(32,2)=496 candidate pool does not scale to omics.",
        "Primary corrected discovery metric allows gate/leaf pair role swaps."
      ]
    }
    out.mkdir(parents=False,exist_ok=False)
    (out/"symmetry_corrected.json").write_text(
        json.dumps(doc,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    with (out/"xor_selected_12.csv").open("w",encoding="utf-8",newline="") as f:
        cols=("group","replicate","class_","D_rank","gate","false_pair",
              "true_pair","interaction_pair","tree_accepted","tree_deltaLL",
              "global_accepted","global_deltaLL")
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader()
        w.writerows(xor_rows)
    for v in table:
        print("H0_FAMILY "+v["group"]+
              f" direct={v['direct']} swapped={v['swapped']} "+
              f"recovered={v['symmetric_exact_recovery']}/6 "+
              f"tree_accept={v['tree_accepted']}/6 "+
              f"tree_mean_deltaLL={v['tree_mean_deltaLL']:+.8f} "+
              f"global_mean_deltaLL={v['global_mean_deltaLL']:+.8f}")
    print("H0_SYMMETRY_EXACT_12_OF_12=PASS")
    print("H0_STAGEG_FROZEN_34560_PATIENT_AUDIT_VERIFIED=PASS")
    print("H0_REUSED_SYNTHETIC_DEVELOPMENT_ONLY")

if __name__=="__main__":
    main()
