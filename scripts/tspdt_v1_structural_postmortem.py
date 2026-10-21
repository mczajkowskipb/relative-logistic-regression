#!/usr/bin/env python3
"""Read-only structural postmortem on all 95 frozen TSPDT-inspired v1 trees.

DOES NOT FIT OR MODIFY MODELS. The result is a structural diagnostic,
not another held-out performance comparison or hyperparameter search.
"""
from __future__ import annotations

import collections
import csv
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
RUN=ROOT/"results/tspdt_recon_19x5_v1"
OUT=RUN/"structural_postmortem_v1"
EXPECTED_METHOD="TSPDT2011_INSPIRED_RECON_V1"


def require(ok: bool,message: str):
    if not ok: raise SystemExit("TSPDT_STRUCTURE_STOP: "+message)


def rows(path: Path):
    require(path.is_file(),"missing "+str(path))
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    require(Path.cwd().resolve()==ROOT,"run from separate TSPDT worktree root")
    require(not OUT.exists(),"refusing to overwrite prior structural audit")
    require((RUN/"status.txt").read_text().strip()=="RUNNER_EXIT=0",
            "frozen TSPDT job has not successfully finished")
    require((RUN/"audit_v1/per_fold.csv").is_file(),
            "independent 95-fit replay audit absent")
    manifest=rows(RUN/"run_manifest.csv")
    require(len(manifest)==95 and all(r["status"]=="OK" for r in manifest),
            "require 95 successful, frozen TSPDT fits")
    out=[]
    per_dataset=collections.defaultdict(list)
    for entry in manifest:
        key=(entry["world"],entry["dataset"],entry["split_id"])
        require(entry["method"]==EXPECTED_METHOD,
                "non-reconstruction model in 95-fit manifest")
        path=(ROOT/entry["model_dir"]).resolve()
        require(path.is_relative_to(RUN.resolve()),
                "model path outside frozen result root")
        node_rows=rows(path/"nodes.csv")
        saved=json.loads((path/"metrics_test.json").read_text())
        require(saved["family"]==
                "TSPDT2011_INSPIRED_RECONSTRUCTION_V1",
                "wrong model family")
        require(len(node_rows)>=1,"missing root")
        internal=[r for r in node_rows if r["is_leaf"]=="FALSE"]
        leaves=[r for r in node_rows if r["is_leaf"]=="TRUE"]
        require(len(internal)+len(leaves)==len(node_rows),
                "unexpected is_leaf encoding")
        require(len(internal)==int(entry["internal_nodes"]) and
                len(leaves)==int(entry["leaves"]),
                "structural manifest differs from saved nodes")
        maxdepth=max(int(r["depth"]) for r in node_rows)
        require(maxdepth<=3 and
                int(node_rows[0]["depth"])==0,
                "depth exceeds preregistered reconstruction")
        used=set()
        for r in internal:
            require(bool(r["split_a"]) and bool(r["split_b"]),
                    "missing gene in nonterminal node")
            used.update((r["split_a"],r["split_b"]))
        selected=saved.get("selected_genes")
        require(isinstance(selected,list) and 2<=len(selected)<=100,
                "selected gene inventory absent or invalid")
        require(used.issubset(set(selected)),
                "tree split uses feature outside train-only selected set")
        rec=dict(
            world=key[0],dataset=key[1],split_id=key[2],
            n_train=int(entry["n_train"]),n_test=int(entry["n_test"]),
            selected_genes=len(selected),
            internal_nodes=len(internal),leaves=len(leaves),
            stump=(len(internal)==0),max_depth=maxdepth,
            distinct_split_genes=len(used),
            min_leaf_train_support=min(int(r["n"]) for r in leaves),
        )
        out.append(rec)
        per_dataset[(key[0],key[1])].append(rec)

    require(len(out)==95 and len(per_dataset)==19 and
            all(len(v)==5 for v in per_dataset.values()),
            "not 19 complete datasets / 95 folds")
    summary=[]
    for (world,dataset),d in sorted(per_dataset.items()):
        summary.append(dict(
            world=world,dataset=dataset,
            folds=5,
            root_only_folds=sum(z["stump"] for z in d),
            mean_internal_nodes=statistics.mean(
                z["internal_nodes"] for z in d),
            mean_leaves=statistics.mean(z["leaves"] for z in d),
            mean_max_depth=statistics.mean(z["max_depth"] for z in d),
            mean_distinct_split_genes=statistics.mean(
                z["distinct_split_genes"] for z in d),
            min_observed_training_leaf_size=min(
                z["min_leaf_train_support"] for z in d),
            selected_gene_count_min=min(z["selected_genes"] for z in d),
            selected_gene_count_max=max(z["selected_genes"] for z in d),
        ))
    OUT.mkdir(parents=True)
    for name,data in (("per_fold.csv",out),("per_dataset.csv",summary)):
        with (OUT/name).open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    aggregates=dict(
        method="TSPDT2011_INSPIRED_RECONSTRUCTION_V1",
        fits=95,datasets=19,
        folds_pruned_to_root=sum(r["stump"] for r in out),
        mean_internal_nodes=statistics.mean(r["internal_nodes"] for r in out),
        mean_leaves=statistics.mean(r["leaves"] for r in out),
        max_tree_depth=max(r["max_depth"] for r in out),
        allowed_max_depth=3,
        max_feature_count_per_fold=max(r["selected_genes"] for r in out),
        feature_selector="MAD training-only top100",
        paper_feature_selector="Relief-F training-only top1000 (2011)",
        pairing_search_space_100=100*99//2,
        pairing_search_space_1000=1000*999//2,
        original_2011_source_implementation=False,
        is_inferential_claim=False,
    )
    (OUT/"aggregate.json").write_text(json.dumps(
        aggregates,indent=2)+"\n",encoding="utf-8")
    (OUT/"README.txt").write_text(
        "Structural postmortem on ALL 95 frozen TSPDT-inspired v1 models.\n"
        "No fitting, no selected-dataset filtering, no threshold tuning.\n"
        "Used to distinguish root-only trees from nontrivial models and "
        "to quantify the MAD f100 restriction.\n")
    print("TSPDT_MODELS=95")
    print("ROOT_ONLY_TREES="+str(aggregates["folds_pruned_to_root"]))
    print(f"MEAN_INTERNAL_NODES={aggregates['mean_internal_nodes']:.4f}")
    print(f"MEAN_LEAVES={aggregates['mean_leaves']:.4f}")
    print("MAX_OBSERVED_DEPTH="+str(aggregates["max_tree_depth"]))
    print("MAX_SELECTED_GENES="+str(aggregates["max_feature_count_per_fold"]))
    print("TSPDT_STRUCTURAL_POSTMORTEM=PASS")
    print("REPORT="+str(OUT))


if __name__=="__main__":
    main()
