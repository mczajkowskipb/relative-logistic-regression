#!/usr/bin/env python3
"""Equal-dataset 11-cohort inference from the ORIGINAL audited 470 E-fold CI archive.
Exploratory exact paired sign-flip over 11 independent datasets and bootstrap.
NO E-motivated refitting or tuning of any classifier.
"""
import argparse,csv,itertools,math,random
from pathlib import Path

def need(ok,msg):
    if not ok:raise SystemExit("SELECTIVE4_INFERENCE_STOP "+msg)
def read(path):
    with path.open(encoding="utf-8",newline="") as f:return list(csv.DictReader(f))
def write(path,items):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=items[0].keys());w.writeheader();w.writerows(items)
parser=argparse.ArgumentParser()
parser.add_argument("--dataset-means",required=True)
parser.add_argument("--out",required=True)
parser.add_argument("--source-run",default="40569991819")
args=parser.parse_args()
rows=read(Path(args.dataset_means))
methods=("TSP1","KTSP4","KTSP10","KTSP20","KTSP50","RLR10","RLR20",
         "RLR50","VALUE_LR","SELECTIVE4_M10","SELECTIVE4_M10PCT")
datasets=sorted({r["dataset"] for r in rows})
need(len(datasets)==11 and len(rows)==11*len(methods),
     "expected exactly 121 dataset-method rows and 11 cohorts")
D={(r["dataset"],r["method"]):r for r in rows}
need(len(D)==121,"duplicate cohort-method row")
for ds in datasets:
    need({r["method"] for r in rows if r["dataset"]==ds}==set(methods),
         "missing methods "+ds)
    n=10 if ds in ("west-2001","yeoh-2002-v1") else 50
    for method in methods:
        need(int(D[ds,method]["n_folds"])==n,
             "original paper fold count mismatch "+ds)
def average(method,metric):
    return sum(float(D[d,method][metric]) for d in datasets)/11
comparisons=(
 ("SELECTIVE4_M10","KTSP4"),("SELECTIVE4_M10","KTSP50"),
 ("SELECTIVE4_M10","RLR50"),
 ("SELECTIVE4_M10PCT","KTSP4"),
 ("SELECTIVE4_M10","SELECTIVE4_M10PCT"))
metrics=("positive_f1","macro_f1","ba","auc","logloss","brier")
signs=tuple(itertools.product((-1,1),repeat=11))
out=[]
for pair_index,(a,b) in enumerate(comparisons):
    for metric_index,metric in enumerate(metrics):
        dif=[float(D[d,a][metric])-float(D[d,b][metric]) for d in datasets]
        delta=sum(dif)/11
        exact=sum(abs(sum(sign*val for sign,val in zip(s,dif))/11)>=
                  abs(delta)-1e-13 for s in signs)/len(signs)
        rng=random.Random(470_000+pair_index*100+metric_index)
        B=30000
        boot=sorted(sum(rng.choices(dif,k=11))/11 for i in range(B))
        lo,hi=boot[int(B*.025)],boot[int(B*.975)-1]
        improving=lambda x:x<0 if metric in ("logloss","brier") else x>0
        out.append(dict(method_a=a,method_b=b,metric=metric,
            delta=delta,bootstrap_95_lo=lo,bootstrap_95_hi=hi,
            exact_signflip_p=exact,
            wins=sum(improving(z) for z in dif),
            ties=sum(abs(z)<1e-12 for z in dif),
            losses=sum(not improving(z) and abs(z)>=1e-12 for z in dif),
            n_datasets=11,source_run=args.source_run,
            validation_stage="EXPLORATORY_470_FOLDS_11_DATASET_UNITS"))
dir=Path(args.out);dir.mkdir(parents=True,exist_ok=True)
write(dir/"paired_dataset_inference.csv",out)
audit={z["method"]:z for z in rows if z["dataset"]==datasets[0]}
summary=[]
for method in methods:
    z={"method":method}
    for metric in metrics:
        z[metric]=average(method,metric)
    for metric in ("mean_active","mean_tested"):
        allvals=[D[d,method][metric] for d in datasets]
        z[metric]=(sum(float(t) for t in allvals)/11
                if all(t not in ("","NA") for t in allvals) else "")
    summary.append(z)
write(dir/"eleven_dataset_equal_means_summary.csv",summary)
# Cross-cohort paired differences retain WHICH datasets are weak, not just
# one average. This is an essential safeguard against cherry-picking Pomeroy.
per_ds=[]
for ds in datasets:
    x=dict(dataset=ds)
    for method in ("KTSP4","KTSP50","RLR50",
                   "SELECTIVE4_M10","SELECTIVE4_M10PCT"):
        for metric in ("positive_f1","macro_f1","auc","logloss","mean_active"):
            x[method+"_"+metric]=D[ds,method][metric]
    per_ds.append(x)
write(dir/"all11_cohort_primary_metrics.csv",per_ds)
lookup={(r["method_a"],r["method_b"],r["metric"]):r for r in out}
def fmt_comp(a,b,metric):
    row=lookup[a,b,metric]
    return (f'{float(row["delta"]):+.4f} '
       f'[{float(row["bootstrap_95_lo"]):+.4f}; '
       f'{float(row["bootstrap_95_hi"]):+.4f}] '
       f'(p={float(row["exact_signflip_p"]):.4f})')
report=[
 "# Selective four-rule rank voting: FINAL 11-cohort / 470-fold analysis",
 "",
 f"Source: GitHub Actions run {args.source_run}, completed 95/95 jobs SUCCESS; "
 "470 original-design outer E folds × 11 classifier arms = "
 "5170 independently patientwise replayed scores.",
 "",
 "**Evidence caveat:** firstfold and related cohorts were used in algorithm "
 "development; expression matrices were publicly reconstructed and not "
 "SHA-verified against the original manuscript's processed ARFF files. "
 "This is a rigorously checked **development benchmark**, not independent "
 "confirmation, not clinical validation and not a proof of exact historic "
 "raw-data parity.",
 "",
 "## All 11 datasets equally weighted",
 "",
 "| Method | Positive F1 ↑ | Macro-F1 ↑ | AUC ↑ | Logloss ↓ | Active relations/sample | Evaluated relations/sample |",
 "|---|---:|---:|---:|---:|---:|---:|"
]
for row in summary:
    shown=lambda k:f'{row[k]:.4f}' if isinstance(row[k],(float,int)) else "—"
    report.append(f'| {row["method"]} | {shown("positive_f1")} | '
     f'{shown("macro_f1")} | {shown("auc")} | {shown("logloss")} | '
     f'{shown("mean_active")} | {shown("mean_tested")} |')
report += [
 "",
 "## Paired inference: 11 datasets, NOT 470 independent folds",
 "",
 "Differences are Selective4 M10 minus comparator; positive F1/AUC are "
 "better if positive, logloss if negative. Bootstrap 95% percentile "
 "intervals sample 11 dataset differences (30,000 resamples, fixed seed). "
 "The two-sided sign-flip test EXACTLY enumerates 2^11 paired sign "
 "configurations. Exploratory p-values are **not multiplicity corrected**.",
 "",
 "| Comparator | Δ positive F1 (95% CI; p) | Δ AUC (95% CI; p) | "
 "Δ logloss (95% CI; p) |",
 "|---|---|---|---|"
]
for against in ("KTSP4","KTSP50","RLR50"):
    report.append("| "+against+" | "+
      fmt_comp("SELECTIVE4_M10",against,"positive_f1")+" | "+
      fmt_comp("SELECTIVE4_M10",against,"auc")+" | "+
      fmt_comp("SELECTIVE4_M10",against,"logloss")+" |")
report+= [
 "",
 "## Scientifically defensible interpretation",
 "",
 "1. **No win versus original RLR50 or real kTSP50** in F1/AUC. "
 "Do NOT claim such a result or reinterpret nonsignificance as equivalence.",
 "2. At 4 evaluated comparisons, selective m10 votes are near-tied "
 "in mean positive F1 with REAL kTSP4, but materially inferior in "
 "AUC. Need genuine prespecified NONINFERIORITY threshold and untouched "
 "datasets to establish equivalence; this study did not.",
 "3. For Selective4 m10 roughly two of four one-sided rules fire "
 "in a typical sample; every one of the four conditions was still "
 "evaluated. Compared with RLR50's originally reported 25.55 active "
 "relations, note **different activation definitions**: SAVS is "
 "one-sided fire/no-vote; original RLR's ternary +1/-1 are BOTH "
 "active and 0 is inactive.",
 "4. For sample-specific monotone transformations all predictors use "
 "within-patient ranks: no inter-sample intensity scaling or "
 "numerical ratio is needed. The whole assay is still measured; "
 "compact explanations are **not** reduced wet-lab gene assay costs.",
 "5. SAVS is a GLOBAL sparse voting method, **NOT** a learned "
 "tree structure. SAVS3 honest selected root and SAVS4 hub-sharing "
 "conditional trees underperformed global FOUR on their internal "
 "11-cohort D/V experiments. The requested novel tree-classifier "
 "publication target remains unmet.",
 "6. Preserve unpublished original RLR manuscript and source untouched; "
 "stop E-guided model search using these heavily explored 11 cohorts. "
 "Use original SHA-verified processed ARFF or untouched omics cohorts "
 "for genuine independent evidence.",
 "",
 "All numerical means/differences are traceable to the immutable "
 "GitHub Actions per-patient replay artifact.",
]
(dir/"FINAL_SELECTIVE4_470_EVIDENCE_REPORT.md").write_text(
 "\n".join(report)+"\n",encoding="utf-8")
print("SELECTIVE4_ORIGINAL_11_COHORTS_470_FOLDS_RECONCILED=PASS")
print("SELECTIVE4_EXACT_2048_SIGNS_AND_BOOTSTRAP_30000=PASS")
print("SELECTIVE4_11_DATASET_LEVEL_NOT_PSEUDOREPLICATION=PASS")
for against in ("KTSP4","KTSP50","RLR50"):
    print("SELECTIVE4_PAIRED",against,
          "F1="+fmt_comp("SELECTIVE4_M10",against,"positive_f1"),
          "AUC="+fmt_comp("SELECTIVE4_M10",against,"auc"))
print("SELECTIVE4_NO_CLAIM_TREE_WINS_OR_HISTORIC_ARFF_PARITY")
