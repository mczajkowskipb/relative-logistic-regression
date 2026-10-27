#!/usr/bin/env python3
"""Frozen 11-dataset inference for full470 rank-margin RLR vs capacity-matched real kTSP.

Run ONLY on independently patient-replayed all11_dataset_mean_matched_budget.csv.
Inferential n=11 biological cohorts, NEVER n=470 overlapping resampling folds.
"""
import argparse,csv,itertools,random,statistics,math
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument("--dataset-means",required=True)
p.add_argument("--out",required=True)
p.add_argument("--bootstrap",type=int,default=30000)
a=p.parse_args()
with open(a.dataset_means,encoding="utf-8",newline="") as f:rows=list(csv.DictReader(f))
datasets=sorted({x["dataset"] for x in rows})
expected=("alizadeh-2000-v1","armstrong-2002-v1","bittner-2000",
 "chen-2002","golub-1999-v1","gordon-2002","pomeroy-2002-v1",
 "shipp-2002-v1","singh-2002","west-2001","yeoh-2002-v1")
arms=("RLR20","MARGIN_RLR20","RLR50","MARGIN_RLR50",
 "KTSP10","KTSP15","KTSP20","KTSP25",
 "KTSP_M20_GLOBAL","KTSP_M20_LOCAL",
 "KTSP_M50_GLOBAL","KTSP_M50_LOCAL")
if set(datasets)!=set(expected) or len(rows)!=11*len(arms) or a.bootstrap<1000:
 raise ValueError("wrong original 11 dataset ×12 method coverage or bootstrap")
D={(z["dataset"],z["method"]):z for z in rows}
if len(D)!=len(rows) or any((ds,m) not in D for ds in datasets for m in arms):
 raise ValueError("duplicate/missing dataset+method")
metrics=("positive_f1","macro_f1","ba","auc","logloss","brier",
         "global_nonzero_pairs","test_mean_active_pairs",
         "train_mean_active_pairs","test_mean_tested_pairs")
summary=[]
for method in arms:
 z={"method":method,"cohorts":11,"folds":sum(int(D[(ds,method)]["n_folds"]) for ds in datasets)}
 if z["folds"]!=470:raise ValueError("not 470 fold dataset counts")
 z.update({key:statistics.mean(float(D[(ds,method)][key]) for ds in datasets)
           for key in metrics})
 summary.append(z)
# Comparisons and endpoints fixed ahead of accessing 470 outcome data.
comparisons=[
 ("MARGIN_RLR20","RLR20","PRIMARY: margin selector effect"),
 ("MARGIN_RLR20","KTSP10","PRIMARY: common kTSP"),
 ("MARGIN_RLR20","KTSP15","PRIMARY: common kTSP"),
 ("MARGIN_RLR20","KTSP20","PRIMARY: common kTSP"),
 ("MARGIN_RLR20","KTSP_M20_GLOBAL","PRIMARY: matched global K"),
 ("MARGIN_RLR20","KTSP_M20_LOCAL","PRIMARY: matched local K"),
 ("MARGIN_RLR50","RLR50","SUPPLEMENT: large K"),
 ("RLR20","KTSP15","REFERENCE: original method"),
 ("RLR20","KTSP20","REFERENCE: original method"),
 ("RLR50","KTSP_M50_LOCAL","SUPPLEMENT: matched large K")]
res=[]
rnd=random.Random(20261027)
for aa,bb,scope in comparisons:
 for metric in ("positive_f1","macro_f1","auc","logloss"):
  dd=[float(D[(ds,aa)][metric])-float(D[(ds,bb)][metric]) for ds in datasets]
  observed=statistics.mean(dd)
  # Valid under null exchangeable paired signed dataset differences;
  # exact enumeration across 11 independent dataset units (2**11).
  extreme=0
  for signs in itertools.product((-1,1),repeat=11):
   s=sum(d*v for d,v in zip(dd,signs))/11
   if abs(s)>=abs(observed)-1e-12:extreme+=1
  pv=extreme/2048
  boots=sorted(statistics.mean(rnd.choices(dd,k=11)) for _ in range(a.bootstrap))
  ll=boots[int(.025*a.bootstrap)]
  hi=boots[min(a.bootstrap-1,int(.975*a.bootstrap))]
  wins=sum(d < -1e-12 if metric=="logloss" else d>1e-12 for d in dd)
  ties=sum(abs(d)<=1e-12 for d in dd)
  res.append({
    "method_a":aa,"method_b":bb,"scope":scope,"metric":metric,
    "delta_mean":observed,"bootstrap95_lo":ll,"bootstrap95_hi":hi,
    "signflip_two_sided_p":pv,"wins":wins,"ties":ties,"losses":11-wins-ties,
    "near_tie_f1_0005":sum(abs(d)<=.005 for d in dd) if "f1" in metric else ""
  })
# Holm across 6 prespecified PRIMARY positiveF1 comparisons.
family=[z for z in res if z["scope"].startswith("PRIMARY") and
        z["metric"]=="positive_f1"]
ordered=sorted(family,key=lambda z:z["signflip_two_sided_p"])
adj=0
for i,z in enumerate(ordered):
 adj=max(adj,(len(ordered)-i)*z["signflip_two_sided_p"])
 z["holm_primary_f1_p"]=min(1.,adj)
for z in res:z.setdefault("holm_primary_f1_p","")
out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
for name,table in [
 ("all11_model_dataset_equal_means.csv",summary),
 ("paired_11_dataset_signflip_bootstrap.csv",res)
]:
 with (out/name).open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=table[0].keys());w.writeheader();w.writerows(table)
with (out/"final_11_dataset_science_inference.md").open("w",encoding="utf-8") as f:
 f.write("# Margin-aware RLR and matched-complexity real kTSP — full 470 folds\n\n")
 f.write("Eleven original paper named publicly reconstructed source datasets; exact processed original ARFF parity not verified. All these datasets were repeatedly exposed during algorithm development. Thus this study is **exploratory**, not confirmatory generalization. One inferential unit is an entire dataset, not a fold.\n\n")
 f.write("## 11 dataset-equal-weighted descriptive performance\n\n")
 f.write("| Method | Positive F1 ↑ | Macro-F1 ↑ | AUC ↑ | LL ↓ | Nonzero beta/global K | Active ternary patient K | Tested rel. K |\n|---|---:|---:|---:|---:|---:|---:|---:|\n")
 for m in summary:
  f.write("| %s | %.4f | %.4f | %.4f | %.4f | %.2f | %.2f | %.2f |\n"%
   (m["method"],m["positive_f1"],m["macro_f1"],m["auc"],m["logloss"],
    m["global_nonzero_pairs"],m["test_mean_active_pairs"],m["test_mean_tested_pairs"]))
 f.write("\n## Paired differences across 11 independent dataset units (unadjusted descriptive except Holm-primary F1)\n\n")
 f.write("| A − B | Metric | Δ mean | 95%% bootstrap | exact p | W/T/L | Holm primary F1 p |\n|---|---|---:|---:|---:|---:|---:|\n")
 for z in res:
  if z["metric"] not in ("positive_f1","auc","logloss"):continue
  adj=z["holm_primary_f1_p"]
  f.write("| %s − %s | %s | %+.4f | [%+.4f,%+.4f] | %.4f | %d/%d/%d | %s |\n"%
   (z["method_a"],z["method_b"],z["metric"],z["delta_mean"],
    z["bootstrap95_lo"],z["bootstrap95_hi"],z["signflip_two_sided_p"],
    z["wins"],z["ties"],z["losses"],
    ("%.4f"%adj if isinstance(adj,float) else "—")))
 f.write("\n**Interpretation:** class-label F1 and probability AUC/logloss are separate scientific goals. No noninferiority/equivalence claim based on p>0.05. RLR locally active relations are not the same as pair tests necessary to determine their status. Matched kTSP K is trained on original RLR global/within-T local size without E feedback; it is not vanilla kTSP K selected independently. Actual CUDA TSPDT not tested; no tree novelty can be claimed from a changed global pair order.\n")
print("MB_RLR_FULL470_DATASET_LEVEL_SIGNFLIP_AND_BOOTSTRAP=PASS")
for x in summary:
 print("MB_RLR_DATASET_MEAN %s F1=%.6f MacroF1=%.6f AUC=%.6f LL=%.6f global=%.3f active=%.3f"%
 (x["method"],x["positive_f1"],x["macro_f1"],x["auc"],x["logloss"],
  x["global_nonzero_pairs"],x["test_mean_active_pairs"]))
