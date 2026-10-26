# RLR tree extensions — six original-paper cohort evidence gate, 2026-10-26

## Scope and scientific question

The original RLR article is **unpublished**. The research goal is a **single enhanced RLR manuscript** with a defensible interpretable TREE contribution for the habilitation achievements, not an unrelated second manuscript. Experiments are limited to original RLR paper's eleven CompCancer cohorts. Here: Chen-2002, Singh-2002, Yeoh-2002-v1, Armstrong-2002-v1, Bittner-2000, Gordon-2002. Original RLR functions, original manuscript lock and main branch unchanged. All data are public Schliep **reconstructions**, not yet byte-parity-verified original paper ARFF files; this is an exploratory prepublication pilot, not independent external confirmation.

## Results 1: genuine ARLET three-way RLR tree — 3/3 all CI SUCCESS

[ARLET workflow run #40366023886](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40366023886). 3 heldout 3-fold experiments × 6 arms = 54 completed method-fold fits. Independent patient-level audits PASS for all 3 cohorts. The root uses 3/5 rel. features, optional deferred route is a separate 7-relation logistic multitest; coefficients at both nodes are trained on all relevant outer TRAIN data, not tiny terminal leaves. Pool original signed switchBox top1000 + about 250–450 seeded pair candidates from all measured genes. Feature selection and EA were fully inside train/validation of each outer fold. EA budget 3 generations ×8 individuals.

| Original paper cohort | RLR BA | kTSP9 BA | TSP1 BA | Greedy TSPDT **proxy** BA | ARLET deterministic BA | ARLET EA BA | RLR LL | ARLET deterministic LL | ARLET EA LL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Singh-2002 | 0.91231 | 0.92211 | 0.82394 | 0.84300 | **0.93137** | 0.92150 | 0.33194 | **0.24149** | 0.27664 |
| Chen-2002 | **0.91291** | 0.87331 | 0.81193 | 0.80801 | 0.90067 | 0.88924 | **0.19818** | 0.27131 | 0.30771 |
| Yeoh-2002-v1 | **0.92446** | 0.86829 | 0.82791 | 0.76660 | 0.74524 | 0.72222 | **0.10106** | 0.25352 | 0.24178 |

Full Macro-F1 and Brier metrics are in GitHub-uploaded `method_metrics.csv` artifacts. **TSPDT proxy is NOT the user's original CUDA TSPDT**; do not advertise real TSPDT victory on the basis of the proxy. Baseline `TSP_1` returns hard 0/1 probabilities, which are not calibrated and yield distorted logloss (do not emphasize them in logloss plots). The ARLET result is positive on Singh but loses vs original RLR on both Chen and Yeoh; EA is worse than deterministic on all three. **ARLET V1 FAILS the extension evidence gate**. Do not tune these exposed E test results and call them confirmatory.

## Results 2: probabilistically certified early-exit RLR tree — 3/3 SUCCESS

[Certified tree workflow #40367860184](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40367860184). Uses the literal unchanged original RLR K50 fitted model. For any prefix S of rank-margin-active `r_j∈{-1,0,+1}`, partial `L=β0+Σ_{S}β_j r_j`; worst-case remaining contribution `B=Σ_{¬S}|β_j|`; original probability is **provably** in `[σ(L-B),σ(L+B)]`. Early exit permitted only when the entire bound lies on one side of 0.5 and has width ≤0.05. Therefore original RLR classification is EXACT; midpoint approximate probability absolute error ≤0.025; its logloss/Brier are **not** guaranteed to improve.

| Cohort | Original selected NONZERO coefficients / assessed pairs | Fixed-order certificate mean pairs assessed | Mean active relations actually contributing | Early exit |
|---|---:|---:|---:|---:|
| Chen-2002 | 20.000 | 12.040 | 9.722 | 99.44% |
| Singh-2002 | 20.000 | 13.157 | 11.394 | 100% |
| Yeoh-2002-v1 | 8.667 | 4.971 | 3.211 | 87.86% |

All 529 distinct patients (3-fold coverage per dataset) had exact original RLR class; independent patient-level bound audits PASS (3×529=1587 patient-method scores). Simple absolute-coefficient order outperformed or tied a small (4×8) EA ordering budget in number of evaluated relations. Note that rank relation evaluation cost reduction does **not** imply laboratory assay cost reduction, because all original gene measurements may still be collected and within-patient ranks computed.

## Results 3: state-dependent conditional next-relation policy — 3 more original-paper cohorts PASS

[True conditional tree workflow #40369948776](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40369948776). Tests three different original-manuscript datasets. At each inconclusive node, select the next gene comparison by minimizing expected future probability-interval width, with TRAIN-only empirical ternary-conditional frequencies shrunk to unconditional frequencies; different earlier `−1,0,+1` outcomes can lead to different next comparisons. **Certificate remains valid independent of how order is chosen.** Models tested with exact original RLR logistic regression.

| Dataset | Original RLR mean nonzero pairs | Fixed order mean evaluated | Conditional tree mean evaluated | BA (all identical to original RLR) | Original RLR logloss | Conditional approximate logloss |
|---|---:|---:|---:|---:|---:|---:|
| Armstrong-2002-v1 | 29.000 | 17.167 | 17.139 | 1.00000 | .01307 | .03053 |
| Bittner-2000 | 39.333 | 30.869 | 30.813 | .81349 | .43535 | .43631 |
| Gordon-2002 | 35.333 | 14.984 | 14.984 | .98152 | .03275 | .05151 |

Independent E-patient certificate audits PASS for all 3, exact class parity + ≤0.025 posterior approximation, 3 disjoint outer E folds per dataset. Earlier failed first CI attempt was caused by missing reporting field (`NA` early-stop fraction), not a mathematical certificate failure. Fixed code rerun success at #40369948776. Conditional feature acquisition was **not meaningfully better** than a fixed descending-|β| order. The tree is a genuine ternary-outcome-dependent test selection policy, but its predictive advantage is mathematically not expected because both represent the original RLR; it must justify learning complexity through acquisition cost savings, which this pilot did not demonstrate.

## Recommendation: what IS and is NOT supported

- **Positive real result:** exact original RLR class on 6 binary original-paper datasets with large reduction in evaluated relation comparisons, and a provable per-patient posterior interval. Original RLR's local inactive relation concept is exploited (active means `r=±1`, inactive `r=0`; do not conflate inactive with statistical irrelevance of β=0).
- **Negative/mixed real result:** the new independently trained nonlinear ARLET root/rescue tree is strongly dataset-dependent, while EA search never beat deterministic in the tested cohorts. No robust claim of higher BA/Marco-F1 vs original RLR; original TSPDT CUDA NOT compared.
- **Scientific novelty is NOT established:** existing sequential feature acquisition, anytime linear-model evaluation and tree representations precede this work. A fixed-order proof-based early exit is not obviously a new learned decision tree. State-dependent next-test policy is closer to a genuine tree but added almost no cost benefit.
- **No manuscript edits / no publication-ready method:** next step should establish whether a hierarchical *relational* policy offers a meaningful unique contribution beyond classical early-exit linear thresholds, and whether the full paper original 11 datasets and 470 splits support it with exact source verification. If not, the existing unpublished RLR should be assessed separately as a complementary relational method rather than adding weak tree novelty just for habilitation.
