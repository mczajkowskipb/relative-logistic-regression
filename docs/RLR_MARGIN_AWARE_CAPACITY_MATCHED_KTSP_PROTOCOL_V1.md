# MBRLR v1 — fair effective-capacity kTSP and margin-aligned RLR pair ranking

**Status:** original complete RLR manuscript remains **UNPUBLISHED** and locked; original RLR function `rg_parent`, `d=2`, `K=20/50`, rank margin `m=10`, `alpha=.5`, `lambda.min` and original switchBox original candidate-pool machinery are NOT altered. New source is isolated on [analysis/rlr-margin-aware-budget-match-v1](https://github.com/mczajkowskipb/relative-logistic-regression/tree/analysis/rlr-margin-aware-budget-match-v1).

## Motivation — compare real complexity, not K labels

Original paper's CompCancer11 dataset-equal weighted means derived from real already-audited 470 outer test folds:

| Original model | Nominal K | Mean nonzero coefficients (global) | Mean nonzero ternary per patient |
|---|---:|---:|---:|
| RLR10 | 10 | **9.31** | **8.41** |
| RLR20 | 20 | **16.63** | **15.10** |
| RLR50 | 50 | **28.09** | **25.55** |

These are unweighted **means of 11 per-dataset means**, not the 470 pooled test-fold means. The original 470 fold-pooled global active means are 9.491, 17.504, 29.819 (skew because 9 datasets have 50 folds and 2 have 10). It would be wrong to claim a 50-pair RLR invokes 50 active weighted relations for every patient. But equally wrong to claim it needs to **evaluate** only ~25.55 tests: the ternary rank state must be assessed to know if a fitted coefficient's term is zero. State complexity along four axes: nominal K, global nonzero beta K, mean patientwise nonzero test K, and gene/probe footprint.

**Context for kTSP50:** Literature is not monolithic: original [Tan et al. 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1987374/) cross-validated odd k up to 10; [switchBox package 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4287945/) historically defaults to 2–10; [Schliep et al. BMC Bioinformatics 2009](https://link.springer.com/article/10.1186/1471-2105-10-453) explicitly selected from k∈{10,25,50}. Thus do not **delete** historical 50 comparisons; move kTSP50 to supplementary sensitivity with other genuinely large-K baselines. PRIMARY publication comparison: TSP1, kTSP5/10/15/20 and capacity-matched kTSP K, original RLR10/20 and extended RLR. K=50 remains interpretable only after reporting model's effective capacity and original manuscript historical protocol. Avoid claiming one chosen K universally best from the already-exposed E.

## First hypothesis A: equal-budget kTSP — choose K using original RLR TRAIN only

On each outer training fold T, fit **ORIGINAL** RLR20 and RLR50 exactly. Record the globally nonzero coefficients `n_beta=Σ1(|β_j|>1e-9)` and within-T local active `a_T=(1/|T|)Σ_sΣ_{β_j !=0}1(r_j(s) !=0)`, where each `r_j(s)` is original exact ternary m10 rank relation. Choose matched budget K from **T only**, `K=max(2,min(50,round(n_beta)))` and `K=max(2,min(50,round(a_T)))` for RLR20 and RLR50 separately. Fit GENUINE switchBox kTSP on the exact same T and E with this K; record **realized** K (can be smaller than requested). Native kTSP predicted labels, not imposed 0.5 voting-threshold reinterpretation, retained. Additional real fixed K10,15,20,25 controls distinguish commonly used small panels from capacity-matched comparisons.

**Caution:** this is a deliberately *RLR-informed benchmark selection protocol*: kTSP's budget is determined by another model fit on training labels. That is a fair predictive-capacity control because it reads no E labels, but is NOT a vanilla independent kTSP default K-selection. For fully autonomous kTSP classical selection, separately add switchBox `krange=2:10` on T with its native default training procedure, avoiding family-size extrapolation. No kTSP optimum should be picked using E results.

## Second hypothesis B: margin-aware feature selection for same original RLR model

The switchBox `SWAP.CalculateSignedScore` candidate ranking is primarily based on signs of gene inequalities; original RLR actually uses three-valued **signed rank margin m10**, where `r_{ij}(s)=+1` if `rank_s(j)-rank_s(i)>=10`, `-1` if ≤−10, otherwise `0`. Reorder the **same original top1000 training switchBox candidate pairs**, no new tests or value normalization, by the frozen train-only score

`Score_ij=|mean(r_ij | class1)-mean(r_ij | class0)| × (0.25+0.75 min(P(r_ij !=0|class0), P(r_ij !=0|class1))).`

Then apply literally original `apply_disjointness(d=2,K=20 or 50)` and original binary glmnet as before. This is a simple, transparent feature-selection ablation; not a novel learned decision tree and not yet empirically proven to improve F1. The support penalty and rank gap are fixed pre-E. Compare original RLR20/50 on EXACT same T/E and record changes in nonzero coefficients, local activity, F1, Macro-F1, BA, AUC, logloss and Brier. Critically compare to kTSP with matched **T-derived** K.

All rank features and selector scores are invariant under any patient-specific **strictly increasing** transform of expression intensities, as long as the same assay feature/probe universe and identifiers are used; there are no learned numeric raw-value cutoffs or out-of-range normalizations. A dedicated R unit checks the feature ranking on monotone-transformed patients and K clipping/rounding. It does not imply cross-platform probes may be missing; source harmonization remains essential.

## Real execution / quality bar

[First-outer-fold 11-cohort GitHub Actions #40588678830](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40588678830): 11 real original-paper-named Schliep source reconstructions; 11 jobs + independent Python auditor, **12 method arms ×11 = 132 auditable E model/dataset scores**, patients identical on all comparators. Model selection and any K matching exclusively from outer T. Everything previously exposed is exploratory/development, NOT independent confirmation. **Original exact processed ARFF SHA/probe parity unverified**, and historical other paper worlds PTB fixed, PTB LOCO and Legacy8 not included.

Go/no-go:
- Do not claim margin-aware RLR beats kTSP from one good cohort; inspect both RLR20 and RLR50 across 11 independent datasets, especially rare-class Pomeroy and previously unreliable West, and compare true independent 11-cohort F1/AUC to matched-K kTSP. If first fold suggests real reproducible gains, run same FROZEN variant on original repeated 470 folds and compute per-dataset paired inferential intervals, not 470 pseudo-independent samples.
- Do not claim newly proposed tree-related novelty for a *global* margin-aware rank logistic model: it is possible **extension of author's unpublished RLR original manuscript**, offering more aligned feature screening and a stronger fair comparator. For user's habilitation-specific tree goal, only promote a genuinely learned conditional algorithm if it outperforms non-tree models on controlled validation. Prior C-RET static certified early-exit and contextual tree had small/no added advantage and must not be renamed a new contribution.
- Full original manuscript code and locked TeX unchanged.
