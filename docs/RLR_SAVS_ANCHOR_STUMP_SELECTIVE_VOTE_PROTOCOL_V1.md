# SAVS V1 — selective ordinal anchor-stump and one-sided abstaining relational votes

**Status:** Implemented in new isolated experimental branch. Actual R test CI and 11 original CompCancer cohort first outer T/D-V diagnostic launched in [#40552370636](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40552370636). Original **UNPUBLISHED** RLR implementation and locked article remain unchanged.

## Distinguish claims that would be easy to confuse

1. **Original RLR's ternary relation** is `r_ij(x)=+1` when `rank(j)-rank(i)>=10`; `-1` when `rank(j)-rank(i)<=-10`; `0` otherwise. Both +1 and -1 can contribute to logistic prediction if the global model fitted a nonzero coefficient; *only 0* makes the term inactive for that patient (and nonzero coefficient selection by glmnet is a second separate sparsity mechanism). It is **not correct** to assert simply `x_i>x_j` false makes the original RLR pair inactive: the reverse ordering can yield a nonzero negative feature.
2. **New SAVS** uses truly *one-sided* rule `h_j(x)=I[s_j(rank(j)-rank(i))>=m]`. A prelearned D-only directional class vote `v_j in {-1,+1}` enters the patient's branch vote `S_b(x)=sum_{j in branch b} v_j h_j(x)`. When the condition does not hold, this rule **abstains**. Class0 votes come from separately selected class0-associated ONE-SIDED inequalities, not treating a failed class1 rule as a class0 vote.
3. **Prediction structure:** root gate counts `target>anchor` comparisons for 3 or 5 D-nominated anchor genes. For each of the two branches, at most **two positive and two negative** gene-disjoint one-sided rules chosen from original D-only ranked switchBox gene-pair pool. No extra per-branch glmnet fitting or minimum minority-class count; class prediction uses sparse vote sign `-1/0/+1` and fixed empirical Bayes shrinkage of positive-class frequencies using D training counts. Absent active rules receives smoothed branch prior; probabilities remain finite even in one-class leaf. Root anchor vote is strictly within patient.
4. **Interpretation:** each prediction can show: root anchor count, branch, list of rules that fired and their class directions, count of 0–4 active local rules, max root 3/5 plus max 4 local comparisons, patient-specific vote margin, resulting posterior. This is a *conditional selected evidence rule*, not a claim that just having fewer total listed pairs makes accuracy better.
5. **Strict per-patient monotone invariance**: for fixed measured assay probe universe and any strictly increasing transformation of the same patient's feature values, ranks and all branch and local decisions are identical. Neither ratio nor value threshold nor cohort-wide scaling required. Rank-margin mode `RLR10` reuses original m10, `RELATIVE10PCT` uses max(10,ceil(.1*(p-1))) where p is number of probes; both fixed before V; these are distinct exploratory variants, not E-selected hyperparams.

## Scientific comparator/gates

Same D/V original first fold of ALL 11 original paper-named primary CompCancer datasets. Two fixed margin worlds × four methods = 8 arms per cohort:
- **CONDITIONAL_ANCHOR_ABSTAIN**: actual sparse model-stump with branch-specific rule panels.
- **GLOBAL_SAME_RULES**: identical exact union of nominated one-sided rule tests (same gene IDs and signs), no root routing, smoothed bucket posterior from full D.
- **GLOBAL_FOUR_RULES**: four best global sparse one-sided vote tests (2 per class), same selection/categorization, no routing.
- **ANCHOR_ROOT_ONLY**: root gate with smoothed class posterior, no child rule votes.

Logloss, positive-class F1, Macro-F1, BA, AUC, Brier plus median/mean patient-activated one-sided rules, available rules, percent all abstaining and percent with at most 2 active. Actual root comparisons **must be included** in test-time computation; do not claim a four-rule stump only evaluates four comparisons when it also evaluates 3 or 5 root comparisons.

**Important limit:** fewer evaluated within-patient *comparisons* does **not** imply fewer genes need to be measured on the assay/platform (the whole transcriptomic profile is already collected). Thus claim is local **explanation size and comparison computations**, not reduced wet-lab diagnostic costs. Sparse voting as such is also prior art in kTSP; adaptive test-time evaluation trees / cost-sensitive feature acquisition have long precedent (e.g. Xu et al., ICML 2013, https://proceedings.mlr.press/v28/xu13.html). Methodological novelty requires **the particular stable rank anchor/reference selection + branch-conditional directed one-sided rule abstention + measurable patient-specific sparsity**, shown experimentally beyond globally identical feature union.

## Direct manuscript structural evidence: not a new RLR contribution, but a baseline to beat

The user's completed **UNPUBLISHED original RLR manuscript** already reports exactly four stages of sparseness: candidate shortlist, nominal K, nonzero global coefficients, and nonzero pair states per specific patient. In its CompCancer structural table (descriptive group-level means):

| Model | Nominal K | Global active mean | Patient local active mean |
|---|---:|---:|---:|
| RLR-10 | 10 | 9.31 | 8.41 |
| RLR-20 | 20 | 16.63 | 15.10 |
| **RLR-50** | **50** | **28.09** | **25.55** |

Original RLR is **not** a 50-active-pair prediction for every patient. Its ternary pair feature is exactly zero only when the absolute signed rank gap is below m (while both far-apart signs remain active), and the fitted logistic coefficient can also be zero. This original claim is already in `manuscript_lock/RLR3_submission_final.tex`, not a newly invented result and not to be presented as specific novelty of SAVS.

The new SAVS V1 has at most four child *voting* rules per patient and 3 or 5 unconditional anchor-root comparisons. Both counts MUST be reported separately. A negative or inactive single-sided inequality does not represent the original RLR's ternary negative relation; it is an explicit no-vote. Do NOT compare only 2 activated votes versus 25.55 original active RLR relations without giving model F1/AUC and true tested comparison count and class coverage.

## Interim actual first completed D/V: Chen

[#40552370636](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40552370636), D=108/V=35:
- Mode RLR10, **actual tree**: F1 .81250, Macro-F1 .82730, AUC .95167, LL .36206, **mean active leaf rules 2.09 per V patient**, 5.7% all leaf rules abstain.
- Global SAME nominated rules: F1 .90323, Macro-F1 .91315, LL .26292, mean 3.60 active rules.
- Four-rule non-tree sparse vote: F1 **.93750**, Macro-F1 .94243, LL **.20551**, mean **1.74 active** rules.
- Root only: F1 .81250, LL .46563. Early data show sparse selective evidence is feasible, but CURRENT conditional TREE routing worsens predictive F1 substantially. The apparent performance winner so far is a **NON-TREE** model, and the exact report MUST say so.
- The alternate predeclared 10%-assay-rank margin gave identical Chen outcomes (possible m10 and m10% margins predict same rule statuses for these selected genes on these patients). No manual λ, root or m adaptation to seen V outcomes.

## GO/NO-GO

- Complete all 11 D/V sources and **independent patient-level replay** of 88 method-arm metrics, including abstentions and true comparison counts.
- Do NOT rerun 470 external E folds solely to showcase a known-negative tree; require training-only evidence across cohorts. Use an explicit D/V criterion: positive F1 and BA vs both GLOBAL_SAME_RULES and GLOBAL_FOUR_RULES, not only original RLR. New candidate root lookahead may be tested in a separate independently versioned protocol (model selection within T only), with outer E data strictly withheld during changes.
- If branch selection remains inferior to equal-rule global, retain selective vote as an **ablation** or a separately useful global 4-rule classifier, not label it a tree contribution.
- The original RLR paper's main already-audited 11×470 result remains original RLR50 F1 .868613, AUC .957076, LL .226026 vs real kTSP50 F1 .883625, AUC .962416, LL .377593; prior CFO-rank local tree F1 .868825, LL .223235. Those external E outcomes were already seen by the research team, so not an untouched independent confirmatory test for SAVS.
