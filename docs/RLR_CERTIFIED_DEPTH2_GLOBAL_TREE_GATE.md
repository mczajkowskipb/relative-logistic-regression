# RLR C-RET depth-two lookahead: bounded global structure optimization, 2026-10-26

## Research motivation and strict scope

Target: **extend the existing unpublished RLR manuscript** with a genuine interpretable tree contribution for the omics habilitation, not an unrelated second article. Current original paper includes 11 **binary** CompCancer cohorts. Current depth-two pilot is confined to `alizadeh-2000-v1` (42) and `pomeroy-2002-v1` (34), the two original-paper cohorts not yet used for prior tree-extension pilots. Reconstructed public Schliep source data have **not** been proven byte-equivalent to the original processed ARFF files.

## Mechanism

- Fit **literal frozen original RLR K50** on each outer TRAIN fold using original switchBox top1000 pair ranking, original gene-disjoint d2, ternary rank margin m10 and elastic net alpha.5 / lambda.min. Preserve code and locked original manuscript unchanged.
- Let `eta=beta0+sum_j beta_j*r_j` with rank evidence `r_j in {-1,0,+1}`. At a prefix, `L=beta0+sum_observed beta_j*r_j` and `B=sum_unobserved |beta_j|`. The original RLR probability lies in `[sigmoid(L-B), sigmoid(L+B)]` for ANY ternary input, regardless of what test-order policy was learned.
- Early class 0/1 terminal is allowed only if **both** the whole posterior interval excludes 0.5 and its width ≤0.05, hence exact original RLR classes and maximum individual posterior absolute error ≤0.025. This remains true even for trees learned from training samples.
- New mechanism `R/experimental/rlr_certified_lookahead_tree_v1.R`: optimize a **restricted depth-two ternary decision-tree prefix** *globally over root choices and branch-specific second feature choices*, with coefficient-descending suffix. The root is selected from the top-10 nonzero |beta| relations; for each of its observed states -1/0/+1, exhaustively score each eligible second relation on the full downstream stopping path. Optimize training mean number of evaluated relations plus an explicit complexity penalty. **Not** unrestricted globally optimal tree induction and **not** a more accurate new classifier.
- The original RLR coefficients are fitted on outer TRAIN. Policy choice and acceptance depend on a nested D/V split of outer TRAIN only. To count as a real branching tree (rather than cosmetic reordering), at least two *observed* ternary states must lead to **different second tests**; and the V mean test-count saving over coefficient order must be ≥0.25. Otherwise fall back to the static coefficient-descending evaluation.
- Fair benchmarks: original RLR full evaluation, certified static |beta| ordering, prior adaptive *myopic* ternary next-test policy, proposed depth2 lookahead. Original class equivalence is guaranteed; the scientific outcome of interest is **how much additional conditional test-cost saving a nontrivial learned tree produces compared to static/mypopic baselines**, NOT BA uplift.

## Reproducible real evidence

[Two-cohort heldout GitHub run #40377950346](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40377950346). Three outer stratified folds each, seed 42. All splits and patient-level posterior intervals independently replayed by `scripts/rlr_certified_lookahead_audit_v1.py`.

**Pomeroy-2002-v1** — GitHub job SUCCESS, 3/3 folds and independent patient predictions PASS:

| Arm | Mean BA | Mean macro-F1 | Mean heldout logloss | Mean evaluated rank relations |
|---|---:|---:|---:|---:|
| Original RLR | .70370 | .72870 | .40323 | 36.667 |
| Static coefficient order | .70370 | .72870 | .39739 | 27.470 |
| Previous myopic context policy | .70370 | .72870 | .39780 | 27.497 |
| New restricted globally optimized depth2 tree | .70370 | .72870 | .39739 | 27.470 |

Static certified early stopping saved **25.1%** evaluated relational pair comparisons against all global-nonzero coefficients, without altering any original RLR classification. However the new depth-two learner **accepted zero of three conditional trees**, because D/V support did not show positive gain over static baseline meeting the preregistered 0.25 threshold. It correctly fell back to static in all three folds. Thus Pomeroy supports original certified **path pruning**, but **not** a tree-structure-specific gain.

**Alizadeh-2000-v1**: original-manuscript input, same frozen plan; when this document was drafted, GitHub's runner was still at original RLR dependency provisioning. **Do not report results until full audit PASS.**

[Exhaustive mathematical 81-state CI #40378410914](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40378410914) completed SUCCESS: enumerated **all 3^4 = 81 possible ternary inputs** for a four-pair model, verified class equality and posterior interval/prediction error for static and two independent genuine three-branch policies; checked train-only policy optimizer. This is a unit correctness proof-check on a finite fixture, **not** biomedical empirical validation.

## Novelty / limitations

Early-exit and cost-aware adaptive feature acquisition **already exist**, including binary-valued decision trees to reorganize weighted boosting decisions (Sun & Zhou, ICML 2013, https://proceedings.mlr.press/v28/sun13.html) and the broad taxonomy of acquisition/policy trees (Aronsson et al., IJCAI 2026, https://ijcai.org/proceedings/2026/860). A ternary relational early-exit *representation* is not automatically a publishable new interpretable tree classifier. Here the construction is strictly a hierarchical *inference-time policy for the original RLR*, not a model improving TSP/kTSP/TSPDT BA. It also does not reduce original biological assay costs: all original gene measurements and full-sample ranks are still required to evaluate any pair.

**Scientific stop rule:** If depth2 fails to beat coefficient-ordered fixed chain on BOTH cohorts, do not use this branch to claim tree-specific methodological novelty. Do not tune thresholds based on exposed E tests and call the result independent. The original unpublished RLR manuscript and core library remain unchanged.
