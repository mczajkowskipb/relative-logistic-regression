# RLR tree extension gate: ARLET V1 vs certified relational early-exit tree

**Date:** 2026-10-26. The unpublished original manuscript is the target for an eventual single, tree-relevant RLR paper, but **no changes to locked original RLR source/manuscript** are authorized by these exploratory tests. All test data are from the original-paper cohort inventory (Chen-2002, Singh-2002, Yeoh-2002-v1), **reconstructed public Schliep source**, not SHA256-confirmed original processed ARFF. Three heldout stratified outer folds per cohort, seed42, repeated-fold 11-world original benchmark **not** reproduced here. Independent E-patient predictions and audits are archived.

## Hypothesis A: genuine three-way conditional ARLET, large pair pool, deterministic vs EA

The new experiment `analysis/rlr-arlet-evidence-tree-pilot-v1` implements a 3-way conditional evidence node: confident class-0, confident class-1, or unresolved patient requiring an additional **distinct** 7-relation multitest. Original `m=10` signed ternary margins. Root test has 3 or 5 relations. Second relation set is nominated on difficult *training* cases, while logistic coefficients at both stages are fitted on the **entire** training part, never tiny leaves. Candidate pool = original top-1000 switchBox supervised pairs plus about 250–450 deterministic seeded gene-pair samples from the full feature universe. The inner D/V split is fully inside outer T. Deterministic roots come from signed ranking; bounded actual EA of 3 generations×8 explores split relations and confidence thresholds only on D/V and uses heldout E only for a single final evaluation.

Controls call untouched original `fit_rlr_train_test`, `fit_tsp1_manual_train_test`, `fit_ktsp_train_test` (`k=9`); the supplied CUDA TSPDT is NOT being reproduced by the new `TSPDT_GREEDY_PROXY`, which is only a simplified 2-level greedy single-relational-tree control. No parity claim with exact C++/CUDA TSPDT.

[ARLET GitHub Actions #40366023886](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40366023886)

**Independent-completed cohorts at this note's time:**

| Dataset | Model | Mean BA | Macro-F1 | Logloss |
|---|---|---:|---:|---:|
| Singh-2002 | Original RLR | .91231 | .91061 | .33194 |
| Singh-2002 | kTSP-9 | .92211 | .92147 | .28230 |
| Singh-2002 | TSP-1 | .82394 | .81980 | 2.84095 |
| Singh-2002 | TSPDT GREEDY PROXY | .84300 | .83610 | .40348 |
| Singh-2002 | ARLET deterministic | **.93137** | **.93141** | **.24149** |
| Singh-2002 | ARLET EA | .92150 | .92070 | .27664 |
| Yeoh-2002-v1 | Original RLR | **.92446** | **.93957** | **.10106** |
| Yeoh-2002-v1 | kTSP-9 | .86829 | .73184 | .45874 |
| Yeoh-2002-v1 | TSP-1 | .82791 | .73208 | 3.39148 |
| Yeoh-2002-v1 | TSPDT GREEDY PROXY | .76660 | .80822 | .25307 |
| Yeoh-2002-v1 | ARLET deterministic | .74524 | .79654 | .25352 |
| Yeoh-2002-v1 | ARLET EA | .72222 | .77397 | .24178 |

**Scientific interpretation:** Positive signal on Singh, substantial failure on Yeoh (both discrimination and logloss). Not sufficiently robust to add ARLET v1 to manuscript. EA did not improve deterministic ARLET on Singh or Yeoh. In addition, TSP-1 raw hard probabilities of 0 or 1 produce very large logloss by design, so its logloss should NOT be used to claim a calibrated superiority; BA and macro-F1 are informative.

**ARLET Chen**: at the moment of this report, GitHub runner remained in apt dependency provisioning; do not fabricate result or count it as validation failure.

## Hypothesis B: exact-class certified RLR evidence evaluation as a structured three-way early-exit tree

Instead of fitting a new tree and giving up RLR quality, derive a sequence of ternary relational tests from *the literal original published RLR coefficients*. Let `eta(x)=b0+sum b_j*r_j(x)`, `r_j in {-1,0,+1}`. After a prefix `S`, partial logit `L=b0+sum_{j in S}b_j*r_j`, remainder absolute bound `B=sum_{j not in S}|b_j|`. Therefore `p_original(x) in [sigmoid(L-B),sigmoid(L+B)]` **for every patient and any future unobserved relation**, independent of phenotype or test data.

An early exit is allowed only if (a) this interval stays on one side of .5, so the **predicted class is exactly identical to the full RLR**, and (b) interval width ≤0.05, so the interval midpoint `p_approx` has a mathematically verified maximum absolute error ≤0.025 to the full RLR probability. The rest follow an additional relation. This produces a tree-shaped stop-class-0/stop-class-1/continue sequence, not an independently learned branching classifier. Order: descending |coefficient|, and an optional 4-generation×8 EA (training-only) minimizing expected number of evaluated relations. No invented individual gene causal explanation from this screening.

[Certified RLR GitHub Actions #40367860184](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40367860184). All 3 datasets SUCCESS; all patient-level audits PASS, 3 outer E-fold partitions per dataset.

| Dataset | Mean originally needed nonzero-coefficient pairs | Mean evaluated (certified) | Early-stopped E patients | BA (both; identical) | Approx logloss | Exact RLR logloss |
|---|---:|---:|---:|---:|---:|---:|
| Chen-2002 | 20.000 | 12.040 | 99.44% | .91291 | .21154 | .19818 |
| Singh-2002 | 20.000 | 13.157 | 100.00% | .91231 | .29283 | .33194 |
| Yeoh-2002-v1 | 8.667 | 4.971 | 87.86% | .92446 | .09878 | .10106 |

Mean evaluated reductions: Chen 39.8%, Singh 34.2%, Yeoh 42.6%. Mean **active** relations evaluated: Chen 9.72 vs RLR 16.24; Singh 11.39 vs RLR 17.49; Yeoh 3.21 vs RLR 5.68. These are patient-conditional rank-margin-active evidence counts, NOT costs of acquiring raw omics input data. Approximate logloss can increase or decrease; the bound guarantees maximum per-patient absolute posterior error, **not lower logloss**. EA reordering was almost identical to descending |b| in early-exit cost, and is not essential.

**Key scientific limitation:** A fixed-order early-exit linear-score computation is not automatically a novel decision tree learning algorithm. Existing research on sequential probabilistic prediction and anytime feature sequencing is relevant. For a hab-aligned method extension, need *genuinely state-dependent branching* or locally learned routing policy built using the ternary inactive state, with comparison against early-exit linear baselines and method-specific interpretability evidence. The strong provable bound and exact-class property make this a stronger mathematical starting point than a fragile new leaf classifier, but originality is **not established**. The current results demonstrate real engineering/interpretability savings, not predictive superiority or published-method novelty.

## Next scientific gate, before any claims

1. Do **not** tune the failed ARLET confidence thresholds/selection rules using the exposed Singh/Yeoh heldout cases and then call that independent proof. Finish Chen if its CI completes.
2. The certified, exact-class structure is the most defensible next idea for extending the current unpublished RLR manuscript, but only if it becomes a genuinely conditional/learned multi-test tree distinct from trivial early-exit linear evaluation; verify published prior art.
3. Any next method prototype should be fixed by training-only rules and evaluated on other **original eleven CompCancer datasets** not used in this ARLET/C-RET pilot, then the complete original published protocol with exact preprocessed ARFF parity if available.
4. Do not claim replicated published results from these reconstructed public-source matrices, nor original TSPDT equivalence to the greedy proxy.
5. Do not modify the locked original RLR manuscript until a genuine, empirically successful and defensibly novel algorithmic extension exists.
