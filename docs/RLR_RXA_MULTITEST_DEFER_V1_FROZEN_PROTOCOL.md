# RXA MultiTest with rank-margin deferral — V1 frozen exploratory protocol

Date 2026-10-31; fixed **before V1 outcome retrieval**. PI's goal is an alternative internal DT test that combines historical RXA multi-test nodes, rank-margin ternary states, and (optionally) global logistic rescue for deferred samples. This is **NOT** a new 'full decision tree' or an unqualified novelty claim. It is a **depth-two cascaded relational node test**, useful as a candidate for DT integration (RMCT / EMTTree+RX) after source and control verification.

## Primary V1 questions, not guaranteed positives

1. Does rank-margin abstention in a fixed **3-pair multi-test at each of two sequential levels** yield better predictive scores than exactly the same pair-voting scheme with binary pair-order tests?
2. At fixed margin and routing, does replacing an uncertain terminal smoothed-prior leaf with **a global, originally trained RLR20 logistic model** improve held-out predictive quality? This global model is trained on *all outer-training cases*, never on the small deferred leaf.
3. Does the cascade reduce evaluated tests (median/mean/p95 per patient), **including** the global logistic cost if rescue is used?

No RLR-DT teacher distillation, no per-branch RLR retraining, no supervised outer E threshold tuning, no Monte-Carlo split rerolling.

## Frozen data / source parity

First pilot on all eleven public Schliep reconstructed, original-paper-named CompCancer datasets, **NOT original ARFF SHA-verified**. 5 stratified outer folds, original seed42, fixed first repeat. All methods receive identical train/test splits, features, and a **single outer-T-only switchBox top1000 pool and K20/d2 gene-disjoint shortlist**. First 6 pairs are split into two fixed ordered 3-pair node multi-tests. No filtering by outer E outcomes.

The owner verified SHA256 parity for the original author-held CompCancer ARFF files on Ubuntu in a prior separate experiment. That does **not** turn public reconstructed data used by CI into original SHA-verified files. V1 is exploratory on previously inspected 11 cohorts. Should V1 show a meaningful gain, independently rerun the **same frozen algorithm** on originals using a new gated runner, then decide whether any exploratory supplementary inclusion is justified.

## Fixed routing and probability outputs

For pair k, let Δ_k(x) be the within-patient rank(B)-rank(A). The binary node state is sign(Δ), with ties producing 0. Ternary state is +1 if Δ>=10, -1 if Δ<=-10, else 0. A pair's class-vote orientation is estimated **only on outer-T** from difference of class-average signed binary ordering, with zeros assigned orientation +1 (deterministic). At each 3-pair node, assign POS if >=2 oriented +1 votes and **zero** -1 votes; NEG if >=2 oriented -1 and **zero** +1 votes; otherwise DEFER (including conflicts, insufficient margin, and ties).

- On node 1 POS/NEG, emit node1 smoothed class frequency estimated only from its T-assigned samples.
- On node1 DEFER, evaluate node2 with *different next 3 gene-disjoint pairs*. On node2 POS/NEG, emit smoothed class frequency estimated only from T samples deferred by node1 then assigned by node2.
- On node2 DEFER, either (a) emit T-global-smoothed class prior (ordinary unresolved leaf), or (b) emit probability of the **full original RLR20** model fit on all outer-T samples, which is a real logistic backstop for the ambiguous region.
- Fixed prior smoothing weight = 10 globally to guard tiny leaves; min 5 training samples for terminal empirical leaf use, otherwise global prior. Do not adjust after inspecting E.
- The standard logistic backstop means each node is a multi-test, not a branch-specific retrained glmnet expert, and no duplication of observations occurs.
- Pair orientation and class priors are fitted once on all T. No outer E labels are read or used before scoring.

## Six comparator arms on identical outer folds

1. ORIGINAL_RLR20: frozen original ternary m10/d2 alpha=.5/lambda.min logistic model, all K20 pairs.
2. BINARY_MT2_LEAF: fixed binary 2-node/3-pair tests, terminal smoothed class-rate leaves.
3. TERNARY_MT2_LEAF: same fixed tests but margin m10 and zero as abstention.
4. BINARY_MT2_RLR_BACKSTOP: same binary nodes, global logistic only for last unresolved state.
5. TERNARY_MT2_RLR_BACKSTOP: margin nodes, same global logistic backstop.
6. TERNARY_SAME6_ENET: elastic-net logistic model on exactly same first six original ternary gene-disjoint pairs, to control for multi-test routing with equal feature set and no tree.

Measured per sample: probability, classification, route at which prediction was made, gate-pair evaluations 3/6, full RLR20 pair evaluations if backstop applied, total pair checks, whether backstop used. **RLR20 constant 20 pair checks** counts model feature evaluations, conservatively even if some coefficients zero. For TERNARY_SAME6_ENET use six pair checks. Since rank-based margin requires the full measured feature universe to rank, no claims of 3/6-gene wet-lab assay cost.

Primary outcome = mean dataset-equal heldout logloss. Secondary: Brier, AUC, F1, mean pair checks and backstop rate. Main comparison TERNARY_MT2_RLR_BACKSTOP vs BINARY_MT2_RLR_BACKSTOP (isolate ternary gating holding logistic availability fixed), and TERNARY_MT2_LEAF vs BINARY_MT2_LEAF (isolate ternary gating without logistic). The other models are sanity/routing controls. Paired inference uses exactly **eleven datasets**, not 55 folds; exploratory exact two-sided sign flips 2^11 for the two fixed main contrasts and two metrics (LL, F1), BH across the four tests. Report all null and negative outcomes as well as positive. All 55 folds + patient-level audit gate PASS required before numeric claim.

## Preventing novelty exaggeration

The original 2014 MTDT and 2025 RMCT/EMTTree+RX already have multi-test nodes (published by same research line). Three-state rules, reject/abstain options, early-exit cascades, partial evaluation and mixture-of-experts are also prior art. This v1 only tests whether a **rank-margin deferral mechanism with coherent logistic rescue** adds validated value over appropriate within-research controls. It is not yet an EMTTree induction algorithm, a truly fuzzy tree or proof of superior performance to full historic GPU RMCT/TSPDT/EMTTree implementations.

## Publication gate

Current CMPB paper (RLR as ternary penalized logistic method) stays **unchanged** while V1 runs. A new tree section requires genuine improvement in a stated tradeoff and explicit comparison against matching binary multi-test controls; do not sacrifice paper submission deadline for negative/ambiguous results.
