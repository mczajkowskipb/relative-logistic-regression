# Anchor Look-Ahead with Child RLR Experts — V1 predeclared protocol

**Aim:** Test whether choosing a rank-invariant anchor root by the *predicted improvement available to downstream child RLR models* is more effective than choosing the root for its immediate class separation.

**Original manuscript and original source RLR unchanged.** This is a development hypothesis on **already inspected** original RLR-paper cohort reconstructions; no independent validation or research novelty claim.

## Invariance contract and architecture

Each internal test is `S_g(x) = sum_{a in A} 1[x_g>x_a] >= t`, with one target hub gene g, reference panel A of 3 or 5 anchor genes, and integer threshold t in 1..|A|. The entire representation and every downstream RLR pair feature are within-patient rank-order tests; the same patient can be independently transformed by any strictly increasing transformation without changing these features, provided that the measured gene/probe universe and rank tie structure are unchanged.

A patient travels down one root branch b in {0,1}. The global original RLR supplies a fixed logit offset `eta0(x)`. The two child experts apply branch-specific regularized relational corrections `delta_b(x)`, so `eta(x) = eta0(x) + delta_{g(x)}(x)`. Each `delta_b` is a small model of at most 8 selected ternary-rank gene-pair predictors `m=10`, gene-usage constraint `d=1`, fit with `glmnet` binomial alpha=0.5 and fixed lambda=0.10, intercept=FALSE and the original RLR as offset. When a local fit is not feasible, the branch predicts literal parent RLR. **This is one decision-stump model tree with two local RLR experts, not an ordinary class-constant leaf tree and not a demonstrated depth-two tree.** The shared original RLR offset is not refit separately inside each leaf.

## Bounded look-ahead and honest partition roles

For each fixed 3-fold outer test E:

1. Split the outer training T into D and V with the pre-existing fixed `rg_inner` protocol. E is not read until all root/validation decisions are made.
2. D-only cross-fit the *original RLR* using three stratified folds. Each heldout D patient's RLR prediction comes from a model trained without that patient, including fresh supervised switchBox pair selection per OOF fold. These predictions supply crossfitted residuals `y - p_OOF`.
3. Nominate at most 12 target genes by occurrences among D top TSP pairs; reference panels of sizes 3 and 5 come from D TSP neighbors and rank-stable D candidates using the existing frozen star heuristic. Candidate thresholds 1..k must yield at least 12 D patients and at least three of each class in each branch.
4. Evaluate **cheap contextual second-order gradient scores**, using D OOF original-RLR residuals and first 160 D top scored ternary rank pair columns: best improvement available from **two separate branch-specific coefficients** minus the approximated improvement from one common coefficient. The approximation is a ranking filter only, not an unbiased out-of-sample estimate. Retain at most **six** root candidates.
5. For each retained candidate, perform *real branch-level supervised pair selection on D patients in each branch*: switchBox top600, disjointness d1, K8, m10. Fit the two child RLR corrections on D, using a parent RLR trained exclusively on D as an offset. Skip roots without two successfully fitted branch models and any nonzero added coefficients.
6. Predict D-independent V with each candidate, choose the lowest V logloss among the six maximum proposals. Compare on exactly the same V subjects against:
   - the D-trained original RLR K20 (d2/m10/alpha0.5/lambda.min);
   - a gate-only correction atop the parent RLR (root split but no conditional rank features);
   - a same-pairs global non-tree elastic-net correction using the union of all child-nominated pairs plus root gate, without branch-dependent coefficients;
   - a **greedy root** selected on D for its own class-separation leaf loss (same root panel heuristic) and fitted with the same child-RLR experts; this comparator is reported even if no lookahead root passes V.
7. **Accept** the lookahead tree only if V logloss improvement is at least 0.008 versus original RLR, 0.004 versus same-pairs global, 0.004 versus gate-only, and V balanced accuracy loss no worse than 0.02. Report V improvement or loss versus greedy even when it fails acceptance. Otherwise **exact original RLR fallback** on E.
8. Refit original RLR and the D-selected branch pair coefficients using all outer T, with root structure and pair identities fixed from D; predict E only once. For a fair comparison the greedy root ALSO nominates its child pairs on inner D, and only their coefficients are refitted on outer T; the identities of greedy and lookahead child pairs are frozen from D. Predeclared outputs include full diagnostics even if 0/6 accepted.

## Computational upper bound

At most six D root expert fits per outer fold, each with two branch correction models (12 bounded child fits) plus one greedy-root comparator. Three D OOF original RLR fits, one D original RLR, one T original RLR and independent same-pair/global control. Candidate root scoring is vectorized over at most 160 precomputed ternary rank pair columns. No exponential feature-set search, EA, model-in-leaf CV or brute-force search of all possible gene subsets.

## Auditable comparator grid

- ORIGINAL_RLR
- GREEDY_ROOT_CHILD_RLR
- GREEDY_SAME_PAIRS_GLOBAL (greedy D-selected child pair union, but one shared global coefficient vector)
- GREEDY_GATE_ONLY (the same greedy root without local gene-pair corrections)
- LOOKAHEAD_CHILD_RLR (unconditionally fitted, diagnostic)
- LOOKAHEAD_SAME_PAIRS_GLOBAL (same child gene-pair union; no branch-dependent coefficients)
- LOOKAHEAD_GATE_ONLY
- V_GATED_LOOKAHEAD_CHILD_RLR (actual inner-selected decision)

Development datasets: already-exposed Chowdary-2006 (104 patients/182 genes) and Singh-2002 (102/339), source-reconstructed public matrices not checksum-confirmed to original processed ARFF. Outer 3-fold stratified E seed42 × one repeat => 6 E partitions, 48 model-fold metric records, 1648 heldout patient-method scores. Primary analysis Δ logloss versus original RLR; secondary Brier, BA, Macro-F1, number of fit local experts and nonzero relational coefficients, accepted root count, comparison against greedy/same-pair-global/gate-only. Independent Python scorer verifies all E rows, original class coding, patient separation/coverage, metric parity, and exact V-gated fallback.

**Stop condition:** A technically successful but unfavorable real-data result is archived as negative. These E folds have already been used repeatedly; do not modify the root screening, lambda, accept cutoffs, or support conditions after examining these heldout results and claim validation. If the child-RLR mechanism is infeasible or has no accepted roots, future research may probe conditional train-only signal and potential independent data, not sweep E. Original finished RLR paper stays intact; real CUDA TSPDT is not part of this experiment.

## Audit amendment — comparison fairness, prior to complete eight-arm audit

An intermediate six-arm run revealed that the greedy-root competitor selected its local pair *identities* anew on full outer T, whereas the lookahead root froze local pair IDs from D. This asymmetry was corrected before the complete eight-arm run: **both freeze candidate root and child pair IDs on D, then refit coefficients on T**. Two further logically necessary greedy-root controls (global union of exactly the greedy pair features, and gate-only) were added without altering lookahead feature nomination, penalties, threshold or V acceptance. The complete eight-arm design thus has 48 fold-method scores and 1648 heldout patient-method predictions, individually independently audited. These are **methodology repairs and controls on already exposed development E**, not an independent pre-registered confirmatory hypothesis.
