# CRIT Stage H2 — oracle-anchored second-pair search, discovery-only

Protocol fixed 2026-10-26. See [Stage H issue #1](https://github.com/mczajkowskipb/relative-logistic-regression/issues/1).

## Why this study

Frozen Stage H0 established **12/12 orientation-invariant planted XOR pair recoveries**, accepted + test gains on a 32-gene exhaustive 496-pair discovery pool. Stage H1 independently confirmed that three fixed-size 320-pair proposals degrade sharply when the number of genes increases: awareness-graph true joint screening top6 was 7/8 at p32, 2/8 at p64, 1/8 at p128, 0/8 at p256 and 1/8 at p512. Stage H1 is a *D-only surrogate screen* under the known synthetic residual y−0.5, not a full RLR model.

It remains unclear whether the hard part at high p is (a) identifying **any** suitable seed pair, or (b) finding its interacting partner after a correct seed is available. Stage H2 isolates (b), using a scientifically explicit **oracle anchor** that would not exist in a blind clinical or omics deployment. It does not claim to solve (a) and is not a deployable estimator.

## Predeclared synthetic design

- **p** = 32, 64, 128, 256, 512, 1024; nD = **120**; **8 independent replicates** each for null and XOR; hence **96 independent synthetic discovery datasets**, **288 anchored scans** (three anchors per dataset).
- New seeds: `20270117 + scenario_index*100000 + p*311 + replicate*37` (scenario indices null=1, xor=2). Data X are independent continuous Gaussian genes; planted signed rank comparisons: `z = sign(X1-X2)`; `r = sign(X3-X4)`; null `logit P(Y=1)=0`, XOR `logit P(Y=1)=2.4zr`. **All D labels only**, no V or external E set.
- Arms, all predeclared: **ORACLE_Z_TO_R** with supplied pair (1,2), target (3,4); **ORACLE_R_TO_Z** with supplied pair (3,4), target (1,2); **SHAM_56_TO_R** with irrelevant anchor (5,6), reference target (3,4). A sham arm is an ability-to-discriminate reference, not a formal null-statistical test.
- For each anchor, enumerate **every gene-disjoint unordered target pair** in blocks, but compute only its sign relation and scalar score `abs(sum_D ((y−0.5)*sign(anchor))*sign(target))`. Max temporary matrix 120×2048. Candidate set size = `(p−2)(p−3)/2`. Select up to six best-ranked candidate pairs with deterministic tie-breaking; report planted target rank, top1, top6 and elapsed time.
- Important computational limitation: **this is blockwise O(nD p²) search, not O(p) or subquadratic.** At p1024 it scans **521,731** candidate pairs per anchor. It avoids storing a p²×nD array but may still be slow. Observed runtime is machine-dependent, and the complete-blind CRIT search still lacks a feasible oracle anchor.
- Analytical baseline: exact probability that a random **320-pair** pool from all `M = p(p−1)/2` unordered gene pairs includes both planted pairs is `320/M × 319/(M−1)`, regardless of outcomes. Compare this baseline to Stage H1 coverage, without equating it to a test-set success probability.

## Scientific interpretation

- If true oracle anchor finds its partner near rank1 at high p while the sham arm fails, the main bottleneck is **finding a first informative rank pair**, not conditional partner completion. Next research must identify seed pairs without using planted truth, outcome leakage, or private data.
- If oracle anchor itself fails with nD=120 for p512–1024, large multiple-comparison burdens and sample size affect even privileged conditional search. Investigate sample-size scaling and control adaptive discovery bias before model integration.
- No original RLR parent fit, original RLR OOF, V acceptance, E performance, patient-level predictions, class error control or high-dimensional generalization claim. This experiment is **only synthetic D-screening diagnostics**.
- No clinical or biomedical data; safe to run on GitHub-hosted actions. Do not attach an unrestricted runner with access to sensitive local data to the current public repository.
- All source and outputs versioned; Stage H1/fix1 archives and original RLR untouched.

## Required evidence and outputs

`results/crit_stageH2_anchored_completion_v1/`: `study_summary.csv` (288 method-dataset rows), `top6_candidates.csv` (1,728 rows), `run.log`, `independent_audit.log`, `SOURCE_SHA.txt`, `PROTOCOL.txt`, `status.txt`; uploaded GitHub Actions artifact.

Required PASS:
- `STAGEH2_96_D_ONLY_SYNTHETIC_DATASETS=PASS`
- `STAGEH2_288_ANCHORED_SCANS=PASS`
- `STAGEH2_288_PAIR_AND_METRIC_PROVENANCE=PASS`
- `STAGEH2_ORACLE_VS_SHAM_TARGET_IDENTITIES_AUDITED=PASS`
- `STAGEH2_BATCH_EXIT=0`

The independent Python audit verifies published study/seed/pair IDs, target rank vs top6, score monotonicity, source cardinalities, and random analytic benchmark. It does **not** claim independent reconstruction of all nD × p² R scores because the entire D matrices are not archived.
