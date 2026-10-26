# Stage H1 v1 — bounded discovery-only interaction screening feasibility

Protocol version fixed 2026-10-26. Parent task: [CRIT Stage H issue #1](https://github.com/mczajkowskipb/relative-logistic-regression/issues/1). Frozen H0 audit source SHA `21986f8ece37ae5056ea1123360e12c74bcf60d3`. Frozen Stage G source `6ca44c30a439dbdb97b3a45f7e8f45432da3bd59`.

## Facts motivating H1

- The Stage H0 independent symmetry audit PASS verified **5 direct + 7 swapped = 12/12** exact planted XOR recoveries in the completed **32-gene full-496-pair** Stage G ablation; all 12 conditional corrections accepted and improved test logloss. This is **development evidence**.
- The original CRIT-v9 `crit7_pool()` uses the fixed proposal seed `cfg$seed+617` for random pairs and omitted the planted gate in all six historical synthetic XOR replicates. Thus the six misses were not six independent random proposal draws.
- Screening every pair of high-dimensional genes, and then every gate×relation combination, scales at least quadratically and may be prohibitive in transcriptomics. A feature-level supervised interaction *retrieval* heuristic is worth studying but must not be confused with a validated predictor.
- Pure XOR can have zero marginal class association for each participating rank pair. **No marginal-only gate ranking can be assumed to retrieve it.**

## Specific questions and scope

1. Under a fixed cap of **320 candidate rank pairs**, how often does uniformly random sampling include *both* exact planted pairs? How much do results depend on a fixed seed versus independent per-study pair sampling?
2. Can a lightweight **D-only second-order rank-covariance graph** increase coverage, without exceeding 320 candidate pairs, and how does coverage deteriorate as dimension grows?
3. Conditional on both planted pairs being included, does a **simplified six-position rank-pair×rank-pair residual interaction shortlist** contain the exact planted pair product in either order? This is *NOT the full `crit7_nominate` routine* (which also constructs two leaf-specific relations and applies support filters), and this result is only a screening feasibility metric, not model identification or acceptance.
4. What are the pair budgets and *observed screening run times*, with original RLR and all historical files untouched?

## D-only synthetic benchmark design

- **Dimensions:** p = 32, 64, 128, 256, 512. **nD=120** independent samples (the typical Stage-G discovery size).
- **Scenarios:** Bernoulli null, and XOR where `z = 1{gene001 > gene002}`, `r = 1{gene003 > gene004}` and `logit P(Y=1)=2.4 z_signed*r_signed` with signed values ±1.
- **Replicates:** 8 predeclared seeds independently generating sample matrices/outcomes for each p×scenario, hence **80 distinct D-only studies**.
- **Methods:** (i) FIXED_RANDOM: exactly 320 unordered gene pairs (uniformly sampled with fixed seed 659 across replicates); (ii) SEED_RANDOM: exactly 320 unordered pairs with predeclared per-study seed; (iii) AWARE_GRAPH: build *per-patient feature ranks*, center each gene-rank across D samples, L2-normalize, take absolute outcome-residual weighted gene covariance `|X_rank^T diag(y−0.5) X_rank|`, accumulate node degrees over top 128 weighted cross-gene edges, retain the 24 highest-scoring genes, include *all 276 within-node pairwise comparisons* and complete to 320 with per-study uniformly sampled other pairs. These are **not the same as original v9's RLR-conditioned candidate pools**. (iv) EXHAUSTIVE32: all 496 rank comparisons *only for p=32* as a gold-standard retrieval reference.
- All methods compute relation signatures `sign(gene_i−gene_j)` within D and rank joint interactions by `|Σ_D (y−0.5) sign(pair_a)sign(pair_b)|`, requiring disjoint genes, considering both gate–relation orientations equivalent. Keep **top 6 unordered interacting pair-pairs**, without any V or E data. The unconditional marginal baseline `0.5` is known for the synthetic generator, **not estimated original-RLR OOF**. This is a deliberately simplified retrieval experiment, and no RLR test, conditional model, calibration or patient E endpoint is evaluated.
- **Workload:** 80 D studies × 3 capped methods + 16 p32 exhaustive controls = **256 screening runs**. Candidate-pool archives and top-six candidate IDs/scores are retained for each run. All method counts and pair budgets checked by an independent Python structural/provenance auditor.
- **Source constraints:** no changes to original published RLR, CRIT-v7–G, H0; no overwriting old results; no use of patient E test labels; deterministic seed manifest and new directory only.

## Scientific interpretations, strictly limited

- Failure of exact pair coverage means the method **cannot identify the planted exact XOR interaction** in that run, regardless of downstream fitting and V gate.
- Coverage but failure in top-six shortlist shows a distinct **screening/ranking** bottleneck under the surrogate D-only score.
- Success means only that the true pair-product entered a *D-only surrogate shortlist*. It is **not evidence that the original RLR-conditioned CRIT, its leaf-selection logic, or V acceptance would recover it**.
- The label-independent null has the same gene-pair references, but these are *arbitrary reference pairs*, not causal truth. Null candidate scores, by themselves, **do not estimate false acceptance**.
- Larger p can cause serious degradation; null/XOR samples are small (8 each). Any successful method must next be integrated into the exact original-RLR offset pipeline with held-out V and E and many more independent null studies.
- Time measurements are machine-specific and not a standalone runtime benchmark; memory, BLAS and resource profiles should be captured separately before larger p.

## Outputs / required gates

New output: `results/crit_stageH1_bounded_screening_v1/`:
`screening_summary.csv` (256 rows), `candidate_pools.csv` (84,736 rows), `top6_pair_interactions.csv` (1,536 rows), `256_screening_runs.log`, `independent_audit.log`, SHA and PROTOCOL.

Runtime gates:

`STAGEH1_80_INDEPENDENT_D_ONLY_SYNTHETIC_STUDIES=PASS`,
`STAGEH1_256_SCREENING_RUNS_COMPLETE=PASS`,
`STAGEH1_SOURCE_POOL_AND_TOP6_PROVENANCE=PASS`,
`STAGEH1_256_INDEPENDENT_PROVENANCE_CHECKS=PASS`,
`STAGEH1_D_ONLY_PAIR_BUDGETS_AND_SYMMETRIC_JOINT_AUDITED=PASS`,
`STAGEH1_BATCH_EXIT=0`.

No scientific selection or algorithm retuning on archived Stage-G E outcomes. H1 is a capacity/failure-mode pilot, not a finalized CRIT-v10 implementation.
