# CRIT Stage H3 v1 — blind XOR relation-pair retrieval via binary-signature LSH

Protocol frozen 2026-10-26. Corresponding [Stage H tracking issue #1](https://github.com/mczajkowskipb/relative-logistic-regression/issues/1).

## Evidence motivating the experiment

- Stage G/H0: at p32 the exhaustive 496-relational pair pool permitted orientation-invariant true XOR structure recovery in **12/12 developmental simulations**, with one-shot V acceptance and favorable test E logloss. This is not an independent biomedical estimate.
- Stage H1: under a 320-pair budget the XOR reference joint pair top-six rate degraded quickly: label-aware graph 7/8 at p32, 2/8 at p64, 1/8 at p128, 0/8 at p256, 1/8 at p512.
- Stage H2: if one **true** pair is supplied as an oracle, the other planted rank pair was ranked **#1 in 16/16** simulations for each of p32 through p1024, nD=120, under a supervised D-only exhaustive partner-pair scan. Sham anchor 0/8 top-six at all tested p. H2 therefore isolates an apparently difficult *first-anchor discovery* step but does not solve it.
- Neither marginal feature-class association nor independent random 320-pair proposals can be assumed to reveal a pure XOR in high dimensions. Naive exhaustive pairs-of-pairs would need `O(p^4)` scored combinations and be too expensive.

## Proposed scientific method — H3 binary signature LSH

A potential **blind** (no planted-gene information during feature discovery) retrieval scheme:
1. On *D only*, construct one binary vector per rank comparison `S_{ab}(i) = 1{X_ia > X_ib}` for all `M=p(p−1)/2` unordered gene pairs, retaining all vectors as a compressed boolean matrix. This construction is `O(n_D p²)`, still not deployable at transcriptome-wide p without prior dimension reduction.
2. Under a **pure XOR**, the observed binary label `Y` approximately satisfies `S_{ab} XOR S_{cd} = NOT(Y)` up to independent label noise (or `=Y` for reversed class sign), for **two disjoint candidate gene pairs**. Thus, if a candidate pair signature `S_{ab}` is known, the partner is close to `S_{ab} XOR NOT(Y)`. Crucially, this can be searched by a **hash join**, not by enumerating all M² rank-pair-pair combinations.
3. Use **24 independently seeded hash tables** selecting **18 D patients** each as binary positions. Hash all candidate pair signatures. For every candidate pair, query bins matching `hash(S_ab) XOR hash(Y)` and also its 18-bit complement; discard self-pairs, duplicates, and pairs sharing a gene. Evaluate exact `|2 * (# matching signs)/n_D − 1|` on all 120 D patients for retrieved candidate interactions. Keep the six highest-scoring unique unordered interacting pair pairs.
4. This is *entirely blind*: the true gene pairs (1,2) and (3,4) are referenced **only in the retrospective recovery evaluator**, never in hash construction, candidate retrieval, ranking, or tie-breaking.
5. Note label-dependent hashing and scoring on the same D sample is vulnerable to discovery overfitting. No V or E is used in H3. Formal validation must be done later with separate partitions and original-RLR OOF offsets.

## Fixed synthetic pilot

- `p ∈ {32,64,128,256,512}`, `n_D=120`, Bernoulli null and XOR `logit P(Y=1) = 2.4 sign(X1-X2)sign(X3-X4)`.
- Four independent generated D datasets per scenario × dimension = **40 D studies**. Seeds `20270211 + scenario_index*100003 + p*317 + replicate*43` for scenario_index null=1, xor=2; NumPy `default_rng`, standard-normal X, then Bernoulli Y.
- `HASH_BITS=18`, `HASH_TABLES=24`, exact candidate ranking through D-only `|corr((2Y−1), signpair_a*signpair_b)|` (binary correlation, not Pearson normalization).
- Save all study-level search counts, true collision and true top-six recovery, top-six actual pair IDs/scores, elapsed time, and memory. No known anchors, no original RLR fit, no D/V/E split, no protected biomedical data.
- GitHub-hosted Actions runs after branch addition; isolated archive with run SHA, scientific protocol, CSV and independent Python audit.
- Independent audit regenerates *all 40 X/Y datasets* from the stated NumPy RNG seeds and recomputes all saved top-six pair-product scores independently of the LSH algorithm, verifies orientations, disjointness, rank, truth flags and provenance. It cannot certify *all* hash collisions were enumerated, so unit tests or differential exhaustive p32 comparison are still necessary before scientific claims.

## Interpretation and limitations

- Successful exact recovery in XOR with few/no recovery events for arbitrary reference pairs under null would demonstrate **initial feasibility of a blind candidate generator** for this artificial high-signal case. This is not type-I error control, not a logloss gain, and not a calibrated test. A random reference pair not entering top-six under null cannot be equated to low false-positive rates.
- Failure is equally informative: 24 × 18-bit hash settings, strong noise, finite nD, dependence among pair features, and hash joins can yield missed true collisions or high-score false relations. Do not tune hash settings to improve the same 40 D outcomes and then claim a new independent test.
- `M=O(p²)` rank-pair signatures; hash joins can have large collision fan-out. Runtime and memory must be measured and bounded. At full omics dimensionality, even `O(p²)` is too large; any future deployable version must use an *independent* label-safe dimension-reduction or candidate-generation strategy.
- H3 is a **synthetic proof-of-concept**, not a validated biological classifier. Do not merge into original-RLR package or retrofit published results. Subsequent method evaluation requires strict nested/independent splits, patient-level test E predictions, real original RLR OOF offsets, explicit multiplicity calibration and independent unexposed cohorts.

## Required PASS

`STAGEH3_40_BLIND_D_ONLY_STUDIES=PASS`, `STAGEH3_24_TABLE_LSH_XOR_PAIR_DISCOVERY=PASS`, `STAGEH3_TOP6_PROVENANCE=PASS`, `STAGEH3_240_TOP6_SCORES_RECOMPUTED_FROM_FRESH_D=PASS`, `STAGEH3_40_STUDY_BLIND_SOURCE_INVENTORY_AUDITED=PASS`, `STAGEH3_BATCH_EXIT=0`.

No external users' or clinical data. All archives remain independent of previous Stage H0/H1/H2 versions.
