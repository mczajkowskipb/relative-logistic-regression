# RLR rank-pair contrast stump — completed audited experiment (2026-10-26)

Real GitHub run SUCCESS: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40421745933
Fitted code SHA: 3d90ac3858f78f88f2a4579d658ad92ff11b52cd
Original RLR source and locked manuscript are unchanged. Development data were already exposed, so these results are exploratory only.

## Independent, rank-preserving root construction

A root compares TWO separate gene-pair rank distances in a single patient:

delta_AB = rank(B) - rank(A);
delta_CD = rank(D) - rank(C);
gate = 1[delta_AB > delta_CD].

This is ARXA-inspired ordinal rank-contrast, NOT an expression ratio, raw expression threshold, or cohort scale threshold. It stays invariant to a strictly increasing transformation of all features inside each sample at fixed feature schema. Root pair candidates sample the full rank-only candidate pool including high-MAD-of-ranks expansion (not raw feature MAD). Screening ranks root-pair contrasts by interactions with one ternary original-RLR m10 partner. Root gene pairs are allowed to share genes in this pilot, limiting claims that every contrast represents four independent genes.

The parent is the existing original RLR configuration K20, d2, m10, alpha0.5, lambda.min. Non-tree same-pair ridge global correction and stump use the same parent offset, and the stump alone adds the rank-contrast gating interaction. A fixed D selection / inner V acceptance process (>=.005 mean logloss vs RLR, >=.002 vs global no-tree same-feature control, BA decrease no worse than .02) selects or returns exact original RLR. No training or acceptance ever accesses outer E.

## Mechanical proof and completed execution

The new R unit test confirms invariance to independent increasing affine transformations of 120 patient profiles, unchanged rank-contrast root and ternary partner design, and valid proposal discovery on an artificial rank-context target. Previous base-R outlier/schema/invariance tests also PASS. 

Actual original-paper Chowdary-2006 and Singh-2002 source reconstructions were downloaded and audited. Six fixed 3-fold heldout partitions total, with 24 method-fold metrics, 824 patient-method test scores. Independent Python replay of logloss, Brier, BA, patient ID uniqueness, entire heldout cohort coverage, class coding and exact fallback passed.

## Per-fold internal V gain versus external E diagnostic

| Dataset | Fold | Inner V gain | V accepted? | Ungated E delta logloss |
|---|---:|---:|---|---:|
| Chowdary | 0 | +.008072 | YES | +.005163 |
| Chowdary | 1 | +.003278 | no | -.005699 |
| Chowdary | 2 | +.000973 | no | +.018979 |
| Singh | 0 | +.007506 | YES | +.000761 |
| Singh | 1 | +.000672 | no | -.004555 |
| Singh | 2 | +.035200 | no (other condition) | -.012106 |

Both accepted stumps worsened E logloss. A large V gain by itself did not guarantee gate acceptance because the same-pair non-tree control and BA constraints were also enforced.

## Six-fold-unweighted mean heldout E delta logloss vs original RLR

| Variant | Mean delta logloss | Improved folds |
|---|---:|---:|
| Original RLR | +.000000 | 0/6 |
| Same-feature global no-tree correction | +.000271 | 3/6 |
| Ungated rank-contrast stump (diagnostic) | +.000424 | 3/6 |
| V-gated rank-contrast stump (actual selected model) | +.000987 | 0/6 |

**Result: technically successful, scientifically negative.** Rank-contrast trees preserve within-sample ordinal invariance but are NOT yet demonstrated to add meaningful predictive benefit to the original RLR. Stop this mechanism; no threshold, lambda, pruning or EA tuning using inspected E folds. Original RLR remains the ready manuscript fallback. This experiment does not include actual original CUDA TSPDT benchmark.

These are source-reconstructed public matrices, not SHA256-confirmed equality with processed original paper ARFF files. Neither six repeated folds nor this negative pilot can establish external generalization or novelty. Literature on ordinal/rank-gap scores and model-based recursive partitioning remains an essential novelty control.
