# CRIT-v9 D-internal candidate ranking, single V confirmation: locked development protocol

Frozen 2026-10-25. **Experimental, non-publication-ready.** Original RLR, v7, v8 and all prior benchmark results are immutable and remain preserved.

## Basis: experimentally verified stages v7–D

- v7 (Chen/Gordon 2×3): accepted two Chen trees, mean 6-fold Δlogloss −0.001923086. Acceptance used the maximum validation score over six discovery-D-ranked candidates, i.e., *selection was performed on V*.
- v8: nominated only discovery candidate #1, conducted one predeclared V test per correction type. All 30 V-arm tests rejected; all 36 test prediction vectors exactly equaled the published original RLR.
- Stage C: archived original data, D/V patient partitions and v7 first D candidate's gate, leaf pairs, V mean, SE and lower score exactly matched v8. V-max chose discovery ranks (3,3,6,1,4,5); the two v7-accepted Chen candidates were ranks 3 and 6.
- Stage D: 36 actual V-patient candidate gain vectors reconstructed and checked against frozen v7 scores, with 30,000 paired artificial-zero-mean simulations and independent audit. The **Chen-fold D-first/V-max6 artificial-null acceptance fractions** were 6.34%/19.52%, 12.38%/41.24%, and 10.10%/36.00%. The approximate Bonferroni-style rule gave 0%, 1.58%, 0.80%. All Gordon folds gave 0% for the three rules. The average increase in these artificial-null acceptances across Chen was about 3.1–3.6 fold. **These are not calibrated type-I error rates or p-values for CRIT**, because the sign-flipped null reuses exposed patient matrices with uncertain exchangeability.

## New CRIT-v9 methodology

**Primary objective:** combine the broader six-candidate proposal opportunity with a single final confirmation decision, avoiding selecting the maximum heldout V score. The proposal and coefficient-ranking algorithm is implemented under `R/experimental/crit_v9_d_crossfit_one_shot.R`, reusing the exact frozen original RLR and CRIT-v7 function primitives.

For each **outer training fold T**:

1. Fit the exact original RLR parent (same switchBox original top1000 pair screening, K=10, m=10, d=2, glmnet α=.5, lambda.min), and get genuine original RLR 3-fold OOF predictions for T.
2. Partition T into **discovery D and confirmation V** via fixed seed+731; both disjoint, test E remains untouched. Fit the exact RLR parent and genuine parent OOF inside D.
3. On **D only**, construct the v7 candidate pool and shortlist up to six gate/leaf-pair candidates in fixed discovery order. Nomination uses D outcomes and D OOF residuals; NO V outcomes influence these features.
4. Rank these at-most-six candidates by *within-D coefficient-crossfit* heldout logloss gains. Use one fixed stratified 3-fold D partition (seed+1907). For each candidate, fit its four ridge correction coefficients on D\fold against D-OOF parent logit offsets, apply on Dfold with its D-OOF offsets, aggregate the D gain mean and standard error, and rank by `mean − 1×SE − 0.002`. Tie break by original discovery order. **Only the highest D crossfit-ranked nominee proceeds to V**.
5. IMPORTANT scientific limitation: the nominated feature pairs were already selected using all of D's labels; and D-wide parent-OOF probabilities are reused in the internal correction crossfit. Therefore these within-D numbers are **not an unbiased nested CV estimate**. They are only an entirely **D-only deterministic candidate-selection heuristic**. The genuinely heldout V nevertheless has no role in candidate nomination.
6. For the single chosen candidate, use five predeclared correction designs (GLOBAL_CAL, GATE_ONLY, ADDITIVE_PAIRS, GLOBAL_INTERACTION, CONDITIONAL_TREE). Each is trained on D and tested **once** on V, independently accepted if `mean_gain − 1×SE − 0.002 > 0.002` (same frozen ridge/threshold/cap as v7-v8). No choosing of the best arm on V is allowed for a confirmatory claim. **CONDITIONAL_TREE is the primary hypothesis; other arms are secondary controls.** There is no multiplicity correction across the five arms, and the 1-SE rule is heuristic.
7. For an accepted arm, refit fixed correction coefficients on ALL T with the exact original RLR OOF T offsets; for rejected arms, return **bit-exact original RLR predictions** for E. Score independent outer test E only after all T-only decisions; the original RLR implementation and archive stay untouched.

## Frozen experiments, under one guarded batch

### A. Synthetic structural tests

- 32 genes ensure original RLR m=10 does not create all-zero predictors.
- Exact original-RLR identity, six arms, D/V disjointness, at-most-six D candidates, max D crossfit rank score, five separately evaluated V models and precise original fallback.
- Full rank invariance under positive affine patientwise gene transforms, RDS bit replay, exact GLM interaction equivalence for accepted tree and repeat-run determinism.
- Forced zero-candidate branch must revert all relational corrections to the original model.

### B. Fresh synthetic development panel

**24 independent train-test studies** (4 predeclared mechanisms × 6 independently seeded replicates, 180 training and 120 test patients each):

- Null: `P(Y=1)=0.5`, independent of all gene ranks.
- Global pair: `logit P(Y=1) = 1.6 r` for planted pair `r=sgn(gene03-gene04)`.
- Gate-only: `logit P(Y=1) = 1.6 z` for planted gate `z=sgn(gene01-gene02)`.
- XOR rank interaction: `logit P(Y=1) = 2.4 z r`.

All features are continuous Gaussian, labels Bernoulli conditional on the explicit mechanism. Seeds are fixed in source before the experiment. All six models produce heldout logloss and Brier; store **17,280 synthetic test patient predictions**. Independent Python independently recomputes metrics, class orientation, exact fallback, alignment, 24×6 model inventory and D ranking provenance. Only six null replications means **this panel cannot determine type-I error at 5%**. No scientific superiority claim from this panel; results are engineering/development feasibility.

### C. Exposed biomedical original-paper parity panel

Exactly Chen-2002 and Gordon-2002 × first 3 published original CV folds × 6 methods = 36 patient-level vectors. Verify original input SHA, individual train/test patient IDs and labels, exact original RLR patient probabilities, frozen 3fold RLR OOF, D-only candidate rank provenance, separate V acceptance flags, original fallback and Python-independent logloss/Brier. **These datasets, methods and splits are exposed and cannot serve as independent validation.** All original RLR, v7, v8 and Stage C/D result files remain unchanged.

## Completion criteria

- `CRIT9_SYNTHETIC_ALL_GATES=PASS`
- `CRIT9_24_FRESH_SYNTHETIC_DATASETS_SIX_ARMS=PASS`
- `CRIT9_17280_SYNTHETIC_PATIENT_PREDICTIONS_AUDITED=PASS`
- `CRIT9_36_METHOD_FOLDS=PASS`
- `CRIT9_EXACT_ORIGINAL_RLR_PARITY=PASS`
- `CRIT9_D_CROSSFIT_36_VECTORS_AUDITED=PASS`
- `CRIT9_BATCH_EXIT=0`.

These are *planned gates*, not achieved results until the Ubuntu run completes. No threshold changes based on Stage D or v9 pilot results.

## Scientific decision after this stage

Use the independent synthetic null and planted-alternative test sets to look for genuine conditional gain over ORIGINAL_RLR and GATE_ONLY, not mere improvements on exposed folds. Inspect acceptance/rejection and loss distributions separately. If v9 continues to reject all valid interactions under planted XOR, revise the theoretical nomination objective **on new development simulations only**, not using the original Chen/Gordon fold E labels. If v9 exhibits reproducible null acceptances, formalize a valid error-controlled confirmation procedure with future heldout datasets and larger simulation samples. Any eventual primary biomedical efficacy claim requires data not used during design or debugging (ideally patient-independent cohort and external batch).
