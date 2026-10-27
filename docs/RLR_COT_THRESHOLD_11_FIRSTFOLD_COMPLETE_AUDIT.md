# Conditional Ordinal Threshold Decision Stump for unpublished RLR — COMPLETE first-fold 11-cohort experiment

**Date:** 2026-10-27. **Technical run full PASS, scientific tree gate NEGATIVE.** The original RLR is NOT yet published, and its frozen article and actual implementation remain unmodified. The target is a single combined RLR + genuinely useful interpretable extension paper against actual k-TSP.

## Experimental definition

The continuous probability model is the exact original RLR50 (K50, d2, m10, original `cv.glmnet` alpha=.5 and `lambda.min`). Preserve all continuous scores and probabilities. Search only for a *new classification rule*, not a new probability model.

1. On original outer T, stratified `rg_inner` split D/V. On D, 3-fold OOF refit original RLR50 and obtain OOF p for **train-only threshold choice**.
2. Fit one global probability decision threshold `t` on D from a fixed grid .10..90 by .05, optimize original article's positive-class F1 (both Macro-F1 and BA also reported).
3. Generate rank-only, original-paper-gene-universe anchor roots `g > t-of-3/5 anchors` with at most one anchor exchange and 240 unique patient partitions. For each find two conditional *thresholds on p_RLR*, optimizing OOF D positive-F1 relative to single D-global threshold, with BA guard and minimum branch support. Retain at most 3 D finalists.
4. Evaluate tree against **both** original RLR threshold .5 and the global D-OOF threshold on V. Accept only if tree improves V positive-F1 by >= .02 and balanced accuracy relative to globally thresholded RLR decreases no more than .02. Otherwise use global threshold. All roots, anchors, thresholds and decisions are chosen without outer E outcomes.
5. Test the frozen rule once on original outer E; report the same original RLR p in all RLR threshold arms. Samplewise strict monotone transformations leave both root and original parent relation features unchanged with fixed gene universe.

## Real 11-cohort first original fold experiment — final auditable evidence

- GitHub run [#40491250870](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40491250870) **ALL JOBS SUCCESS**, code SHA `991db0c138dcab31d4cb05b79aa7a6bce312f9e7`.
- Public-source reconstructed matrices for the SAME original 11 CompCancer paper cohort names, original stratified 5-fold seed42, **first external E fold only per cohort**; source SHA manifests retained, original processed **ARFF checksum parity NOT established**.
- Real `switchBox::SWAP.KTSP` K50 (not a proxy), original K50 RLR, original RLR .5, global threshold, ungated conditional threshold stump, inner-V-gated conditional stump.
- 11 × 5 = **55 patient-matched method-fold records**, independently recomputed from archived heldout patient probabilities, labels and split manifests in Python. AUC, positive-F1, Macro-F1, BA, Brier, logloss, exact fallback and E ID/truth parity PASS.
- **Zero of eleven trees accepted in inner V**.

| Arm | Mean positive-class F1 over 11 datasets |
|---|---:|
| real k-TSP50 native label | .865621 |
| original RLR50 at .5 | **.867924** |
| original RLR50 with D-OOF global threshold | .859316 |
| ordinal conditional threshold stump without V gate (diagnostic only) | .874773 |
| ordinal conditional threshold stump WITH mandatory inner-V gate | **.859316** |

Note: Ungated conditional stump has numerically highest mean F1 but was rejected on all 11 V splits; **it is not an accepted model or evidence of reproducible gain**. It is an overfitting diagnostic, not a positive comparative result. The global threshold **worsened** average F1 compared with original RLR .5. The complete score vector (and thus AUC, logloss and Brier) is **bit-identical** among all four RLR threshold variants, as independently asserted.

### Fold-specific real E positive-class F1

| Cohort | E patients | Real kTSP50 | Original RLR50 .5 | D-OOF global t | Inner-V gated conditional threshold |
|---|---:|---:|---:|---:|---:|
| Alizadeh | 10 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| Armstrong | 15 | .9474 | 1.0000 | 1.0000 | 1.0000 |
| Bittner | 8 | .5714 | .7500 | .7500 | .7500 |
| Chen | 36 | .8000 | .9032 | .9375 | .9375 |
| Golub | 15 | .9091 | .9091 | .9091 | .9091 |
| Gordon | 37 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| Pomeroy | 7 | 1.0000 | .6667 | .6667 | .6667 |
| Shipp | 16 | .6667 | .7500 | .7500 | .7500 |
| Singh | 21 | .8182 | .9091 | .9000 | .9000 |
| West | 10 | .9091 | .9091 | .8333 | .8333 |
| Yeoh | 50 | .9000 | .7500 | .7059 | .7059 |

Only the Chen global-threshold variant improved (F1 from .9032 to .9375 and BA from .9190 to .9524); in West, Singh and Yeoh the global threshold lowered positive F1. All selected conditional trees were rejected.

### Mechanism diagnostic: actual error complementarity in SAME 225 E patients

Independent local post-analysis of the 11 archived patient-level artifact CSVs (purely descriptive E-label analysis, **no fitted new classifier**):
- total E patients **225**
- both classifiers correct **197**
- both classifiers wrong **10**
- original RLR50 correct but kTSP50 wrong **12**
- real kTSP50 correct but RLR50 wrong **6**
- both classifiers disagree on **18/225** predictions.
- A hypothetical oracle that picks whichever of the two is correct would fix up to six of RLR's 16 errors, reaching 215/225 accuracy rather than original 209/225; **this is not a feasible learned model and provides only a dataset-specific upper bound**.
- kTSP-only correct cases occur mainly in Yeoh (4) and Pomeroy (1), with one Chen; RLR-only correct are spread across Chen, Singh, Yeoh etc. Hence no robust general-purpose evidence of a large error-complementarity opportunity in this sample.
- **These are repeated E results**, not to be used to select, tune or validate a new hybrid algorithm on the same patient folds.

### Pre-E-run algorithm clarification

The run above uses the `ct_global` training F1 ranking without an explicit BA relative-to-0.5 constraint. During engineering, a further **train-only safeguard** was subsequently added on this branch to require global-t BA at least original .5 BA minus .02; this is a **new code revision** and was NOT part of fitted run SHA `991db0c...`. Do **not** present the table above as results of the revised safeguard code. All current code differences are explicit. Any rerun on these repeatedly seen E samples remains exploratory, not a new confirmatory 11-cohort test. The same branch includes synthetic invariance tests; do not assert their latest run succeeded until corresponding workflow status is PASS.

## Research conclusion

1. Stop the stand-alone conditional-threshold stump at this evidence stage: zero accepted trees in 11 internal V checks, no general F1 advantage over a plain original RLR50 decision rule.
2. A decision threshold is routine postprocessing, not a substantive new interpretable rank-tree mechanism. It cannot independently justify a tree-method publication.
3. **Prioritize a newly constructed conditional RELATION SELECTION criterion** rather than repeatedly searching gate thresholds: current child models select pairs by class-separation switchBox even though the original RLR50 already explains that main signal. Select candidate child pairs by **cross-fitted residual partial score** `U_j = sum_{i in child} z_ij(y_i - p_RLR50,OOF,i)`, with Hessian shrinkage `H_j = sum z² p(1-p)`, then fit an actual sparse correction with original RLR50 as offset. Compare against a no-tree model that uses the identical feature union and coefficients rather than attributing plain additional features to routing. This *has not yet been validated*.
4. Fixed disjointness `d=2` original baseline is preserved and exact original paper source kept untouched. Train-only conditional pair shortlist may relax shared genes between **different branches** while penalizing overlap along each patient path. This is more aligned with the author's original anchor-gene hub hypothesis.
5. The actual final paper must use **all original 11 × 470 repeated outer folds** on exact processed ARFF if available and report dataset-level paired differences, RLR K10/20/50, real switchBox kTSP K10/20/50, value-LR, kTSP calibrated Platt, and the selected compact tree versus same-feature global/no-tree controls. With no verified tree gain, submit original completed UNPUBLISHED RLR paper rather than exaggerate a new contribution.
