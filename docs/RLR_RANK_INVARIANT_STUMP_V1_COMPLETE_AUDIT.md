# RLR Rank-Invariant Stump V1 — protocol and audited outcome

Date: 2026-10-26. Branch: analysis/rlr-rank-invariant-stump-v1.
Actual run: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40420778168
Code SHA of fitted run: f514d1a1c4ba1b013f587a2fcdd2cf01a60a36af

## Purpose

Test a bounded, fully within-patient rank-based decision stump as an extension of the already completed but unpublished original RLR. No original RLR core or locked manuscript was changed. This is exploratory development on previously inspected original-paper cohorts, not external validation.

## Invariance and features

Only patient-wise ranks are used in new prediction tests; no raw expression values, ratios, cross-sample scale thresholds or raw-value MAD candidate filter. The two binary root variants are (i) ORDER: rank(B)-rank(A)>0; and (ii) WIDE: rank(B)-rank(A)>=max(10,ceil(0.15*(p-1))) where p is the fixed measured gene universe. The partner feature is the original signed ternary m10 RLR encoding. Candidate expansion uses rank MAD computed over the D training patients, not expression-value MAD.

Any strictly increasing transformation of all genes of one patient leaves the frozen-rank representation unchanged. This assumes fixed and identical gene/probe schema, no gene-specific transformations and no new ties. For one perturbed gene, an unrelated pair rank gap changes by at most one if there are no ties; a threshold-adjacent decision can still flip.

Mechanical GitHub R tests PASS: 120 independent patient-specific affine and nonlinear transforms; candidate expansion, pair gap and stump design exactly invariant; fixed gene schema enforced; outlier bound; synthetic XOR candidate formation. The synthetic result is a software test, not biomedical classifier validation.

## Algorithm

Parent is exact original RLR: switchBox top1000, original disjoint d2, original ternary margin m10, K20, elastic-net alpha .5, lambda.min. A single rank root gate g and a ternary rank partner z define two ridge-corrected offset comparisons:
GLOBAL_MAIN_CONTROL logit(p) = logit(p_RLR) + a*g + b*z.
RANK_STUMP logit(p) = logit(p_RLR) + a*g + b*z + c*g*z.
Both use glmnet ridge lambda .10, shared parent offset and no additional intercept; the interaction gets penalty factor 2.

Selection is inside outer-training T: class-stratified split into D and V. D computes parent RLR and generates at most eight residual-interaction candidates, satisfying minimum 12 examples and 3 per class in each branch; V is never used for feature nomination. The highest-V candidate is accepted only when V logloss beats original RLR by >=.005 and the same-root/same-partner no-tree control by >=.002, and V BA decreases by no more than .02. Otherwise the heldout E prediction is byte-exact original RLR. Accepted designs are refitted on full outer T before E prediction.

D gradient uses IN-SAMPLE RLR residuals (optimistic as a screening statistic) rather than cross-fitted residuals. V selection among eight proposals has its own selection variance. This should not be marketed as unbiased interaction discovery.

## Real-data experiment

Previously exposed original-paper cohorts: Chowdary-2006 (104 patients, 182 filtered gene features) and Singh-2002 (102, 339); source-reconstructed Schliep public matrices, NOT proven identical to the manuscript's processed ARFF hashes. Three fixed stratified outer folds per cohort, seed42, one repetition. Exactly 6 fold studies, 24 method-fold metrics, 824 heldout method-patient predictions.

Completed GitHub Actions PASS. All E patient logloss, Brier, BA, coverage, positive-class coding, fallback parity and manifest alignment recomputed in independent Python audit PASS. All purely rank-invariance mechanical unit tests PASS.

### Inner V acceptance and outer E

- Chowdary fold0: V gain +.00163, rejected; raw ungated E delta logloss +.000529.
- Chowdary fold1: no balanced eligible rank root; exact original RLR fallback.
- Chowdary fold2: V gain +.00174, rejected; raw ungated E delta - .010927.
- Singh fold0: V gain +.003589, rejected; raw ungated E delta - .007269.
- Singh fold1: V gain +.00150, rejected; raw ungated E delta - .001852.
- Singh fold2: V gain +.02575, accepted; raw ungated E delta +.000450.

Accepted truly contextual stumps: 1/6; that sole E evaluation got WORSE than original RLR. The other 5/6 gated fits returned exactly original RLR.

### Dataset/fold-unweighted mean E logloss differences (variant minus original; negative better)

| Method | Mean E delta | Improved folds |
|---|---:|---:|
| Original RLR | +0.000000 | 0/6 |
| Global same-pair main-effect control | -0.002151 | 5/6 |
| Stump ungated (diagnostic only) | -0.003178 | 3/6 |
| V-gated selected stump | +0.000075 | 0/6 |

Differences are very small and subject to model-selection optimism; no statistical claims are justified from six non-independent, repeatedly used folds. The global control had more favorable fold-count direction than the tree, implying that any occasional nominal benefit need not be tree-specific.

## Decision

**Scientifically negative for this specific stump mechanism.** It passes technical and rank-invariance tests, but there is no coherent E evidence of advantage from the learned tree itself. Do not tune the current gates, penalties, support cutoffs or lambda using these exposed E results and then label them fresh validation. If continuing the tree research, formulate a genuinely distinct feature-selection or model-induction mechanism and predeclare its controls. Keep completed original RLR as a publication fallback.

Actual original CUDA TSPDT has NOT been executed here; do not compare against the earlier greedy proxy as if it were TSPDT.
