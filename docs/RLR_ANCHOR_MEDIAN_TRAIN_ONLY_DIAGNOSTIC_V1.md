# Anchor majority — training-only median identity and single-pair collapse diagnostic

**Date:** 2026-10-27. **Status: completed real R run SUCCESS.**
[GitHub Actions #40423626886](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40423626886), fitted SHA `615dc91ab886fce626a7333c4806cd44180aedb7`.

This is a new read-only **TRAINING-ONLY diagnostic**, not a revised classifier, not a retuned method, and it never inspects new heldout outer E labels or probabilities. Input cohorts (Chowdary-2006, Singh-2002) are already-exposed public reconstructed original-paper source matrices, not SHA-verified original processed ARFF.

## The mathematical concern

The reference vote `g > majority(a1,...,ak)` equals `g > median(anchors)` for odd k. If the identity of the median reference never changes among patients, the group is **exactly one simple gene-pair classifier**, no matter how many anchors appear. When the median identity changes, the group can differ from a single fixed pair, but may remain practically redundant.

The diagnostic chooses an anchor root separately on each outer **TRAINING partition** only, using the frozen algorithm code from the parent Anchor Majority Tree v1. It records how often the central (median) anchor identity differs, how many group votes match a *fixed* most-frequently-median anchor, and how much the hub target's normalized within-sample rank shifts between classes relative to anchors.

## Verified actual run on six outer training partitions

| Cohort | Fold | Anchors | Distinct median identities | Most common median identity | Majority agreement with single modal anchor | Target normalized rank shift | Median reference shift | Anchors with shift < half target |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Chowdary | 0 | 3 | 3 | 39.1% | **97.1%** | .6011 | .0293 | 3/3 |
| Chowdary | 1 | 5 | 5 | 26.1% | **95.7%** | .5729 | .1283 | 5/5 |
| Chowdary | 2 | 3 | 3 | 35.7% | **97.1%** | .6010 | .0347 | 3/3 |
| Singh | 0 | 3 | 3 | 41.8% | **89.6%** | .2119 | .1982 | 0/3 |
| Singh | 1 | 3 | 3 | 44.1% | **97.1%** | .4163 | .0321 | 3/3 |
| Singh | 2 | 3 | 3 | 40.6% | **94.2%** | .2017 | .1520 | 0/3 |

Unweighted six-fold average training **agreement of full panel vote with single fixed modal-reference gene: 95.12%**, despite an average of **3.33 distinct genes taking the median identity** across patients. The reference panel changes which anchor is its actual median, yet usually returns the same target comparison as an ordinary fixed pair.

**Interpretation:**
1. There is direct training-only rank evidence of a hub with substantially larger class-conditional normalized rank shift than references in three Chowdary folds and one Singh fold, supporting the intuition for at least some candidate signatures. But rank shifts do not demonstrate raw-expression shifts or the clinical molecular mechanism.
2. In two Singh folds, none of the chosen references had class-rank drift < half the hub's drift. Calling them stable reference anchors would be misleading.
3. 95.1% decision agreement with one reference suggests severe *functional redundancy*, even though median-gene identity itself is not constant. This explains why naive majority voting can fail to add meaningful predictive capacity compared with TSP or original RLR.
4. An alternative theoretical method could require **reference reliability across patient perturbations, majority-vs-best-single-pair incremental gain, and truly diverse conditional responses**, assessed within D, but must not be tuned to already-seen E outcomes and later needs untouched independent validation. No such method has yet been shown to win.
5. Nothing in this diagnostic demonstrates generalization advantage or novelty of the anchor group. The previous completed original RLR is still the manuscript fallback.

Machine log assertions: `ANCHOR_MEDIAN_6_ORIGINAL_PAPER_TRAIN_FOLDS_DIAGNOSED=PASS` and `ANCHOR_MEDIAN_ALL_NEW_E_PATIENT_LABELS_UNREAD=PASS`.
