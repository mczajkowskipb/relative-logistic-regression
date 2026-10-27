# FINAL: SAVS selective four one-sided rank votes — 11 CompCancer × 470 outer folds
27 October 2026. The original RLR article is **UNPUBLISHED and LOCKED**; no changes to the original RLR manuscript or implementation.

## Immutable, independently audited evidence
- [All 470 original-design outer folds, 11 methods each, run #40569991819](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40569991819): **95/95 jobs SUCCESS** (94 repeated 5-fold blocks and independent audit). Exact test-patient assignments per repetition, label/class correspondence, each probability and 5,170 method-fold F1/BA/AUC/LL/Brier predictions independently replayed. No test labels used for training or feature-selection hyperparameters.
- [Independent EXACT 2^11 paired 11-dataset statistical inference, run #40583400575](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40583400575): **PASS**; seeded 30,000 bootstrap dataset-level resamples and exact two-sided sign flip enumerating all 2,048 signs for each paired comparison. Original source script [rlr_savs470_final_inference_v1.py](../scripts/rlr_savs470_final_inference_v1.py). Five output files archived in GitHub artifact `RLR-selective4-470-evidence-inference-40583400575`.
- [First outer fold all11 #40569478318](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40569478318) 12/12 PASS independently, but DO NOT generalize firstfold performance, which overestimated four-rule advantage.
- **Input caveat:** paper-named public Schliep CompCancer matrices reconstructed from public source files with manifests, but the original author's processed ARFF byte/feature identity has **NOT** been established. Prior V/E results of these cohorts were seen in research. This is a **reproducible development result**, not untouched confirmatory evidence or external clinical validation.
- Original design 9 cohorts × 10 repetitions × five folds + West/Yeoh 2×5 folds = **470** folds per method; each of **11 datasets is equally weighted**, never 470 pseudo-independent tests.

## Main dataset-equal-weighted means
| Method | Positive F1 ↑ | Macro-F1 ↑ | AUC ↑ | Logloss ↓ | Mean active relations | Mean pair comparisons evaluated |
|---|---:|---:|---:|---:|---:|---:|
| TRUE kTSP4 | 0.814501 | 0.854853 | 0.935723 | 0.727450 | 4.000 | 4.000 |
| TRUE kTSP10 | 0.860293 | 0.887541 | 0.954313 | 0.378239 | 10.000 | 10.000 |
| TRUE kTSP20 | 0.871680 | 0.894456 | 0.958268 | 0.346522 | 20.000 | 20.000 |
| **TRUE kTSP50** | **0.883625** | **0.900266** | **0.962416** | 0.377593 | 49.273 | 49.273 |
| **ORIGINAL RLR50** | **0.868613** | 0.897200 | 0.957076 | **0.226026** | **25.548** | 28.093 |
| Original value LR | 0.829073 | 0.869138 | 0.953675 | 0.269454 | — | — |
| **NEW SELECTIVE4 m10** | **0.812587** | 0.844936 | 0.901404 | 0.315791 | **1.879** | **4.000** |
| NEW SELECTIVE4 m10pct | 0.797536 | 0.835267 | 0.892204 | 0.328864 | 1.788 | 4.000 |

- The original RLR's reported **25.55 active patient-specific signed ternary rank pairs** means ±1 contributions, with inactive if rank gap is too small (feature 0) or coefficient exactly 0. New SAVS `h=1[sign*rankgap>=m]` has a **different one-sided notion of active**: absent rule vote even if the opposite large gap is present. Do not conflate these percentages or imply biology used fewer assayed genes. RLR50 original global active mean is 28.093, vs nominal K=50; kTSP50 actual realized mean is 49.273.
- **NEW selective4 has 4 tested, ~1.88 firing relations** per patient; true kTSP4 also evaluates 4. The new model's smaller *active explanation*, not test count, differentiates it from actual kTSP4.

## Inference properly paired by independent datasets (NOT CV folds)
Positive delta F1/AUC helps model A; negative delta logloss helps model A. All p two-sided unadjusted exploratory; percentile 95% CIs from 30k 11-dataset bootstraps.

| SELECTIVE4 m10 minus | Δ positive F1 [CI], p | Δ AUC [CI], p | Δ logloss [CI], p |
|---|---|---|---|
| **TRUE kTSP4** | **−0.0019** [−0.0221, +0.0166], **p=.8594** | **−0.0343** [−0.0518, −0.0183], **p=.0029** | −0.4117 [−0.7677,−0.1552], p=.0029 |
| TRUE kTSP50 | −0.0710 [−0.1278,−0.0193], p=.0342 | −0.0610 [−0.0900,−0.0331], p=.0020 | −0.0618 [−0.1153,−0.0093], p=.0576 |
| ORIGINAL RLR50 | −0.0560 [−0.0965,−0.0213], p=.0098 | −0.0557 [−0.0787,−0.0326], p=.0010 | +0.0898 [+0.0449,+0.1308], p=.0049 |

**No equivalence or non-inferiority was proved** by near-tied F1 vs kTSP4: a priori equivalence margin and an independent validation dataset were NOT employed. New selective4 AUC is appreciably inferior to actual kTSP4; F1/AUC are inferior to original kTSP50 and RLR50. Such findings are evidence AGAINST treating sparse abstention voting as better main classifier.

## Heterogeneity and failure points
Per-dataset positive-F1 SELECTIVE4 m10 versus kTSP4 (means across each dataset's folds):
- Alizadeh .810 versus .813; Armstrong .985 versus .981; Bittner .798 versus .756; Chen .863 versus .867; Golub .887 versus .925; Gordon .968 versus .951; Pomeroy .475 versus .452; Shipp .750 versus .795; Singh .879 versus .865; West .869 versus .832; Yeoh .655 versus .722.
- kTSP50 particularly dominates in Alizadeh (.962 vs .810), Pomeroy (.712 vs .475), and Yeoh (.890 vs .655).
- m10pct looked good in some firstfold inner D/V cohorts, but the full external repeated splits yielded lower F1 .798 and lower AUC .892 than original m10; a warning not to pick complexity thresholds on already-explored cohorts.
- Original RLR50 retains excellent LL .226 and strong F1 .869; it is itself an important unpublished contribution and should not be delayed to manufacture a positive tree.

## Tree ambition: fail honestly
- Root-first SAVS1 tree (inner V all11) F1 .7669 vs global 4 vote .8731 at m10pct; honest train split SAVS3 and root-child shared-gene SAVS4 failed to improve upon the global control.
- SAVS4 hub sharing is mechanically sound and strictly within-patient rank invariant ([standalone CI #40568133585](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40568133585) PASS), but shared hub **tree** m10pct F1 .7742, honest nonhub tree .7848, global four .8731 (inner V 11).
- Thus the desired **new high-accuracy decision-tree extension of original RLR is NOT achieved**. The global four-rule system is an interpretability ablation or exploratory sparse model, not a "tree".

## Research-integrity and concrete next work
1. **Preserve and prepare to submit original UNPUBLISHED RLR article.** Keep a distinct idea-development branch for tree research; do not modify RLR source/manuscript.
2. Stop E-driven variant search on the same 11 heavily explored CompCancer reconstructions; any new tree performance claim requires SHA-matched author's original processed ARFF (historical parity) plus truly untouched external omics cohorts (new confirmation).
3. Specific next tree hypothesis should directly address kTSP50's better rank discrimination: keep a **strong global RLR base** and use a rare branch correction only when OOF evidence supports a new class boundary. Alternative class-imbalance-adjusted *training-only decision threshold* should be benchmarked against original RLR separately, not sold as a tree novelty. Both require predeclared selection and untouched final tests.
4. Original publication other worlds PTB fixed/LOCO and Legacy8 are not completed in current GitHub environment because the original processed ARFF matrices are absent; do not say "all original manuscript data" tested.
