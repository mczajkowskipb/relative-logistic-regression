# ORIGINAL (frozen) DeCAF versus original RLR-20 on PDAC CAF subtype labels — results

Date: 2026-10-30. **Full clinical-source model comparator, not a proxy**. Real GitHub Actions [#41095501198](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41095501198) ended **SUCCESS**. Independent Python patient-level audit (recomputed from exported out-of-training probabilities) also PASS.

## Sources and exact original predictor validation

- **Author's original frozen DeCAF Rdata and `decaf_functions.R`**, unchanged: https://github.com/jjyeh-unc/decaf.
- Author's original combined expression (18 official classifier genes), original `SCISSORS_CAF` labels, validation cohort membership and published predictions: https://github.com/jjyeh-unc/decaf_manuscript/tree/main/data/public_PDAC_clean/combined_data.
- Model's predictions exactly reproduce **1324/1324 author-published probabilities**, maximum absolute difference **7.77e-16**, mean absolute difference **1.93e-16**. Thus frozen model and assay gene matrix were matched to the reference.
- Original publication task: **proCAF vs restCAF**, not arbitrary CompCancer class labels. The `SCISSORS_CAF` labels are CC/SCISSORS subtype assignments, not independent pathology gold-standard diagnoses.

## Same clinical task and patient groups

Original model was pretrained by authors on CPTAC, Dijk, Moffitt_GEO_array, TCGA_PAAD studies. Train RLR20 from original authors' CAF subtype labels on those same four study groups, **372 eligible training specimens**. RLR uses ONLY 18 published DeCAF marker genes; same original RLR code `rg_parent`, fixed K=20 with strict d2 disjoint pairs, ternary rank gap m10, glmnet alpha=.5. The 18-gene panel constrains RLR to at most nine gene-disjoint pairs and is already selected to favor DeCAF's biology. No RLR tuning on independent validation labels.

Original combined author input included **52 repeated biological specimens between `PACA_AU_array` and `PACA_AU_seq`**. Explicitly removed duplicate seq copies (array preferred), pooled remaining `PACA_AU` into **one biological study unit**. External validation: **553 distinct specimens** in **6 independent validation study units**; no training vs external specimen ID overlap. All six units with both classes, no selective exclusion by performance. Per-study n and counts in artifact.

## Unweighted mean external validation study scores

| Model | Logloss ↓ | Brier ↓ | AUC ↑ | F1 ↑ |
|---|---:|---:|---:|---:|
| **Official original DeCAF** | **0.289858** | **0.088163** | **0.960308** | 0.897191 |
| Original RLR20 restricted to 18 DeCAF genes | 0.390312 | 0.101095 | 0.900144 | **0.900182** |

Original DeCAF beats RLR for mean logloss/Brier/AUC; RLR has a negligible mean F1 advantage. Two-sided exact signflips on six **independent study** units: mean difference RLR minus DeCAF logloss=**+0.100454** (RLR worse), **p=0.28125**, NOT significant. This is **exploratory** and does not justify a general DeCAF superiority claim either.

### Each validation study: logloss

| Study | n | RLR20 | Original DeCAF | Better |
|---|---:|---:|---:|---|
| Grunwald | 15 | 0.195906 | **0.146372** | DeCAF |
| Hayashi | 90 | **0.266559** | 0.323723 | RLR |
| Linehan | 47 | 0.415617 | **0.214945** | DeCAF |
| Olive | 114 | 0.600339 | **0.217109** | DeCAF |
| PACA_AU | 84 | 0.741756 | **0.530939** | DeCAF |
| Puleo_array | 203 | **0.121695** | 0.306059 | RLR |

PACA_AU counts after deduplication, all original samples/labels and relative model probabilities retained in artifacts.

## Two independent checks

1. Real author frozen DeCAF predictions reproduce published author predictions for **all 1324 input columns** (maximum absolute error <1e-15).
2. Separate external Python replayer recomputed logloss, Brier, AUC and positive-class F1 for **all 12 study/method summaries from all 553 held-out prediction rows**. Verified no duplicate validation sample IDs and no train/test specimen overlap. Reproduced six-study exact signflip p=0.28125.

Original R code and all outputs: [GitHub Actions #41095501198](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41095501198), artifact `original-decaf-rlr-headtohead-41095501198` (30-day retention). Files: `external_patient_predictions.csv`, `external_study_metrics.csv`, `dataset_equal_summary.csv`, `paired_exploratory_logloss.csv`, `original_decaf_reference_replay.csv`, `source_identity_audit.csv`.

## Crucial interpretation for CMPB

- Separate previously completed CompCancer matched-pair encoding ablation, on **11 original SHA-confirmed ARFF datasets**, favored RLR20 over two binary refits. Its statistical interpretation remains a **controlled representation ablation** on a different set of tasks.
- **It is INCORRECT to claim RLR beats original DeCAF**, given original clinical-task results above. The DeCAF-pretrained weights solve a different task and have a biological gene set specifically optimized for CAF subtyping, and these controls are not generic classifiers.
- The direct original DeCAF comparison is limited to its published 18-gene panel, different study-level training/validation sample composition and six units; the exploratory logloss result lacks statistical evidence. It should be treated as an important robustness counterexample, not as a confirmatory superiority comparison.
- Original PurIST has only frozen model and original inference code publicly at https://github.com/naimurashid/PurIST; a rigorous direct comparison needs independently defined *tumor* Moffitt basal-like/classical labels and a disjoint external cohort. **Do not use the `PurIST` prediction column in DeCAF subtype CSV as ground truth.**
- No changes to frozen RLR manuscript, main branch algorithms or source data. Do not transfer these results into paper as a general original-model victory. If reporting binary pair ablation in main, call it that and include the real DeCAF comparison as countervailing evidence when germane.

## Status

**ORIGINAL_DeCAF external same-task head-to-head: COMPLETE, AUDITED, NOT RLR-SUPERIOR.**
**ORIGINAL_PurIST same-task head-to-head: NOT YET SUPPORTED BY AN INDEPENDENT LABEL SET.**
