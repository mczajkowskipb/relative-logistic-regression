# Original PurIST / DeCAF versus RLR — genuine clinical-task benchmark protocol

Date: 2026-10-30. Exploratory design fixed before running any biological-task comparison below. Existing RLR CompCancer comparison results known; **this is not preregistered confirmatory analysis**.

## Scientific question and required comparability

**PurIST (Rashid et al. 2020)** predicts basal-like versus classical PDAC *tumor* subtype. **DeCAF (Peng et al., Cell Reports Medicine 2026)** predicts proCAF versus restCAF *fibroblast* subtype. They solve DIFFERENT label tasks. Applying their frozen coefficients to the class labels of 11 unrelated CompCancer cohorts would be scientifically invalid. Their published AUCs are not directly comparable to the RLR CompCancer AUC.

All original-model claims require unmodified **authors' source and frozen Rdata classifiers**, complete published input normalization and independently defined target labels. Authors' repositories state non-profit research use only.

Original sources:
- PurIST: https://github.com/naimurashid/PurIST (trained `fitteds_public_2019-02-12.Rdata`, inference code `functions.R`; public repository contains NO full generic learning pipeline)
- DeCAF: https://github.com/jjyeh-unc/decaf (trained `decaf_classifier.Rdata`, inference code `decaf_functions.R`)
- DeCAF full original training scripts and public study data: https://github.com/jjyeh-unc/decaf_manuscript

## Phase 0: immutable model identity (implemented)

Reproduce authors' prediction calls on authors' original example datasets with the exact official pretrained objects and code. This is an integrity test only, not head-to-head performance.

## Phase 1: ORIGINAL DeCAF / RLR same-outcome external comparison (implemented pilot)

Use authors' `combined_classifier_gene_exp_DeCAF.csv` (published 18 classifier gene measurements) and `combined_subtype_DeCAF.csv`, plus `combined_classifier_calls_DeCAF.csv` as an optional version check. The **true outcome** is `SCISSORS_CAF` (unsupervised CAF subtype, proCAF/restCAF), **NOT** the DeCAF prediction column. Use original authors' four training study names: CPTAC, Dijk, Moffitt_GEO_array, TCGA_PAAD. Use remaining independent named studies as external validation units; require both outcome classes in each reported study. Training labels are limited to `CAF_whitelist=TRUE` and valid `SCISSORS_CAF`, without accessing validation labels during model fit.

DeCAF is the **ORIGINAL FROZEN** pretrained author model, with unmodified `apply_decaf`. RLR is the original frozen RLR20 procedure, but trained on **the same 18 predefined DeCAF classifier genes only**, not on the wider array expression universe. Its K20,d=2,m=10,alpha=.5 constraints are held fixed from the established RLR algorithm. Consequently this experiment is an exploratory comparison restricted to the author's selected 18-gene panel, **not evidence that RLR beats DeCAF across general cancer gene-expression tasks**. In particular, this favors DeCAF feature selection and constrains RLR disjoint pairs to <=9.

Endpoint primary descriptive mean validation-study logloss, secondary Brier/AUC/F1, study-equal aggregation (not treating patients across studies as iid). Compute exact paired two-sided sign-flips over independent validation studies only when enough eligible study units. Never tune parameters after observing external labels. Record all cohorts including ineligible/missing-class and version-discrepant ones; do not select favorable cohorts.

Potential caveat: source-authored `combined_classifier_gene_exp_DeCAF.csv` may be a transformed subset, check official predicted probabilities against published `combined_classifier_calls_DeCAF.csv`; if they disagree, report this and investigate before any manuscript claim. Study-labeled CC subtype is not necessarily an independent clinical gold standard and may reflect the source pipeline.

## Phase 2: ORIGINAL PurIST / RLR same-outcome comparison (blocked pending independently defined labels)

Use published, author-independent **Moffitt tumor-intrinsic basal-like/classical** labels on PDAC cohorts not used for training PurIST. DO NOT use `PurIST` column of DeCAF datasets as ground truth: that is the model's own predicted label and would make the original model tautologically perfect. RLR must train from suitable training labels (prefer PurIST's original clinical training cohorts) and be evaluated on the SAME samples with the SAME target and source processing as frozen PurIST. Historical sample overlap with either original model training cohorts must be documented and excluded from claimed external validation. Need recover consistent original independent Moffitt CC labels or original patient label tables before proceeding.

## Publication decision

Original-ARFF CompCancer ablation result (RLR logloss 0.2379 versus binary ENET 0.2904; BH q=.0391) is valid as exploratory **binary rank-feature ablation**, NOT original-PurIST superiority. No claim of direct original DeCAF/PurIST superiority in MAIN/Abstract until genuine matched-label external comparisons, exact source-method fidelity checks, all cohorts and adverse results. If confirmed, legitimate main claim distinguishes methods and phenotypes precisely.
