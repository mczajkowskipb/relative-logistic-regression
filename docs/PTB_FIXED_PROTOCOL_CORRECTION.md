# PTB fixed-split protocol and ORA correction note

Date: 2026-09-10

## Status

**CLOSED.** No RLR, k-TSP, CompCancer, PTB LOCO, or legacy8 predictive model was rerun. The canonical PTB fixed reference state remains unchanged:

`release_reference/ptb_fixed_published_state`

## Fixed-split protocol correction

The audit established that the canonical PTB fixed results are **not** based on GSE83456 as an external test cohort.

- `01_train_data.arff`: 103 samples (69/34)
- `02_test_data.arff`: 44 samples (24/20)
- `03_valid_data.arff`: 106 samples (61/45)
- canonical `predictions_test.csv`: 44 rows for all 9 model IDs
- 103 + 44 = 147, matching the pooled discovery cohorts GSE19439 + GSE19442 + GSE19444
- GSE83456 corresponds to the 106-sample external-validation cohort in the source benchmark, but it is not the canonical PTB fixed test set used for the reported RLR/k-TSP metrics

The manuscript wording is therefore corrected from an “external cohort split” to a **predefined discovery-data fixed split (103 training / 44 test)**. PTB LOCO remains the separate descriptive cross-cohort stress test over all four cohorts.

## ORA correction

The previous ORA background was invalid for a training-only claim because `01_train_valid_data.csv` contained 209 samples. The audit found only 874 shared genes between its MAD-2000 universe and the true training-only MAD-2000 universe (Jaccard-style overlap 0.2796).

The corrected ORA uses only:

`data/cohort_norm/01_train_data.arff`

and reconstructs the training-only MAD top-2000 universe.

Canonical selected/active pair files are read from the frozen PTB fixed reference state; no predictive models were refitted.

## Gene-set source

Because `msigdbr 26.1.1` could not download its first-use data from Zenodo (HTTP 504), the analysis uses official MSigDB 2026.1.Hs symbol GMT files downloaded from Broad and stored under:

`reference_outputs/ptb_ora/gmt/`

SHA256:

- Hallmark: `eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596`
- GO Biological Process: `9be09dd06d6652566eb52eed530d62e6dfecc4365c1e81afd6f0b7f2e86dd4f9`
- Reactome: `5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341`

Corrected ORA output:

`reference_outputs/ptb_ora`

Final run reached `DONE`.

## Corrected ORA summary

Training-only universe: 2000 genes. All extracted case-study genes are present in that universe:

- k-TSP-10: 20/20
- RLR-10: 16/16
- k-TSP-20: 40/40
- RLR-20: 26/26

FDR-significant enrichment rows:

- k-TSP-10: 12 (2 GO BP, 3 Hallmark, 7 Reactome)
- k-TSP-20: 20 (8 GO BP, 3 Hallmark, 9 Reactome)
- RLR-10: 5 (3 Hallmark, 2 Reactome)
- RLR-20: 3 (3 Hallmark)
- total: 40

The corrected interpretation is deliberately conservative: all four gene lists retain an immune/host-response theme with recurrent interferon-related enrichment, while the smaller RLR active sets yield fewer significant terms. The difference in enrichment breadth is **not** interpreted as evidence of superior biological specificity because the model-derived gene lists differ in size.

## Manuscript files produced

- `RLR3_locked_v3_ptb_corrected.tex`
- `Supplement_RLR_locked_v3_ptb_corrected.tex`

These supersede the PTB-specific wording and ORA content in the v2 locked files. The v2 files should not be used for final submission because they still contain the incorrect fixed-split interpretation and obsolete ORA values.
