# Data acquisition and preprocessing

This repository does not silently redistribute processed biomedical datasets. The expected processed inputs and their SHA256 hashes are listed in `checksums/DATA_SHA256_EXPECTED.txt`.

## CompCancer

The 11 primary datasets are the filtered CompCancer benchmark matrices distributed by the Schliep Lab supplementary site for de Souto et al., *Clustering cancer gene expression data: a comparative study* (BMC Bioinformatics 2008, 9:497). Source landing page:

`https://schlieplab.org/Static/Supplements/CompCancer/`

Dataset inventory page:

`https://schlieplab.org/Static/Supplements/CompCancer/datasets.htm`

The publication benchmark uses:
`alizadeh-2000-v1`, `armstrong-2002-v1`, `bittner-2000`, `chen-2002`, `golub-1999-v1`, `gordon-2002`, `pomeroy-2002-v1`, `shipp-2002-v1`, `singh-2002`, `west-2001`, `yeoh-2002-v1`.

Final schedule: the first nine datasets use 10 repetitions x 5 folds; West-2001 and Yeoh-2002-v1 use 2 repetitions x 5 folds, for 470 outer test folds per model. The two execution shards are directories `data/CompCancer_truth_9x10/` and `data/CompCancer_truth_2x2/`; `data/CompCancer/` is the combined 11-dataset directory used by reporting/supplementary scripts.

The exact ARFF files used in the locked analysis are identified by SHA256 in the checksum manifest. Conversion from the Schliep text matrices to these ARFF files should be verified by checksum before a public reproduction is declared exact.

## Pulmonary tuberculosis (PTB)

Source benchmark: Wu et al., *Applications of gene pair methods in clinical research: advancing precision medicine*, Molecular Biomedicine 6, 22 (2025), DOI `10.1186/s43556-025-00263-w`; accompanying code DOI `10.5281/zenodo.14948408`.

GEO cohorts: GSE19439, GSE19442, GSE19444, GSE83456. The source study reports limma 3.54.0 background correction/normalization, probe-to-gene-symbol conversion, intersection of genes across cohorts, a 7:3 discovery split of the first three cohorts, and GSE83456 as an independent external validation cohort.

For the canonical RLR PTB fixed case study in this paper:
- `data/cohort_norm/01_train_data.arff`: 103 discovery samples;
- `data/cohort_norm/02_test_data.arff`: 44 discovery samples;
- `data/cohort_norm/03_valid_data.arff`: 106 external-validation samples (GSE83456) and **not** the test set used for canonical PTB fixed metrics.

PTB LOCO uses one whole cohort as holdout at a time from `data/cohort_prelim/*.arff`. Both PTB protocols recompute the MAD top-2000 feature restriction using training data only.

## legacy8

Expected processed ARFF files: `adenocarcinoma`, `breast2`, `cns`, `colon`, `dlbcl`, `leukemia`, `ovarian`, `prostate`, under `data/comparison/`. The frozen protocol is 10 repetitions of stratified 5-fold cross-validation.

The precise upstream redistribution source for these historical processed ARFF files was not recoverable with sufficient confidence from the stage-6 provenance package. Their exact SHA256 hashes are therefore supplied, and the repository should not claim source-level re-download reproducibility for legacy8 until this provenance is resolved or the processed files can legally be redistributed. This limitation does not affect the frozen numerical results.

## PTB ORA / MSigDB

The corrected ORA uses official Broad MSigDB 2026.1.Hs symbol GMT files for Hallmark, GO Biological Process, and Reactome. Exact files are included under `reference_outputs/ptb_ora/gmt/`; their SHA256 values are recorded in `checksums/MSIGDB_2026_1_HS_SHA256.txt`.

## legacy8 provenance decision for this release

The historical legacy8 processed inputs are retained as a historical reproducibility boundary. Their exact processed-file checksums are recorded in this release. No additional claim is made that these processed inputs constitute a newly curated or newly redistributable source dataset, and no further provenance search is required for this manuscript release.
