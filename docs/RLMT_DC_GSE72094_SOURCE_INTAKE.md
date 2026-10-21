# RLMT-DC — independent omics cohort intake: GSE72094

**Status: candidate / provenance-audit-only. NOT eligible for RLMT-DC training yet.**
Prepared after the completed synthetic DEV v3, without inspecting
any RLMT-DC performance on this cohort.

## Why this dataset is a useful prospective candidate

Primary source:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72094

- Study: Eschrich/Schabath et al., gene expression and mutation biology
  in resected lung adenocarcinoma, GEO accession `GSE72094`,
  BioProject `PRJNA292925`; public since 2015-10-21.
- GEO lists **442 samples** of human lung adenocarcinoma,
  platform `GPL15048` (custom Affymetrix HuRSTA array),
  processed expression per sample and a Series Matrix file.
  Source describes log2 IRON-normalized expression and correction
  for RNA quality; raw CEL archive is ~2.1 GB, **not** required for
  initial metadata discovery.
- **Proposed biological task:** predict *assessed KRAS mutation status*
  (mutated versus wild type) from transcriptomic expression only.
  This is an expression-to-genotype association task, not cancer diagnosis.
- A subsequent published analysis describes **398 usable patients**,
  **139 KRAS-mutant / 259 KRAS-wildtype**:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10292897/
  **These are literature figures, NOT yet verified from downloaded
  primary accession-specific metadata**. The difference 442–398
  cannot automatically be attributed to KRAS missingness, because
  downstream studies may also apply other clinical exclusion criteria.
- The candidate was not part of the published RLR CompCancer, legacy8
  or PTB cohorts as listed in this repository. **Overlap with those
  cohorts at patient level remains to be checked** (different study
  identifiers alone do not prove participant independence).

## Mandatory source checks BEFORE any model evaluation

1. Retrieve GEO `GSE72094` Series Matrix or primary metadata and
   record the official source URL, retrieval time and file SHA256
   once the **complete** original file is downloaded.
2. Confirm that every retained expression column maps to a unique
   `GSM` and a defensible unique **patient** key. If multiple samples
   per patient, predefine grouped train/validation/test splits.
3. Retrieve KRAS mutation labels from the **primary submission or an
   explicitly traceable supplementary mutation table**. Verify
   exact identifier joins and all categorical values; distinguish
   assayed wild type from missing/unknown/NA. Do not derive labels
   from expression or model predictions and do not invent them.
4. Count labeled, unlabeled and ambiguous cases independently.
   Verify the 139/259 figures only if the same cohort/exclusion
   definitions apply; discrepancies require a documented explanation,
   not manual relabeling.
5. Inspect complete expression matrix, feature IDs, sample
   duplication, numeric/nonfinite cells, quality-control probes
   and any failed or repeated profiles; audit preprocessing and
   distribution before choosing a fixed transformation.
6. Verify public use/redistribution conditions with the original
   source and study's data-sharing terms; citation is required.
   Do not redistribute private patient identifiers or clinical data.
7. Record source-to-processing code, hashes and intended group
   splitting. Any feature filtering (e.g. unsupervised top-2000 MAD)
   must be trained inside *each outer-training fold* only; **not**
   performed across all samples before CV.
8. Check whether the cohort contains diagnostic/technical confounders,
   batch effects or class imbalance that require a predeclared,
   non-leaky handling policy.

A candidate's GEO sample count is not equivalent to the number of
valid labeled patients.

## Minimal intake, no expression download and NO training

Script: `scripts/rlmt_dc_geo_header_intake.py`.

- Stream only the header of the public GEO Series Matrix,
  stopping before the expression table; count `GSM` accessions and
  enumerate sample metadata fields/characteristics.
- Audit whether a *direct* KRAS metadata field even exists. Do not
  use free-text guesses to assign class labels.
- Create a report in a NEW `results/` directory, without
  modifying the primary data, manuscript, frozen method, or pilots.
- A header-only inspection **cannot** establish the SHA256 of the full
  compressed matrix, completeness of expression values, a true
  patient map, mutation truth, reuse permissions or cohort overlap.

After intake, decision must be either:
`ACCEPT_SOURCE_PENDING_FULL_AUDIT` or `BLOCK_SOURCE`.
No performance evaluation may be triggered by the intake script.

## Scientific evaluation constraints (not launched)

A confirmatory superiority claim cannot follow from a single
new cohort. Seek further truly independent datasets with a coherent
biological endpoint and group-aware labels. Freeze cohort inventory,
feature processing, base RLR, original RLMT-v1, RLMT-DC-v0, external
test handling and failure policy *before* predictive scores are viewed.
Compare paired test logloss at the **dataset** level, with Brier and
worst-case tail loss, and report zero-split outcomes without cherry-pick.
This remains an exploratory candidate until an evidence-backed
source/label/subject audit and a prospective comparison plan pass.
