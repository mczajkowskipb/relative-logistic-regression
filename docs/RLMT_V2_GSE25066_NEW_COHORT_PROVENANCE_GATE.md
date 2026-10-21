# GEO candidate dossier — GSE25066, NOT cleared for prediction yet

## Public-source verification (2026-10-21)

Primary official record:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066

Official NCBI GEO overview describes it as a SuperSeries comprising
GSE25055 (DISCOVERY) and GSE25065 (VALIDATION); platform GPL96,
Affymetrix HG-U133A, human breast cancer and neoadjuvant
taxane–anthracycline response/survival. The 2011 series has
508 GSM public samples in total (310+198).

Sources:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25055
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25065
- https://www.ncbi.nlm.nih.gov/geo/info/download.html
- PMID 21558518 (original study).
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4148546/
  illustrates pCR vs RCB class outcome ambiguity and reports
  that a later analysis used 286 rather than the 310 discovery
  series cases. DO NOT silently substitute 286 for 310.

**IMPORTANT** The source is independently nominated in a new study,
not a never-before-researched GEO accession globally. Original
discovery/validation are SUBSERIES OF THE SAME STUDY, not two
fully independent clinical institutions by definition. Public
GSM uniqueness is NOT proof of unique patients or no clinical
overlap. This validation cohort is an intended heldout group
with verification obligations, NOT an already established
external patient-independent test for our newly developed v2
algorithm.

Different studies report different usable sample counts due
to clinical outcome availability and exclusions. The official
series record gives 310/198 physical GSM arrays; it does NOT
prove 310/198 complete binary pathological-response labels.
The exact primary outcome is NOT frozen until clinical fields
and publication source are inspected. In particular:
- A predictive model output called `dlda30_prediction`
  is NOT the observed clinical pCR truth.
- A predicted RCB outcome is NOT the true observed RCB.
- Distant relapse after followup cannot be directly recoded
  into pCR/non-pCR without new endpoint semantics.
- RCB-0 alone and RCB-0/I are DIFFERENT endpoint definitions.
- ER/HER2 status restriction and treatment protocol must
  agree with the chosen clinical endpoint before inclusion.
- Probable submitter preprocessing is MAS5/log2; verify
  expression scale from GEO Sample_data_processing rather
  than applying a new global normalization using test data.

## Implementation prepared

\`scripts/rlmt_v2_gse25066_source_gate.py\` downloads TWO exact
submitter processed series matrices, records compressed SHA256,
checks expression column GSM IDs, sample counts, within- and
between-GSM overlap and the GEO metadata fields, then creates
\`PROVENANCE_METADATA_ONLY.json\`. It records *no output labels*,
does *no model fitting* and does *not merge* the two samples.

The script intentionally forbids overwriting a partial/incomplete
download. If it writes \`INCOMPLETE_DO_NOT_USE\`, do not treat
the cohort as available; determine source cause and open a
NEW named evidence root for a retried acquisition.

Discovery source:
\`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25055/matrix/GSE25055_series_matrix.txt.gz\`

Validation source:
\`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25065/matrix/GSE25065_series_matrix.txt.gz\`

Both URLs follow official NCBI GEO FTP series matrix
layout. Network access or gzip format on operator Ubuntu
has NOT been independently verified by this repository
commit. No cloud API/token access is required.

## Pending mandatory clinical/independence gate

Before even requesting the *first* actual biomedical v2 fit:

1. Match subject identifiers (NOT just GSM) against
   GSE25055/GSE25065 sample clinical tables and original paper.
2. Identify direct observed pCR/RCB truth, its exact definitions,
   available/missing labels and discard/exclusion counts.
3. Verify overlap with any previously run patient material,
   discovery vs validation and duplicate biopsies.
4. Register preanalytic processing choices and train-only
   probe selection / gene mapping. RLR uses within-sample
   rank relations and should not implicitly consume
   cross-patient test-normalized features.
5. Document public use/access terms and original authorship,
   platform/processing/reuse rights. GEO is publicly
   downloadable; do not infer a specific permissive
   license from openness.
6. Freeze cohort/sample manifests, binary label mapping,
   grouped train/test identities, SHA256 and primary
   logloss + Brier/BA/recall/calibration analysis **before**
   model comparison outcomes.
7. Maintain complete fixed results regardless of positive
   or negative performance; no heldout-driven retuning.

**No biomedical predictive training is authorized by this
metadata-only preparation stage.** External scientific
claims remain blocked until all checks are completed.
