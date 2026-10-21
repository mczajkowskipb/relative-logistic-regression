# RLMT-DC — publication-dataset exploratory benchmark, smoke v1

**Pre-outcome registration** for a technical / exploratory reanalysis of
the exact data used in the RLR publication. This is **not** an independent
validation or an extension of GATE-4 Gordon's frozen eligibility screen.

## Motivation and status

The PI prefers first to study RLMT-DC on the **same biomedical matrices
used in the RLR paper**, rather than immediately studying a new GEO
cohort. Accordingly, the GSE72094 GEO-header intake remains preserved
but **paused**; no data download, KRAS labels or training from
GSE72094 are part of this plan.

The RLR publication `main` is immutable:
`4079a5ac400f1d2cfded523c7000f655470925ce`.
The frozen original RLMT-v1 outcomes and code remain independent.
RLMT-DC is a separate opt-in experimental method on
`feature/rlmt-dc-v0`; its synthetic DEV-v3 gate has passed.

## Data world — **all** originally published CompCancer + legacy8

Dataset inventory is predefined without evaluating any RLMT-DC
predictions. All inputs are original processed ARFF files, validated
individually against the **publication-main** SHA256 manifest,
`checksums/DATA_SHA256_EXPECTED.txt`; no random replacement files.

**CompCancer11**:

- Publication's `data/CompCancer_truth_9x10/` shard (9):
  `alizadeh-2000-v1`, `armstrong-2002-v1`, `bittner-2000`,
  `chen-2002`, `golub-1999-v1`, `gordon-2002`,
  `pomeroy-2002-v1`, `shipp-2002-v1`, `singh-2002`.
- Publication's `data/CompCancer_truth_2x2/` shard (2):
  `west-2001`, `yeoh-2002-v1`.
- Canonical hashes in the `data/CompCancer/` section of the
  publication manifest must match their shard copies.

**legacy8** from publication's `data/comparison/` (8):
`adenocarcinoma`, `breast2`, `cns`, `colon`,
`dlbcl`, `leukemia`, `ovarian`, `prostate`.

**Total: 19 distinct historical benchmark datasets.**
No dataset is selected/excluded by observed predictive results.
Neither the earlier GATE4 filter `n>=80, min_class>=25` nor any
posthoc easing is applied in this DIFFERENT, explicitly exploratory
development stage. Inadequate leaf support should yield the method's
existing unsplit fallback; numerical failures must be kept in the
manifest rather than hidden. Biomedical subject independence is not
established by anonymized ARFF sample IDs; legacy8 raw provenance is
incomplete as already documented in `docs/DATA_ACQUISITION.md`.

**PTB**, also included in the RLR paper, will be treated as a
**separate stage** after the same-code CompCancer/legacy8 smoke because
the original PTB fixed and LOCO protocols involve distinct cohort
boundaries and training-only MAD top-2000 preprocessing.
Their original configurations are
`configs/paper_frozen/ptb_fixed.json` and
`configs/paper_frozen/ptb_loco.json`.
Do not silently mix PTB samples into ordinary CompCancer stratified CV.

## Stage A — one fixed outer fold per each of 19 datasets

This is a **code-compatibility and exploratory outcome** pilot,
not a statistical comparison suitable for journal reporting.

- Frozen config:
  `configs/rlmt_dc/published_data_smoke_v1.json`.
- Original published stratified split generator:
  `make_repeated_stratified_folds`, **seed 42**,
  `rep00_fold00` (first fold of 5), one repeat.
- Exactly 19 unique outer holdouts and 57 model fits
  (3 comparators × 19 datasets).
- A fresh, shared 1000-pair **outer-training-only** globally
  supervised pool is computed **once** for each dataset/split.
  This is the same representation family, and method-specific
  discovery OOF p=1000 selection is rebuilt inside its own
  training folds for RLMT-DC.
- Comparison models: published RLR `K10/m10/d2/pool1000,
  alpha=0.5,lambda.min`, original RLMT-v1 with max 2 leaves,
  RLMT-DC-v0 with max 2 leaves.
- The methods receive identical outer training/test arrays.
  Use existing original-model R implementations without patches.
  L1/unsplit exact RLR probability parity is a hard gate.
- Failures are recorded in `run_manifest.csv` as `ERROR`
  with error message and `error.json`, including where small
  minority-class folds make a method nonestimable. Never
  silently exclude, change folds, tune class support or rerun
  only favorable methods.
- All metrics are descriptive:
  paired heldout logloss (primary), Brier, AUC, F1,
  realized leaves, accepted splits, class support and failures.
  No p-value or claim of '19 independent external cohorts'.
- New result root:
  `results/rlmt_dc_published_data_smoke_v1/`.
  Existing paper, GATE3/GATE4 and DEV results are never overwritten.
- A run is **technically complete** only when 57 fit records
  exist and each completed method has intact predictions,
  metrics and repeatable out-of-sample logloss. One failed fit
  makes it **incomplete**, but all intermediate evidence must
  be preserved and analyzed.

This stage is registered after the RLMT-v1 GATE3/4 outcomes were
inspected. Thus the familiar CompCancer cases (including
Chen/Singh/Gordon) are all **development datasets** for RLMT,
and reusing their fold splits is NOT a prospective, fully independent
test of a tuned algorithm. Their inclusion is driven by the PI's
desire for continuity with the publication, not by expected
performance gains.

## Stage B — multi-fold comparison, not yet launched

Only after Stage A technical audit, freeze the complete code SHA,
all failure handling and data checks. Then run, at minimum, all
five originally defined outer folds for **all 19 datasets**,
seed 42, without dataset/quality cherry-picking. For further
publication-style exploratory replication, original CompCancer
used 10x5 for 9 datasets and 2x5 for two; legacy8 used 10x5.
Do not reinterpret folds as independent datasets, and do not
select model variants or hyperparameters from the same test outcomes.

PTB fixed and LOCO should be executed with their canonical
separate evaluation/feature-filtering rules after Stage A.

## Stage C — genuinely independent confirmation (future)

Independent omics cohorts remain necessary for claims of predictive
superiority or generality. Collect case/patient-level provenance,
labels, disjoint groups, preprocessing rules and complete SHA before
making inferential claims. Do not replace the exploratory publication
reanalysis with the GSE72094 KRAS task without a new study plan,
or treat it as a validation of tumor-type classification.

## Implementation and non-authorization boundaries

- NEW runner:
  `R/new/44_rlmt_dc_published_pilot.R`.
- NEW preflight / launcher:
  `scripts/rlmt_dc_published_pilot_launch.py`.
- No existing RLR runner or `R/new/28_rlmt_utils.R` source changes.
- No model selection based on already examined outer-test
  `sample_92` / `sample_80`, Singh, Gordon or GSE72094.
- No '0/0 lookahead = 0 failed candidates' misinterpretation:
  candidate feasibility and split acceptances are separate.
- No automatic train/test dataset redistribution, new GEO fetch,
  or scientific conclusion without visible completed outputs.
