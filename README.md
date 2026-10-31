# Relative Logistic Regression (RLR)

Publication reproducibility release for:

**Relative Logistic Regression: A Probabilistic Refinement of Hard-Voting Pair-Based Omics Classification**

Marcin Czajkowski, Anna Czajkowska, Krzysztof Jurczuk, Marek Kretowski.

## What is frozen

The primary RLR models use an upstream pool of 1000 ranked relations, rank margin `m=10`, strict gene disjointness `d=2`, target shortlist sizes `K=10,20,50`, elastic-net `alpha=0.5`, `lambda.min`, and seed 42. `d=1` variants are technical bridge/control configurations, not primary paper models.

Evaluation roles:
1. **CompCancer** — primary 11-dataset inferential benchmark.
2. **PTB fixed** — descriptive structural case study on a predefined discovery-data 103/44 train/test split.
3. **PTB LOCO** — descriptive four-cohort cross-cohort stress test.
4. **legacy8** — historical positioning benchmark.

PTB fixed is **not** external validation: GSE83456 is the source benchmark's 106-sample external-validation cohort but is not the canonical test set used for the PTB fixed results reported here.

## Repository layout

```text
R/new/                  canonical pipeline/reporting/supplementary analysis code
configs/paper_frozen/   portable copies of canonical truth-rebuild configs
docs/                   protocol, provenance and data documentation
environment/            R/system/package-version capture
checksums/              original-source, publication-source, data and MSigDB hashes
reference_outputs/      selected corrected PTB ORA outputs and GMT inputs
scripts/                 SHA256 preflight, reproduction, independent score audit and smoke helpers
```

## Primary reproduction

Place processed data under the layout described in `docs/DATA_ACQUISITION.md`, verify SHA256, then from repository root run:

```bash
bash scripts/smoke_test.sh
bash scripts/reproduce_primary.sh
```

The CompCancer state is intentionally assembled from two execution shards: nine datasets at 10x5 and West/Yeoh at 2x5. The collector scans the merged output tree to produce final 11-dataset summaries.

## Original-ARFF integrity preflight (Supplement S13)

The public repository intentionally does **not** contain old submission manuscript PDFs or the third-party processed gene-expression matrices. The current manuscript and supplement are supplied to the journal separately. Archived code and data provenance remain under version control.

To check exact byte identity against all **45 entries in the frozen source manifest**, and additionally verify the **11 CompCancer 9×10 / 2×2 shard copies** against the SHA256 values of the corresponding original ARFFs:

```bash
python3 scripts/verify_original_arff_preflight.py --data-root /path/to/processed/RLR/data
```

This read-only preflight fails when a file is missing or has a different SHA256. Exact ARFF byte identity also establishes that the original feature and class schema is unchanged.

## CompCancer ternary-vs-binary ablation (Supplement S7.1)

The **55 outer-fold (11 datasets × 5 folds, three models)** fixed-pair controlled ablation is reproducible with:

```bash
bash scripts/reproduce_compcancer_ternary_binary_ablation.sh /path/to/processed/RLR/data
```

This separate benchmark runs the original RLR-20 ternary margin representation against binary ENET and MCP **methodological controls**, not against the full published clinical PurIST or DeCAF classifiers. It requires R packages `foreign`, `switchBox`, `glmnet`, `ncvreg`, `pROC`, `jsonlite`; source ARFF files are checksum-gated and not redistributed. Original model code and 470-fold primary results are unchanged. The runner writes held-out per-sample predictions and independently replays logloss, Brier, AUC, F1, ECE and exact 11-dataset sign-flip/BH inference. See `docs/COMPCANCER_TERNARY_BINARY_ABLATION_REPRO.md`.

## Corrected PTB ORA

After reproducing or supplying the canonical PTB fixed output tree:

```bash
bash scripts/reproduce_ptb_ora.sh
```

The ORA uses only `01_train_data.arff` (103 samples), a training-only MAD top-2000 universe, and the bundled official MSigDB 2026.1.Hs GMT files. No RLR/k-TSP model is refit by the ORA script.

## Verification

`paper_numbers_manifest_v3.csv` / `.xlsx` map displayed paper numbers to source outputs. `scripts/check_release.py` validates configs, hashes and known protocol invariants. `scripts/verify_reference.py` can compare selected paper-table/statistics CSVs between two supplied reference trees.

## Reproducibility boundary

No additional benchmark reruns are required for the numerical lock. The software is released under the MIT License. For the historical legacy8 benchmark, this release retains the explicit checksum-based processed-input provenance boundary documented in the repository; no further source recovery is required. The clean publication repository is the authoritative public reproducibility release associated with the manuscript.
