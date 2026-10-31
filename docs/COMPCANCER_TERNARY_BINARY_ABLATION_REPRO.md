# Supplement S7.1 — original-CompCancer-ARFF controlled encoding ablation

**Status:** independent, owner-executed run on the original frozen CompCancer ARFF snapshot passed SHA256 and all fit/patient audit gates on 2026-10-30. This page documents that completed historical computation and makes its code reproducible in the publication release. No new local re-fit is claimed for this GitHub branch.

## Objective and scope

Isolate the role of **margin-aware ternary pair encoding** in RLR relative to conventional **binary pair order** indicators under matched **training-fold-only** pair selection. The binary controls use the same exact K=20, d=2 shortlist and differ only in encoding and/or model penalty:
- Original RLR-20: signed rank-margin ternary states {−1,0,+1}, m=10, elastic net alpha=0.5, lambda.min.
- Binary-ENET-20: binary pair order indicator, elastic net alpha=0.5, lambda.min; same pair set.
- Binary-MCP-20: binary pair order indicator, ncvreg MCP (alpha=0.05), same pair set.

The binary controls are **NOT** the clinical pretrained PurIST or DeCAF classifiers. These predict different disease-specific subtype labels. This is a matched **representation/penalty ablation**, not a clinical head-to-head comparison with their original implementations.

## Exact data identity, reproducible invocation

Data are **not** stored in GitHub. The author-owned source tree has 45 frozen file SHA256 entries, including eleven CompCancer main ARFF files and their original 9×10 / 2×2 execution-shard counterparts. Execute from the repository root with an authorized copy of this data directory:

```bash
python3 scripts/verify_original_arff_preflight.py --data-root /path/to/processed/RLR/data

bash scripts/reproduce_compcancer_ternary_binary_ablation.sh /path/to/processed/RLR/data
```

Prerequisites: Python 3, R with foreign, switchBox, glmnet, ncvreg, pROC, jsonlite; Git and sha256sum. Both commands are standalone child processes and do not modify the original ARFF files. The runner verifies every original SHA256 before starting and each target + shard once again before fitting. It uses original RLR functions in `R/new/07_world_protocol_utils.R` and an unchanged copy of the original `rg_parent()` helper embedded in `R/new/28_compcancer_ternary_vs_binary_ablation.R`. No `R/experimental/` dependency is required.

Expected outputs under `results/rlr_binary_outer5_original_arff_v2/` include each dataset's `fold_metrics.csv`, `patients.csv` and `fold_manifest.csv` from five external test folds. Under `results/rlr_binary_outer5_original_arff_audit_v2/`, the independent Python auditor writes all 165 method-fold metric summaries, all heldout predictions, 33 method-by-dataset means, 11-weighted method means, exact paired 2^11 sign-flip p-values and eight-test exploratory BH q-values, plus a source-hash attestation. The runner refuses to claim PASS if the independent patient-level audit or original data SHA256 recheck fails.

The additional experiment uses **one five-fold outer repeat per dataset** (11×5=55 outer folds; ×3 methods = 165 fits), not the full 470-fold primary benchmark; inference uses **11 independent dataset means**, never 55 pseudoreplicates. All train-specific pair nomination and regularization selection occur within the corresponding outer training portion.

## All 11 original-dataset equal-weighted means (owner-completed audit)

| Model | Logloss ↓ | Brier ↓ | AUC ↑ | Positive F1 ↑ |
|---|---:|---:|---:|---:|
| Original RLR20 | **0.237943894** | **0.062637676** | **0.961673479** | **0.876179704** |
| Matched binary ENET20 | 0.290442282 | 0.077383727 | 0.948486850 | 0.823985217 |
| Matched binary MCP20 | 0.291233444 | 0.085454758 | 0.947842990 | 0.813640710 |

Paired exact two-sided tests on **11 dataset-level means**. Delta = RLR minus binary control; lower is better for logloss/Brier, higher for AUC/F1.

| Comparator | Metric | Mean delta | Exact p | Exploratory BH q |
|---|---|---:|---:|---:|
| Binary ENET | Logloss | −0.052498 | 0.0341796875 | 0.0390625 |
| Binary ENET | Brier | −0.014746 | 0.0087890625 | 0.0234375 |
| Binary ENET | AUC | +0.013187 | 0.0781250000 | 0.0781250 |
| Binary ENET | F1 | +0.052194 | 0.0312500000 | 0.0390625 |
| Binary MCP | Logloss | −0.053290 | 0.0263671875 | 0.0390625 |
| Binary MCP | Brier | −0.022817 | 0.0019531250 | 0.0156250 |
| Binary MCP | AUC | +0.013830 | 0.0312500000 | 0.0390625 |
| Binary MCP | F1 | +0.062539 | 0.0078125000 | 0.0234375 |

These exploratory q-values relate only to this **separate eight-test** family, not to the 15 Wilcoxon tests of the primary full-470 matched-kTSP benchmark. Seven of eight exploratory BH q-values are <0.05; binary ENET AUC is not significant. A few MCP fits emitted `Model saturated; exiting...` warnings. The result is **not** evidence of superiority over the original complete DeCAF/PurIST models; an authentic frozen-DeCAF comparison on its own CAF-subtyping task did not show RLR superiority. Because these cohorts have previously been analyzed extensively, inference is exploratory, not preregistered independent confirmation.

## Provenance and scientific boundary

- Main manuscript results and original RLR source are unchanged by this release update.
- Original processed data SHA256 verification was carried out **locally** by the authors; Github does not expose those matrices.
- Existing historical 470-fold CI development runs used public reconstructed Schliep source data, not bytes independently checked against the author's original local processed files at the time. Do not conflate those activities.
- The old RLR3 manuscript and stale PDF hash were intentionally removed from the current repository tree to avoid presenting an obsolete document as the journal submission. Older snapshots remain recoverable from Git history.
