# Binary rank-pair comparators: 11 cohorts × 5 held-out folds — full v2 pilot

Date: 2026-10-30. **Real GitHub Actions #41064335985: 12/12 SUCCESS** (11 cohorts × 5 outer folds × 3 models, plus independent patient-level auditor). [Run and artifacts](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/41064335985).

**Data provenance:** these v2 numbers use PUBLIC SCHLIEP RECONSTRUCTED matrices; original publication processed ARFF SHA256 parity was **NOT VERIFIED** at this pilot stage. One full outer five-fold repeat, not historical 470-fold paper schedule, and 11 datasets are inferential units, not 55 folds.

The three models share *identical training-only switchBox pair candidates* (top1000, K20 gene-disjoint pairs, d=2):
- Original unmodified RLR20: ternary rank-gap m10 + glmnet elastic net.
- Binary elastic net (PurIST-inspired methodology): binary rank/order indicator + glmnet.
- Binary MCP (DeCAF-inspired methodology): binary rank/order indicator + ncvreg.
Controls are NOT the complete clinical original PurIST / DeCAF algorithms. The original full DeCAF method has an additional gene screen involving independent cohorts, and original PurIST is trained for different disease-specific subtype labels.

## Audited dataset-equal means (55 outer folds / 165 fits)

| Method | F1 ↑ | AUC ↑ | Logloss ↓ | Brier ↓ | ECE ↓ | Active coefficients |
|---|---:|---:|---:|---:|---:|---:|
| Original RLR20 | **0.876180** | **0.961673** | **0.237944** | **0.062638** | **0.092092** | 16.818182 |
| Binary ENET20 (PurIST-inspired) | 0.823985 | 0.948487 | 0.290442 | 0.077384 | 0.110474 | 16.745455 |
| Binary MCP20 (DeCAF-inspired) | 0.813641 | 0.947843 | 0.291233 | 0.085455 | 0.134420 | 11.527273 |

## Exact 2^11 two-sided sign-flip tests, across 11 dataset-level 5-fold means

RLR minus binary comparator. Lower is better for logloss/Brier, higher for AUC/F1. BH is **exploratory across all eight tests**.

| Comparator | Metric | Δ RLR − comparator | Raw p | BH q |
|---|---|---:|---:|---:|
| Binary ENET | logloss | −0.052498 | 0.034180 | 0.039062 |
| Binary ENET | brier | −0.014746 | 0.008789 | 0.023438 |
| Binary ENET | auc | +0.013187 | 0.078125 | 0.078125 |
| Binary ENET | positive_f1 | +0.052194 | 0.031250 | 0.039062 |
| Binary MCP | logloss | −0.053290 | 0.026367 | 0.039062 |
| Binary MCP | brier | −0.022817 | 0.001953 | 0.015625 |
| Binary MCP | auc | +0.013830 | 0.031250 | 0.039062 |
| Binary MCP | positive_f1 | +0.062539 | 0.007812 | 0.023438 |

Seven of eight exploratory differences have BH q < 0.05 (all except binary ENET AUC). **This is not preregistered confirmation**; earlier pilot on first outer training fold already showed encouraging differences, the same general cohorts have been used repeatedly, and exact processed paper ARFF parity not yet established.

## Independent audit and full results

Run #41064335985 successfully audited:
- 165 fold-method metric summaries replayed from 11 cohort held-out per-patient probability artifacts.
- Every patient appears exactly once in test data per method within its dataset.
- Logloss, Brier, F1, AUC and ECE replayed independently from probabilities.
- Exact sign-flip enumerations (2048) and exploratory BH correction.
- Original frozen model source unmodified.

Audit artifact file: `rlr-binary-outer5-all11-audit-41064335985` (GitHub Actions artifact, retention 30 days), containing `all11_summary.csv`, `all11_paired_exact_inference.csv`, `all11_165_outer_fold_metrics.csv` and `all11_outer_patient_predictions.csv`.

## Next reproducibility gate

User reports the **canonical original processed source directory** exists on Ubuntu at `/opt/2026/RLR/data/CompCancer/` and the paper verification manifest at `/opt/2026/RLR/relative-logistic-regression/checksums/DATA_SHA256_EXPECTED.txt`. Prepared read-only hash gate and **separate local original-ARFF benchmark**, requiring checksum parity before any fitting:
- `scripts/rlr_original_compcancer_sha_preflight_v2.py`
- `R/experimental/rlr_binary_outer5_original_arff_v2.R`
- `scripts/run_rlr_original_arff_outer5_server_v2.sh`

These are prepared source files; **server execution has not yet been performed or verified**. Do not yet amend the submission manuscript with novel superiority claims; report full binary control method differences and exact source provenance first.
