# RLR teacher compressibility: original11 exact-capacity diagnostic — FINAL

**2026-10-29 | Technical PASS; teacher-specific benefit NOT established.**

## Execution / evidence

- [GitHub Actions #40945618900](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40945618900): **12/12 jobs SUCCESS**, 11 Schliep source-reconstructed original-paper cohorts + independent 99 method-by-cohort patient-score audit, frozen execution SHA **8d5a61ba72b4e2633fa2e9fdd034e67109ad19be**.
- [Independent synthetic brute-force oracle check #40945200602](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40945200602): PASS exact global train loss for 0–3 nodes, max 2 comparisons/visited path, rank invariance under within-patient strictly monotonic transforms, teacher-free inference.
- Outer-first T -> disjoint D/V. Every trained selector/predictor uses only D, outer E excluded from fitting/evaluation (frozen original outer stratified partition uses original labels for split generation). Existing full-D original RLR20 compared on V; D-only teacher is genuine 3fold OOF RLR20 with **3 explicitly labeled ridge-logistic surrogate OOF folds in Pomeroy only**.
- Original frozen RLR source and unpublished locked TeX unchanged. Original publication processed ARFF probe/SHA parity with public Schliep reconstructions **NOT VERIFIED**. These eleven previously viewed cohorts are **development**, not independent external biological confirmation.

## Exactly what was tested

A strictly depth-2, at most 0, 1, 2 or 3 internal-node tree from the **same D-only switchBox top20 TSP rank-pair tests**, finite exhaustive optimum. Two independent induction targets:
1. **TEACHER_CAPk:** minimize sum within-leaf squared deviations of cross-fitted original RLR20 **logits** on D.
2. **LABEL_CAPk:** minimize real-label cross entropy on D at the same capacity and same candidate pool.

**For both induced tree topologies**, compute two D-fitted leaf predictions: smoothed true-label frequency for heldout F1/AUC/logloss, and logit(mean D-OOF teacher logits) for heldout probability **fidelity** to full-D original RLR20. The teacher projection on LABEL trees is therefore also trained from D-only OOF RLR, but **does not change LABEL topology**. No parent model runs at autonomous V student inference. Maximum real pair tests per patient is two. Full within-patient transcriptome rank universe is still needed; do not equate this to wet-lab assay costs.

**Train teacher logit SSE fraction** measures retained *in-sample* OOF-logit variance relative to cap0 within THIS D-only top20 candidate pool; it is an empirical optimistic train-fit ratio, **NOT an analytical upper bound on generalization, on different pair pools or other tree architectures**.

## All eleven equally weighted cohort means

| Model | Positive F1 | AUROC | Logloss | Heldout RLR probability fidelity Brier | D train teacher-logit SSE fraction | Actual tested pairs/patient |
|---|---:|---:|---:|---:|---:|---:|
| TEACHER_CAP0 | .321133 | .500000 | .625041 | .196752 | 1.000000 | 0 |
| TEACHER_CAP1 | .707153 | .808167 | .444391 | .096218 | .336213 | 1.0000 |
| **TEACHER_CAP2** | **.809116** | .866362 | .364846 | **.063138** | .274718 | **1.2252** |
| TEACHER_CAP3 | .792276 | .854571 | .391609 | .077477 | **.261321** | 1.4545 |
| LABEL_CAP0 | .321133 | .500000 | .625041 | .196752 | 1.000000 | 0 |
| LABEL_CAP1 | .814901 | .861845 | .377709 | .075045 | .394379 | 1.0000 |
| LABEL_CAP2 | .834371 | .876472 | .369099 | .064907 | .347436 | 1.2140 |
| **LABEL_CAP3** | **.834660** | **.880495** | **.358958** | .064657 | .347338 | **1.3601** |
| Frozen original full-D RLR20 | **.908225** | **.944545** | **.260502** | 0 (self-reference) | — | 17.5455 |

Fidelity Brier is **mean squared difference between projected student-teacher probabilities and original RLR20 probabilities**, NOT Brier vs true class labels. The F1/AUC/LL student metrics use **true D-label leaf probabilities**. Complexity is evaluated rank comparisons, not number of measured probes.

**Primary diagnostic:** increasing TEACHER_CAP2 to TEACHER_CAP3 lowers D train OOF-logit SSE fraction .274718 -> .261321, yet **worsens heldout RLR probability fidelity** .063138 -> .077477 and positive F1 .809116 -> .792276. This shows a genuine train/V disconnect, compatible with overfitting or mismatch between OOF-D and full-D parent, but cannot prove which mechanism explains it.

## Proper paired statistical interpretation (11 cohort units)

| Contrast A − B | Outcome | Mean delta | Exact 2^11 signflip p | 95% bootstrap |
|---|---|---:|---:|---:|
| TEACHER_CAP1 − LABEL_CAP1 | F1 | −.107748 | .242188 | [−.271381,+.039118] |
| TEACHER_CAP2 − LABEL_CAP2 | F1 | −.025255 | .656250 | [−.118494,+.074894] |
| TEACHER_CAP2 − LABEL_CAP2 | Fidelity Brier | −.001769 | .890625 | [−.041693,+.029442] |
| TEACHER_CAP3 − LABEL_CAP3 | F1 | **−.042384** | .531250 | [−.164622,+.070285] |
| TEACHER_CAP3 − LABEL_CAP3 | Fidelity Brier | +.012820 | .582031 | [−.029796,+.051047] |
| TEACHER_CAP3 − LABEL_CAP3 | Actual path comparisons | +.094406 | .500000 | [.000000,+.234266] |

All are **exploratory** and unadjusted; differences are neither statistically demonstrated improvements nor equivalence. Strong CAP3 training gains over CAP0 constant trees are trivial structural improvements shared by LABEL models; they are **not evidence of RLR-specific benefit**. Never treat 99 patient-method groups as independent biological units.

**Heterogeneity, teacher−label F1 at cap3:** Bittner +.3333, Chen +.0534, Singh +.1474; Alizadeh, Armstrong and Gordon 0; Golub −.1091, Pomeroy −.1667, Shipp −.5000, West −.1333, Yeoh −.0912. Bittner's apparent improvement is not generalizable confirmation.

## Science gate

**DIAGNOSTIC COMPLETE; NO new publishable RLR-tree has emerged.** Small relational trees can fit a substantial component of D-OOF teacher logits, but their V performance and RLR-specific benefits are unstable. Do **not** retune depth, candidate pool, logits, thresholds, leaf smoothing, or v4/v5 loss weights based on repeatedly viewed V/E; do not run 470 outer folds for this diagnostic.

**Next materially different step:** original processed ARFF checksum/probe/sample parity and authentic author's CUDA TSPDT (not the TSP-like proxy) on exactly matched patients, along with novelty check vs TSPDT/GTSPDT/EMTTree+RX. Original archives absent in GitHub and require authorized private local access. If a truly new RLR-specific inferential tree mechanism cannot be scientifically defended, the original unpublished RLR (with native kTSP10/15/20 and matched-K comparisons) remains more mature for publication. No manuscript edit; no automation.

## Permanent machine-readable evidence

- [All 99 original11 model–cohort metrics](results/RLR_CAPACITY_V1_ALL11_METRICS.csv)
- [All preregistered 11-cohort exact signflips + bootstrap](results/RLR_CAPACITY_V1_PAIRED_INFERENCE.csv)
- [Frozen diagnostic method](RLR_TEACHER_COMPRESSIBILITY_FRONTIER_PROTOCOL_V1.md)
- [Complete independently audited Actions run](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40945618900)
