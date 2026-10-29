# RLR-TREE — next scientific decision gate, after v5 (2026-10-29)

**DECISION TYPE: research strategy / evidence inventory.** This is NOT a new classifier, not v6 and not a claimed methodological advance. No hyperparameter changes, no new E evaluation, no changes to locked original RLR/manuscript. PI's newer tree objective supersedes 2026-10-26's older Issue #6 body.

## Result that constrains all next work

Real audited GitHub runs: [WT v4](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40704733286) 12/12 PASS and [GLVT v5](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40919146672) 12/12 PASS. On the exact 11 original-named, **public Schliep-reconstructed D/V development cohorts**, v5 RLR-logit global depth2 tree: **positive F1 0.826367, AUC 0.878139, LL 0.364068, 1.3082 rank-pair comparisons/patient**. Same optimizer, same TSP-top20 candidates but **no RLR**: F1 0.834660, AUC 0.880495, LL 0.358958, 1.3601 checks. Simpler greedy label TSP-like: F1 0.840313, LL 0.349497, 1.2987 checks. Original RLR20: F1 0.908225, LL 0.260502, 17.5455 comparisons; real native switchBox kTSP5: F1 0.863363, LL 0.289605, 5 checks. One RLR-worsening F1 cohort (Yeoh), ten ties; no teacher-specific systematic gain. Full report: [RLR_GLVT_V5_ALL11_COMPLETE_NEGATIVE_SCIENCE_GATE.md](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-global-logit-variance-tree-v5/docs/RLR_GLVT_V5_ALL11_COMPLETE_NEGATIVE_SCIENCE_GATE.md).

It is **NOT** established that a <=2-pair decision tree can represent original RLR20 faithfully while maintaining its patient-level calibration. Inadequate expressive capacity and inadequate induction are distinct hypotheses; do not equate apparent failures with impossibility, but do not assert novelty because a tree is interpretable.

## Historical novelty ledger — claims we must NOT repeat

| Mechanism | Already present / tested | Consequence |
|---|---|---|
| Tree nodes testing gene-expression pair orders | Original author's **TSPDT (2011)** | Not novel |
| Global evolutionary search, variable tree, multi-pair votes in a node | Author's **GTSPDT (2013)** | Not novel |
| Global optimization, mixed relational tests, complexity penalties | Author's **EMTTree+RX** unpublished/published-in-habilitation series; inspect actual manuscript before novelty claim | Not sufficient |
| OOF RLR probability guidance for pair tree splits/leaves | In-house RLR-DT v1/v2/v3 negative/weak | No stable RLR-specific gain |
| Two original RLR beta-weighted ternary relations per node | WT v4 negative | No stable value of beta weight magnitudes |
| Exact global small-tree optimum with OOF RLR logit within-leaf variance | GLVT v5 negative | Global optimizer/teacher is not the solution |
| Certified sequential partial RLR evaluation | C-RET already investigated and strong fixed ordering control | Do NOT recycle as a novel autonomous tree |
| Distilled trees, probabilistic scoring lists, optimal decision trees | Independent published prior work | Need genuinely new identifiable RLR claim, not generic distillation/optimization |

Primary source checks:
- TSPDT: https://pubmed.ncbi.nlm.nih.gov/21431543/ (DOI 10.1007/978-1-4419-7046-6_3).
- GTSPDT: https://www.aragorn.wi.pb.edu.pl/~mczajk/eurogp2013.pdf ; https://evostar.org/2013/cfpEuroGP.html .
- Generic reproducible model distillation: https://link.springer.com/article/10.1007/s10994-024-06597-w .
- Probabilistic scoring lists: https://doi.org/10.1007/s10994-024-06705-w .
- Original unpublished RLR logistic signed ternary m10 relations are specific to the parent's design, but merely including their beta weights in tree nodes (v4) is insufficient.

## Provenance gates — found now versus missing

1. **Original ARFF / within-probe parity:** original processed `data/CompCancer_truth_9x10` (9 sets) and `data/CompCancer_truth_2x2` (2 sets) absent from GitHub analysis branch; current 11 cohorts are publicly reconstructed Schliep matrices. Existing read-only inventory: `scripts/rlr_original_paper_arff_sha256_preflight_v1.py` emits SHA256/size/path (no patient data). Running it on an author-owned original data tree is necessary before an *exact original-paper reproduction* claim. Original PTB fixed, LOCO and legacy8 source ARFF also currently absent from GitHub; do not confuse them with 11 CompCancer sets.
2. **True historical TSPDT CUDA:** code archive `cudatsp4(1).zip` and prior adapter v2 are described in old research handoff, but **not present as runnable source in the present GitHub branch**. Never call the greedy binary TSP-like control the historical CUDA method. Requires author-controlled CUDA files in trusted local environment and real heldout same-patient predictions, with no unsafe private code publication.
3. **Exposure:** all 11 cohort identities/results have been repeatedly consulted in earlier tree development. Cannot claim blind independent validation even after recomputing with original split. No fresh claim without unexposed data or prospective locked-split comparison.
4. **Public manuscript visibility:** the public GitHub repository's `main` branch contains `manuscript_lock/RLR3_submission_final.tex` (and PDF). Original article is **not necessarily journal-published**, but its full source is publicly accessible; PI should check journal preprint/anonymity and release strategy before submission. Do not alter/remove locked file without explicit authorization.

## Immediate research order / exit criteria

**GATE A — Comparator/data parity (next actual engineering task).** Locate/provide the original processed ARFF and CUDA archive *within an authorized private working environment*, verify SHA/probe/sample/positive-class/split mapping, run authentic CUDA TSPDT vs original RLR20 and native switchBox kTSP on exactly matching heldout patients. Do not route unlicensed/private CUDA into the public repo. Minimum: signed summary of source availability, reproducible build environment, genuine TSPDT predictions and identical fold patient IDs (or explicit BLOCKED_MISSING_SOURCE).

**GATE B — Feasibility before v6 algorithm.** Compute a **training-only, bounded teacher-compressibility frontier** (diagnostic, NOT a new classifier or publication claim): for fixed training D, candidate rank-pair pools and budgets of 0/1/2/3 internal nodes, evaluate the *best attainable* surrogate fidelity/label error under equal tree capacity; compare actual TSPDT and small label trees. Use OOF parent RLR20 predictions where they define teacher targets; include honest heldout V solely for reporting. An upper bound of available fidelity and direct capacity-matched controls can distinguish representation insufficiency from poor induction. No selecting architecture on old exposed E. Keep CPU pilot small.

**GATE C — New RLR-specific falsifiable mechanism, not parameter search.** After A+B, propose a single formally defined rank-logistic decision mechanism whose property is absent from TSPDT, GTSPDT, EMTTree+RX and C-RET. Minimum proof obligations: why original ternary logistic coefficients, not any classifier, are essential; what the tree can do better than the strongest identical-budget no-RLR control; how the mechanism generalizes under train-only splits. Reject speculative names/formulas without an explicit ablation and a measurable primary endpoint.

**GATE D — Research release.** Use independently verified original publication data and (ideally) unexposed cohorts; choose a single design and freeze it before tests. 470 folds ONLY after a positive structural/scientific pilot. Do not stall original unpublished RLR manuscript submission indefinitely if GATE C fails.

Current status: **A missing original source and authentic TSPDT benchmark; B unexecuted; C not yet justified; D blocked.** V1–v5 source and negative evidence preserved. **No new v6 training launched; no automation.**
