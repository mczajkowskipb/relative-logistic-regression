# RLR-WT v4 complete 11-cohort D/V SCIENCE GATE — NOT a predictive improvement

**Completed 2026-10-28.** [Actual 11-cohort GitHub Actions run #40704733286](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40704733286): 11 cohort jobs + independent patient prediction auditor = **12/12 SUCCESS**; 77 real held-out patient model×cohort summary metrics recalculated from CSV, outer E UNREAD, D/V IDs disjoint, true 0/2/4 weighted-node comparisons checked, exactly seven arms per cohort. [RLR-WT synthetic test #40704847792](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40704847792) SUCCESS, including a constructed 2-pair sign-conflict where true β magnitudes and sign-only weights take different branches, plus nonlinear per-patient monotone invariance.

**Data are PUBLIC reconstructed Schliep expression matrices under original paper cohort names; original processed ARFF SHA/probe identity NOT VERIFIED, and all cohorts were repeatedly exposed during research development.** These are exploratory D/V experiments, not independent or prospective validation.

## 11 equal dataset means (first original outer T, inner D/V only)

| Arm | Positive F1 ↑ | Macro F1 ↑ | AUROC ↑ | Logloss ↓ | Mean actual pair tests/patient |
|---|---:|---:|---:|---:|---:|
| WT_BETA_LABEL | .749911 | .783075 | .811437 | .488595 | 2.4304 |
| WT_BETA_OOF | **.757869** | .790682 | .819922 | .469701 | 2.4304 |
| WT_SIGN_LABEL | .701958 | .785712 | **.829490** | **.448918** | 2.4818 |
| WT_SIGN_OOF | .701958 | .785712 | .826460 | .455898 | 2.4818 |
| **LABEL_BINARY_TSP_TOP20 (same tiny tree limits)** | **.840313** | **.864467** | **.881017** | **.349497** | **1.2987** |
| Real native switchBox kTSP5 | .863363 | .881064 | .937638 | .289605 | 5.0000 |
| Frozen original RLR20 | .908225 | .925915 | .944545 | .260502 | 17.5455 |

**Global size:** WT_BETA_OOF has average **1.455 internal nodes**, **2.455 leaves**, **2.545 distinct original ternary pair relations** and **5.091 distinct genes** across all 11; actual mean patient path 2.430 pair tests. The label-only TSP binary tree has 1.636 internal nodes, 2.636 leaves, 1.636 distinct relations and 3.182 distinct genes (1.299 mean tests). The two-relation node is **worse in both mean F1 and mean evaluated relations** than the simpler one-pair control. Both approaches still require full assay ranks to determine ranks.

## Predeclared paired contrasts; exact 2^11=2048 signflip on 11 cohort units + seeded 8k dataset bootstrap (exploratory)

| Contrast A − B | Outcome | Mean delta | 95% bootstrap | Two-sided exact p |
|---|---|---:|---:|---:|
| WT_BETA_OOF − WT_BETA_LABEL | F1 | **+.007958** | [.00000,+.02387] | **1.0000** |
| WT_BETA_OOF − WT_BETA_LABEL | AUC | +.008485 | [.00000,+.02545] | **1.0000** |
| WT_BETA_OOF − WT_BETA_LABEL | LL | −.018894 | [−.05668,.00000] | **1.0000** |
| WT_BETA_LABEL − WT_SIGN_LABEL | F1 | +.047953 | [−.03517,+.15556] | .46875 |
| WT_BETA_LABEL − WT_SIGN_LABEL | AUC | −.018053 | [−.08909,+.04649] | .609375 |
| WT_BETA_LABEL − WT_SIGN_LABEL | LL | +.039676 | [−.09204,+.18516] | .625 |
| WT_BETA_OOF − LABEL_BINARY_TSP_TOP20 | F1 | **−.082444** | [−.18161,+.00203] | .2500 |
| WT_BETA_OOF − LABEL_BINARY_TSP_TOP20 | AUC | **−.061095** | [−.16116,+.00888] | .2500 |
| WT_BETA_OOF − LABEL_BINARY_TSP_TOP20 | LL | **+.120204** | [−.01562,+.31548] | .21875 |

**Most important negative finding:** WT_BETA_OOF improves vs same weighted β label-only on F1 **in CHEN ONLY** (+0.0875), **ten other datasets tied** exactly. Its positive average effect is misleading; signflip p=1 and no independent evidence of teacher benefit. WT_BETA_LABEL vs SIGN_LABEL is heterogeneous: a large rare-small Pomeroy +.5 F1 contributes to positive mean, but in Bittner −.1 and Chen −.17. The main beta+OOF model underperforms the simpler label-TSP control in Alizadeh (−.4444), Bittner (−.2667), Shipp (−.25), and is roughly tied in most of the rest. No statistical proof of equivalence or consistently beneficial weights.

## Scientific decision

1. **DO NOT promote this RLR-WT v4 as a confirmed improved classifier or new method.** Two-pair weighted RLR nodes mathematically function and are standalone, but carry more relation tests and less mean accuracy than simple TSP-like binary tree. Original RLR remains the valid unpublished paper foundation.
2. **DO NOT run first outer E or full 470 on v4**: D/V validation failed the prespecified competitive-complexity check. Do not tune β, α, pool size, K, thresholds, branch sign or penalty using these known validation outcomes.
3. **DO NOT claim novelty** of two-pair weighted nodes versus author's TSPDT/GTSPDT, 2020 Tree-based ARXA and EMTTree+RX. Those already use varied multirelational node tests and optimization. RLR-to-tree student benefit is unproven.
4. The testable potential research opening is **different from simply expanding node capacity**: select conditional rank evidence on a patient path by a *globally stable*, fully probabilistic objective with an explicit uncertainty/abstention guarantee; benchmark against C-RET certified fixed-order inference, original RLR, and a strictly capacity-matched TSPDT-like control. But C-RET already investigated certified early exit and contextual tree did not beat fixed ordering; a new formulation needs a clearly distinct mechanism and fresh validation before code.
5. If there is no convincingly distinct, testable habilitation tree extension, **finalize the original UNPUBLISHED RLR with fair real kTSP10/15/20 and matched T-derived capacity**, rather than fabricate algorithmic novelty. Exact original ARFF/probe parity and independent cohort confirmation remain blockers for strong publication claims.

**Implementation and audit:** `R/experimental/rlr_weighted_pairpair_tree_core_v4.R`, `test_rlr_weighted_pairpair_tree_v4.R`, `rlr_weighted_pairpair_first11_DV_v4.R`, `scripts/rlr_weighted_pairpair_first11_DV_audit_v4.py`. Original locked TeX and original RLR training code untouched. [Predeclared protocol](https://github.com/mczajkowskipb/relative-logistic-regression/blob/analysis/rlr-weighted-pairpair-tree-v4/docs/RLR_WT_V4_TWO_PAIR_WEIGHTED_TREE_FROZEN_DV_PROTOCOL.md). This is a **negative science gate**, not a CI failure.
