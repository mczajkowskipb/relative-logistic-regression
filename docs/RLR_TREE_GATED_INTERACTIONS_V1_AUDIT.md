# Ternary Tree-Gated RLR V1 — original-manuscript biological cohorts, pilot result

**Decision date 2026-10-26. Original unpublished RLR code and locked manuscript unchanged.** The research objective is to construct an innovative, high-performing interpretable *tree-structured* extension that can genuinely belong to the habilitation cycle on tree models for omics. The newly tested method is a **three-way tree-gated logistic interaction model**, not independent logistic models in tiny leaves.

## Methods and scientific controls

- Only `golub-1999-v1`, `shipp-2002-v1`, `west-2001`: all from the **11 original binary CompCancer datasets** in the RLR paper. These three were not used in the previous ARLET 3-cohort / C-RET 6-cohort pilots, but they are **not an independent external cohort** and publicly reconstructed filtered source inputs are **NOT yet checksum-equivalent to original processed ARFFs**.
- Three outer stratified folds, seed42, 3 patient-disjoint test folds per dataset; 6 methods compared within each fold.
- Original RLR K10/K50, original `fit_tsp1_manual_train_test`, original `fit_ktsp_train_test(k=9)` used exactly from unchanged `R/new/07_world_protocol_utils.R`. `TSPDT_GREEDY_PROXY` is an intentionally simplified two-level greedy relational stump tree, **NOT original user's CUDA TSPDT**.
- T-GRLR: global main `K=10` original rank-pair ternary representation `m=10`, `d=2`; three-way root gene-pair gate with outcomes `−1,0,+1`; each root branch can activate up to 3 additional m10 ranked pair relations. The original switchBox top1000 candidate pool is broadened with up to 550 seeded pairs spanning all assayed genes. A *single global elastic-net* logistic fit learns all main/root/conditional feature coefficients, not separate tiny-leaf experts; tree structure is an explicit hierarchical conditional interaction basis.
- All root and branch pair nomination is supervised on **inner D only**. Inner V is used to choose among up to 9 root proposals, with prespecified complexity penalty and `≥0.008` validation logloss gain to accept; if not accepted, fall back to literal RLR K10 (verified for parity). The final selected model refits once on whole outer T. Outer E is held out for scores.
- This proposal tests a different mechanism from prior unstable leaf-experts and node-offset trees: **globally regularized coefficients plus branch-gated relational interactions**. It is not automatically a novel decision tree; hierarchical interaction logistic models already exist in literature (e.g. Bien et al. 2013, doi:10.1214/13-AOS1096), and further prior-art work is mandatory.

## Heldout results — verified by independent patient-level audit of archived CSVs

Each entry = **mean over three 3-fold outer heldout split-level scores**. Logloss is lower-better; BA/macro-F1 higher-better. 3 cohorts (72, 77, 49 patients) × 6 methods; 1,188 archived patient×method rows. Full independent replay of BA, macro-F1, logloss, Brier (raw probabilities for Brier, clip only for logloss) and disjoint patient IDs PASS for all 3 datasets.

| Original paper cohort | Model | BA | Macro-F1 | Logloss | Brier |
|---|---|---:|---:|---:|---:|
| Golub | RLR K10 | .95023 | .95348 | .16435 | .0369 |
| Golub | RLR K50 | .93981 | .93940 | .15640 | .0424 |
| Golub | kTSP 9 | **.96065** | **.95588** | .17377 | .0429 |
| Golub | T-GRLR | .91319 | .92256 | **.11887** | .0369 |
| Shipp | RLR K10 | .75591 | .77857 | .45357 | .1200 |
| Shipp | RLR K50 | **.89916** | **.89916** | **.31199** | **.0642** |
| Shipp | kTSP 9 | .85368 | .80105 | .33590 | .1063 |
| Shipp | T-GRLR | .84361 | .85801 | .47307 | .0830 |
| West | RLR K10 | **.87731** | .87532 | .67406 | .1241 |
| West | RLR K50 | **.87731** | **.87665** | **.36755** | **.0943** |
| West | kTSP 9 | **.87731** | .87532 | .37776 | .1134 |
| West | T-GRLR | .85880 | .85558 | .41971 | .1205 |

**Interpretation:** T-GRLR provides one positive *probability-only* signal on Golub (.1189 vs RLR50 .1564 logloss), but reduces Golub BA from .9398/.9502 to .9132. On Shipp and West it underperforms original RLR50 in all BA/F1/logloss. It is also below kTSP9 in BA on all three (Shipp .8436 vs .8537, West .8588 vs .8773, Golub .9132 vs .9607). The three branches do yield local interpretability: the ternary **zero/inactive** root outcome occurred in **27%–40% of outer training observations** in Golub's selected root across folds; this is an empirically meaningful state, but is NOT proof of higher class accuracy. Across cohorts the accepted trees have 17–18 distinct global relational features and ~13 assessed along a patient route, not necessarily fewer than the *nonzero* features in the original RLR50, which would require direct measured comparison.

**Important validation observation:** the same predeclared inner V acceptance criterion admitted the tree on all 9 outer folds, yet external E discrimination was not consistently improved. This strongly suggests small-sample validation selection optimism: do not tune further to these observed E scores and claim independent validation.

**Original TSP1 raw 0/1 scores** yield misleading large logloss. Interpret classification BA/F1 comparisons normally, but do not draw calibrated probability claims from TSP1 logloss. Original CUDA TSPDT has not been executed and no TSPDT victory claim is supported.

## Status and audit repair

- Primary 3-cohort GitHub Actions workflow: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40374157979
- Fixed independent Python audit runner commit `d6edaf294a02e7895ca6117f72a924970cf1d57b` corrected a **pure audit metric bug**: first auditor clipped raw hard TSP1 probabilities for Brier while original R correctly clips only logloss. This caused original Shipp/West CI audit false failures; model fitting and patient archives were complete.
- Corrected run: https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40374362795 . Golub completed full GitHub audit PASS; at initial report write, Shipp/West GitHub jobs were still provisioning dependencies. Independently, archived **Shipp/West** output was re-audited offline against the corrected raw-probability formulas, including 3-fold E coverage and all four metrics: **PASS**. They should not be labelled GitHub-hosted PASS until the rerun actually reports success.
- Verified archive files in the working research environment: `rlr_tgr_golub_v1.zip`, `rlr_tgr_shipp_preaudit.zip`, `rlr_tgr_west_preaudit.zip`. Keep GitHub Actions artifact IDs for permanent retrieval.

## Methodological decision and research gate

**Do not insert T-GRLR V1 into the unpublished RLR paper as a confirmed better classifier.** The nested-gated interaction approach was scientifically valid to investigate, but did not beat strong RLR50 and kTSP9 on heldout BA. A stronger, journal-worthy tree would have to beat existing baselines in BA/macro-F1 or yield a genuinely substantial interpretability/certification benefit while maintaining comparably good prediction quality.

Potential mechanistic issue: conditional feature fitting on ~25–50 training cases is unstable even when *globally* penalized; a single 4-fold inner validation partition was optimistic on all nine folds. Mere enlargement of the pair universe and addition of EA do not address this data-limited variance.

Existing **C-RET** remains the only consistently positive *tree-like* result: formal certificate preserving original RLR decisions with 21–58% fewer evaluated rank relations across six original-paper cohorts, but fixed-order sequence is **not automatically a novel learned decision tree**. Existing conditional policy generated almost no advantage over fixed coefficient ordering. A genuinely novel tree contribution would need theoretical distinction and a clear reproducible improvement over both original RLR and simple sequential evaluation; not a cosmetic tree wrapper.

Scientific conclusion: preserve transparent negative findings. No test-driven hyperparameter sweep on the just-exposed 3-cohort E splits. No original manuscript editing until a stronger extension is supported.
