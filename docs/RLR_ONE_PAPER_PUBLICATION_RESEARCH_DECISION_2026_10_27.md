# RLR paper integration and decision-tree research target — 2026-10-27

## Non-negotiable research objective

**One unpublished paper, one integrated story**: original Relative Logistic Regression (RLR) plus one genuinely improved, interpretable within-sample rank-based extension. The original RLR publication lock remains untouched until a new model passes scientific gates. Use original CompCancer **11 identical dataset IDs** and original **470 patient-wise split specifications** for final comparison, plus the original PTB and historical analyses as applicable. The actual original processed **ARFF checksums must be verified** on the research server before claiming an *identical* numerical reproduction. Reconstructed public Schliep matrices are valid exploration but not paper-original artifacts.

We are NOT requiring a tree to beat original RLR on every metric. Original RLR already beats classical k-TSP substantially on probability quality (CompCancer RLR50 logloss .226 vs kTSP50 .378), but has less compelling discriminative means (RLR50 AUC .957/F1 .869 vs kTSP50 AUC .963/F1 .884). The publication target is to preserve RLR's major logloss/Brier gain while closing the F1/AUC gap. Strong value-based LR is a secondary benchmark, not a necessary win.

RLR and its tree extension can be reported together, with a strict distinction between their predictive and structural contributions. A tree that only changes logloss by a negligible amount and has no benefit over equal-feature no-tree controls is NOT a separate scientific contribution.

## Pilot matched benchmark (not the full 11/470)

### Actual GitHub run, identical 6 development folds and all real baselines

- [Actions #40473159026](https://github.com/mczajkowskipb/relative-logistic-regression/actions/runs/40473159026) **SUCCESS**.
- Cohorts: previously inspected original-paper identifiers Chowdary-2006 (104 patients) and Singh-2002 (102), public source-reconstructed matrices.
- Actual R original SWAP.KTSP trained at k=10/20/50 (not a manually simulated equal-vote proxy), TSP1, original RLR10/20/50, original value LR from train-only MAD top1000, and greedy anchor-gated conditional RLR with both original K20 and K50 offsets. A 3-fold outer split × one repetition each, 6 same E folds.
- Ten aligned arms ×206 original patients = **2060 heldout patient-method probabilities**, 60 method-fold metrics; independently reproduced AUC, **both positive-class F1 and Macro-F1**, BA, logloss and Brier in independent Python; native SWAP.KTSP class labels retained (k-TSP F1 must not be imputed from vote probability threshold). Patient coverage and T/E disjointness verified. **These six folds are NOT fresh post-selection confirmation.**

| Method | Mean AUC | Mean Macro-F1 | Mean logloss |
|---|---:|---:|---:|
| TSP1 | .89014 | .88918 | 3.71629 |
| kTSP10 | .97065 | .93534 | .34995 |
| kTSP20 | .95871 | .92551 | .38207 |
| kTSP50 | .95050 | .92007 | .28129 |
| RLR10 | .97074 | .93553 | .24137 |
| RLR20 | .97131 | .92867 | .31777 |
| RLR50 | .96959 | .94002 | .25337 |
| Greedy anchor + two child RLR20 corrections | .97240 | .92861 | .30742 |
| **Greedy anchor + two child RLR50 corrections** | **.96959** | **.94002** | **.24705** |
| Value LR | .95482 | .95558 | .33387 |

Note that in the small pilot **RLR50 already exceeds kTSP50 on all three summary metrics**, while the tree improves RLR50 logloss by .00632 but does NOT improve AUC or Macro-F1. LR wins on Macro-F1 .95558. These numbers are **not** interchangeable with original locked 11-cohort CompCancer aggregates; in the original CompCancer, RLR50 is worse than kTSP50 on both AUC and F1.

## Evidence from previous tree searches — decision

The existing many negative experiments matter:
- Per-pair adaptive rank margin, graded margin, RLMT/CRIT interaction corrections and sparse-leaf logistic tree variants mostly did not beat RLR, often failing on small node support.
- Group-anchor majority and ordinal 0..K vote thresholds were mostly redundant with a single reference TSP pair: **95.12% average train-only agreement** with best fixed median reference across six folds; greedy class-proportion-leaf stumps had no consistent advantage and failed to induce real depth2.
- First OOF residual-driven **lookahead root screening** ranked six panels constructed from globally selected TSP hubs, and required TWO working child RLR fits; valid positive child fits in only 2/6. The one V-accepted lookahead tree worsened E. Looking ahead by **reranking frozen global anchor panels** is not actual reference-optimization.
- A *greedy root + two child RLR* showed one encouraging but exploratory loss signal: `ΔLL=-.010343`, 4/6 E folds improved (K20 parent) versus global same-pair control `-.001335`, and gate-only `-.000562`. Mean BA not improved. Equivalent RLR50-offset greedy root gave `ΔLL=-.00632` on six E folds, no Macro-F1 change.

There is no evidence yet for a new publishable method from tree extensions. All reported results remain development-stage and should NOT be massaged with further small-cohort threshold scans.

## Strongest next scientific candidate: Selective Anchor-Gated RLR (SA-RLR)

Consider a compact **one-root model tree** with an original global RLR50 plus a sparse *selective* relational correction on ONE branch, falling back to RLR50 on the other, or TWO only when both are supported.

For patient x, let `z_g(x)=1[sum_{a∈A}1[x_g>x_a] >= t]`, a count/ordinal-within-sample gate. Conditional model:

`logit P(Y=1|x)=η_RLR50(x)+1[z_g(x)=0]Δ_0(x)+1[z_g(x)=1]Δ_1(x)`

where `Δ_b` are at most 3–8 branch-specific signed ternary m10 relation features with strong coefficient shrinkage. **Allow exactly one viable branch**: if Δ_other lacks stable informative signal, set it to 0. Forcing 2 successful leaf models discards useful asymmetric disease subtypes and was a major feasibility failure of lookahead V1. Keep exact parent offset as the safe fallback; model is within-sample invariant (assuming fixed probe universe).

Do not nominate local corrections by globally discriminative TSP rank alone. Instead use OUT-OF-FOLD original-RLR residuals to score **branch-conditional improvement**. Strong controls must use the SAME nominated gene-pair features as global RLR, and the root-only correction without branch-dependent interactions, plus no-change original RLR. Use Bayesian/shrinkage intuition to regularize intercept changes; avoid independent unpenalized logistic leaves.

### Truly different look-ahead anchor optimization (if above can exhibit train-only signal)

Optimize the anchor SET, not only select among globally optimized TSP anchor panels. Given a hub gene g, exchange one reference a for a candidate a'; maintain an inner beam of width ~5, <= 2 swap steps and 3/5 reference anchors. Use a cheap **crossfitted gradient/Hessian conditional feature score** as stage one; score finalists against **inner heldout child RLR logloss**, with penalties for support imbalance, duplicate anchor votes, unstable reference rank positions and number of unique genes. Cache within-patient ranks and all pair indicators. This gives a bounded finite search rather than naive exponential combinations or premature EA. A single branch may remain uncorrected, preserving parent RLR. Only if such a method reliably finds valid nontrivial conditional experts across datasets should EA-driven global structure be considered.

### The most neglected, high-value discrimination diagnostic

**Class-label threshold must be audited, independently of probability estimation.** RLR logloss can be better while its 0.5 decision F1 is below k-TSP. Train-only, cross-fitted threshold selection for positive-class F1 or Macro-F1 may close this gap *without changing rank predictors, AUC, calibrated probabilities or logloss*. This is not in itself novel; report it as a fair decision-rule comparator and distinguish it from a tree-specific effect. Never use an outer test fold to fit the threshold. Train-only class-prevalence threshold is a simple fixed comparator. Ensure k-TSP native classifier label F1 and RLR thresholded F1 are both separately reported, preserving their manuscript conventions.

## Final-paper comparator requirements, no substitution

1. Classic TSP1 and **actual** switchBox k-TSP 10/20/50, with their class labels and vote scores as in original manuscript; also calibrated kTSP Platt from the original supplementary protocol.
2. Locked original RLR10/20/50, and candidate **RLR50 + selective tree**, same exact outer patient fold/feature schema/positive label. At minimum K20/50 for the tree may be useful as a sensitivity.
3. Real value-based LR with training-only MAD1000; acknowledge that only relational models enjoy within-sample monotone invariance.
4. Relevant 2024+ ordinal-tree comparator: `ranktreeEnsemble` rank-based boosted/random forests (Lu et al., *Journal of Translational Medicine* 2024); if feasible compare to original author's **actual CUDA TSPDT**, and TSPG/top-scoring pair of groups. Do not present greedy TSPDT proxy as the original TSPDT.
5. Primary statistical unit is DATASET (11), not 470 repeated folds; report F1 (positive-class as manuscript), Macro-F1, AUC, logloss, Brier, calibration, runtime, active rank relations/patient and compactness (number of nodes and genes).
6. Explicitly examine differences and the paired CI, no multiple-testing fishing, no posthoc selection of K based on E. Use stronger per-dataset validation before optimistic global claims. Preserve total method effort fairness.

## Stop/go before editing manuscript

**GO** only if a frozen SA-RLR or other new method shows a consistent, quantitatively meaningful improvement in at least one predeclared manuscript-critical performance metric (preferably classification F1 or AUC vs kTSP and an honest structural benefit) on the exact 11 original study cohorts under the original 470-fold protocol, with logloss not materially degraded relative to RLR50, and stronger interpretability than established ranktreeEnsemble/TSPDT alternatives.

If improvements are absent or uncertain: **publish the already completed RLR manuscript as the original standalone contribution**; negative supplementary tree experiments can remain separate. Do not hold up an already prepared paper indefinitely waiting for a positive tree result.

## Relevant novelty literature

- Tan et al. 2005, Bioinformatics, k-TSP multiple disjoint pair hard-vote classifier.
- Xu et al. 2007, TSPG (Top-Scoring Pair of Groups), prevents claiming multi-anchor majority as unprecedented.
- Lu, Yin & Chen, 2024, *J Transl Med* 22:140, **ensemble rank-based trees for single-sample gene-expression classification**; boosted and forest rank trees reportedly outperform kTSP on 12 cancer transcriptomic datasets. An ordinary rank-tree ensemble cannot be marketed as novel.
- Existing TSPDT (user's 2011), ARXA (ICCS 2020), and evolutionary relational multi-test trees from the user's prior work require an explicit technical comparison and clear differentiating contribution.

**No claimed journal acceptance or publication-quality new tree efficacy exists yet.** The current reproducible deliverables are implemented benchmark code, real CI, independent patient-level audit and a research plan grounded in the original locked manuscript.
