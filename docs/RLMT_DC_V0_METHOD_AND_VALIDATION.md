# RLMT-DC v0 — mathematical specification and validation contract

**Research branch** `feature/rlmt-dc-v0`; **method identifier** `RLMT-DC-v0`
(Relational Logistic Model Tree — Discovery/Confirmation).

**Status:** new experimental algorithm, not yet scientifically verified,
not authorized to replace published RLR or frozen RLMT-v1. Historical
GATE3/GATE4 outcomes and the 2026 publication `main` are immutable.

## 1. Motivation, no retrospective optimization

GATE3 (Chen/Singh) with fixed top1000 split candidates produced one
accepted split on Singh, but outer logloss was worse than frozen RLR
(0.51731647 vs 0.29002105); two minority-class test cases in its
low-support leaf dominated the harm. In GATE4, all 5000 supervised
top1000 pairs from five Gordon folds met the candidate child-size gate,
but **none** had >=6 of each class on both sides. The existing RLMT
uses pairs selected for class discrimination **also** as root split
proposals, which can conflict with class-mixed, locally estimable leaves.

**Research hypothesis H1** (not a finding): a bounded, label-agnostic
proposal universe of unordered feature comparisons, *filtered for
class support and discovery-sample residual heterogeneity*, will find
admissible split proposals that top-ranked class-discriminative
pairs cannot. This does **not** imply improved prediction.

**Research hypothesis H2** (not a finding): using a *disjoint*
confirmation partition to decide whether a single nominated split
improves held-out logloss reduces candidate-selection optimism
relative to assessing and picking among the same screened candidates
on their discovery data. This does **not** establish unbiased
test-set improvement or statistical significance.

The two ideas must be assessed separately in later ablations.

## 2. Object, fixed choices and notation

Given **outer-training** data `T={(x_i,y_i)}`, binary
`y_i∈{0,1}`, and completely untouched outer test `U`:

1. Stratify `T` once into **discovery** `D` and **confirmation**
   `V`, with one validation fold of three (`seed=42`; sizes
   about 2:1). All pool selection, predictions and residual screening
   used to nominate the split are computed **only on D**.
2. Leaf prediction uses the previously defined locally re-ranked RLR:
   K=10, m=10, d=2, top-1000 globally selected ordered pairs,
   elastic-net alpha=0.5, lambda.min, stratified leaf `cv.glmnet`.
   No switch to a new base learner.
3. `P_pred(D)` is the **supervised** top-1000 global RLR pair pool
   computed **only on D**, never on all T at discovery stage.
   For the final refit, `P_pred(T)` is computed **only on T**,
   as in the frozen RLR reference.
4. Draw at most `M=4096` **unordered, distinct, random gene pairs**
   from the predictor names on D using a deterministic seed.
   These are potential **split features**, not the RLR predictor
   shortlist. Sampling uses only feature names and RNG, **never y or
   outer test**. Do not enumerate or store all p² relations.
   Orientation is canonical by feature-name order; reversing the
   orientation merely exchanges child identifiers and does not change
   the two-region partition.
5. Require every proposed child on D to have
   `n_child>=16` and **at least 6 of each class**, unchanged
   from frozen RLMT. No class-support relaxation.
6. If no feasible proposal exists, use exact frozen RLR trained on T,
   without examining V labels.
7. On D, compute 3-fold out-of-fold parent RLR probabilities
   `p_i^{(-k)}` (D-only supervised predictor pool) and residuals
   `r_i=y_i-p_i^{(-k)}`. Nominate **one** feasible pair `s*`
   maximizing
   `S(s)=sqrt(n_left*n_right/n_D) * |mean_D(r|s=1)-mean_D(r|s=0)|`.
   Ties break deterministically by candidate draw order.
   This score is a proposal heuristic, **not a CV improvement claim**.
8. Fit on D a parent RLR and two child RLR models (using
   `P_pred(D)`, re-ranked locally in children) for **only s***.
   These fits **do not see V labels**. If fitting fails,
   choose no split and save the failure reason.
9. For each `i∈V`, compute the **paired** loss improvement
   `d_i=ell(y_i,p_parent,D(x_i))-ell(y_i,p_children,D,s*(x_i))`,
   where
   `ell(y,p)=-y log(clip(p))-(1-y)log(1-clip(p))`,
   `clip(p)∈[1e-15,1-1e-15]`.
   Let `g=mean_V(d_i)`,
   `se=sd_V(d_i)/sqrt(n_V)`.
   Accept the single split **iff** `n_V>=20`,
   `g>0`, `g-se>0` and all probabilities are finite/in [0,1].
   Exactly one candidate may ever be tested on V for a
   single outer-training run. **Do not** iterate over candidates
   using V feedback.
10. If accepted, refit both localized models on **all T**, using
    the outer-training global pair pool `P_pred(T)` and keep the
    chosen split orientation. Before training final leaves, verify
    unchanged leaf-support constraints on all T. If either fit
    fails, classify the run as failed rather than silently
    trying another candidate.
11. If not accepted, use the existing `rlmt_root_output()`
    path, which invokes the *unchanged* RLR implementation.
    `max_leaves=1` must reproduce RLR exactly. If accepted, export
    explicit `tree.json`, leaf artifacts, routes, predictions,
    and metrics with the new method identifier.
12. `U` must be accessed **only after** discovery, confirmation
    and final refit are frozen. `U` is *never* used for ranking,
    support checking, calibration, pruning or fallback.

This procedure is **discovery/confirmation sample-splitting**, not
the stricter 'honest leaf estimation' of Athey--Imbens, because
accepted leaves are finally refitted on all T, including D and V.
Do not label it 'statistically honest trees' or quote a coverage theorem
for classification logloss. One internal confirmation score is
still subject to variance and the subsequent full-T refit changes
the learning sample.

## 3. Partition safety invariants

- Every outer sample belongs to exactly one of D or V.
  `D ∩ V=∅`; `D ∪ V=T`.
- No functions producing the proposed pair are passed `y_V`.
  Unlabeled predictor names are acceptable; validation features
  are used only for validation predictions, not proposal creation.
- `P_pred(D)` is recomputed **inside D**; using `P_pred(T)` to
  fit discovery or confirmation models is forbidden.
- Proposal feature generator is independent of the labels.
- `V` selects no hyperparameters, candidate pool sizes or second
  candidate after rejection. Confirmation is a fixed 1-SE heuristic;
  no inferential/p-value claim.
- Complete code SHA, synthetic smoke outputs, RNG seeds, source
  cohort provenance and train/test assignments are recorded.
- No analysis, rerun or retrofitting on Singh/Gordon may be presented
  as independent evidence for RLMT-DC.

## 4. Stage gates and validation

### DEV gate (synthetic only)

- Same seed/predictor names -> identical proposed pair set;
  unique unordered pairs, valid endpoints, no p² allocation.
- Candidate generation insensitive to label permutation.
- Correct deterministic discovery/confirmation split; no index overlap.
- Class support reject/pass corner cases and ties in proposal score.
- Confirmation loss, SE and decision checked with hand-computed vectors.
- Canary: flip/permutate y only on V => nominated D candidate unchanged.
- Root L1 parity regression (exact RLR class, predictions and
  coefficients) on **synthetic** data, plus a no-split recovery test.
- Saved accepted leaf predictions replay from RDS. A positive
  synthetic split demonstration is a functional test, not scientific
  performance evidence.

### SCI gate (requires separate pre-outcome approval)

Before running any RLMT-DC scientific benchmark, register
the cohort list, licenses, independent subject IDs and cohort
overlap, raw checksum and class definitions; ensure no patient
appears across training and evaluation splits. Seek **new
independent omics cohorts**, not repeated splits of the already
examined Singh/Gordon matrices. Use a multi-dataset, grouped
or external-cohort evaluation where possible; the dataset, not
each fold, is the inferential unit. Fix complete code SHA, exact
algorithm versions, seed schedules, environment and a failure policy.

Comparators: frozen published RLR; unchanged RLMT-v1;
RLMT-DC-v0. Fixed ablations: original discriminative proposal pool
versus label-agnostic proposal pool; same discovery/confirmation
holdout schedule. Do not select the best ablation on outer test.
Primary per-dataset paired out-of-sample logloss difference,
secondary Brier, AUC, F1, accepted-split rate, retained
minority leaf counts, worst-case per-sample logloss and
training failures. Report unsuccessful/no-split outcomes.
Nested selection belongs inside each outer-training partition.

**No** new external cohort is currently registered or approved by
this design. Do not claim the methodology proves a predictive advantage
or constitutes a new published contribution without an independent
literature novelty review.

## 5. Literature and novelty boundaries

This is a new **proposed combination** in the RLR context,
not a claim that residual-guided splits, sample splitting or
validation based acceptance are generally novel:

- Ambroise C., McLachlan G.J. (2002). Selection bias in gene
  extraction on the basis of microarray gene-expression data.
  PNAS 99(10):6562-6566. DOI 10.1073/pnas.102102699.
  Motivation: reselect supervised features inside a training fold.
- Athey S., Imbens G.W. (2016). Recursive partitioning for
  heterogeneous causal effects. PNAS 113(27):7353-7360.
  DOI 10.1073/pnas.1510489113.
  Context for sample splitting; its formal causal-inference
  guarantees **do not transfer** to this logistic tree.

- Landwehr N., Hall M.A., Frank E. (2005). Logistic Model Trees.
  Machine Learning 59:161–205. DOI 10.1007/s10994-005-0466-3.
  Prior art for logistic models in decision-tree leaves.
- Zeileis A., Hothorn T., Hornik K. (2008). Model-Based Recursive
  Partitioning. JCGS 17(2):492–514.
  DOI 10.1198/106186008X319331.
  Prior art for score-/instability-guided model-tree splitting.
- Schlosser L., Hothorn T., Zeileis A. (2019). The Power of Unbiased
  Recursive Partitioning: A Unifying View of CTree, MOB, and GUIDE.
  arXiv:1906.10179. Prior work on the importance of model-score
  and residual choices for split-variable screening.

A future publication must benchmark against suitable **nonrelational
logistic model trees/MOB baselines** or explicitly narrow its scientific
claims to the *relational representation* and small-sample behavior.
Residual screening and logistic leaves alone are well-established;
neither is presented here as a novel invention.
