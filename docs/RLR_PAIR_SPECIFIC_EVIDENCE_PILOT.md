# RLR rank-evidence extension: publication-dataset component pilot

Scientific role: **exploratory real-data ablation on 2 of the original article's 11 CompCancer binary-classification datasets**; not additional synthetic cohorts, not a new set of clinical studies, not an independent fresh validation. The original unpublished RLR paper remains the primary publication.

## Why this rather than new trees

The original paper's central signal is improved probability quality over hard-voting k-TSP, notably logloss/Brier/ECE. Already published RLR defines within-patient signed gene-rank gaps `Δ_ij = rank(j) - rank(i)`, encoding each relation `-1/0/+1` under a single margin m=10, then regularized logistic regression on the top K selected gene-disjoint pairs. The existing supplementary ablation already varies the **one shared** m∈{1,5,10,20,50} and K∈{10,20,50}; therefore repeating global m is not a new result.

Earlier CRIT/RLMT experiments (full RLR in terminal leaves and even 2/5 local predictors) did not yield reproducible V-accepted improvement. Multirelation RLR in root nodes was similarly not confirmed under the tested small-sample constraints. Stop that family; do not keep retuning tree hyperparameters on exposed E outcomes.

## One narrowly controlled question

Does the **pair-specific signed rank separation** carry useful *graded predictive information* beyond the frozen `±1/0` step function, while preserving within-sample rank invariance and an interpretable original-RLR sparse coefficient model?

Three predeclared variants, **same selected pairs and identical elastic-net fitting settings**, altering ONLY encoding:

- `TERNARY_M10` = original exact published encoding: for each signed gap Δ, +1 when Δ≥10, −1 when Δ≤−10, zero otherwise.
- `ADAPTIVE_Q25` = per-pair threshold `m_j = max(1,round(Q_.25^train(|Δ_j|)))`, learned without labels from training data for that pair. Each sample gets sign(Δ) if |Δ|≥m_j, else zero. This creates pair-adaptive inactivity without V/E access.
- `GRADED_M10_Q75` = preserve original zero region |Δ|<10 and encode `sign(Δ)*clip((|Δ|-10+1)/s_j,0,1)`, with `s_j = max(1,Q_.75^train(max(|Δ_j|-10+1,0)))` estimated without labels. This graded rank-strength proxy is *not a calibrated measurement-uncertainty estimate*.

Inputs restricted to `chowdary-2006` (104,182) and `singh-2002` (102,339) from original CompCancer benchmark. Both are **previously analyzed developmental cohorts**. Predeclared three stratified 3-fold outer splits seed42 for each cohort. For each outer split calculate a single original switchBox supervised top1000 pair-pool **using outer training only**, then original strict `d=2` greedy shortlist K=10,20,50; frozen `alpha=.5`, `lambda.min`, CV deviance and seed42. Train-specific encoding Q values are never fit on E. Exactly 54 matched fits across 2 datasets×3 folds×3K×3 encodings, with patient-level E logloss/Brier, active relational structure, and probability archives. Independent Python audit recomputes all patient E scores, verifies no repeated E patients across folds, base ternary parity and printed Δlogloss.

## Stop criteria

- If no coherent, reproducible heldout improvement even in this small **exploratory** pilot, **stop**. No other arbitrary sigmoid/margin tweak sweep on the same exposed E data.
- If there is a plausible signal, evaluate the one chosen frozen variant on **all 11 already-in-the-paper CompCancer datasets** using their *existing* original 5×10 / 5×2 benchmark splits, and re-evaluate PTB fixed and LOCO if claimed. Do not introduce new datasets or silently extend protocol.
- Rank-gap magnitude or adaptive pair thresholds alone are not novel inventions (literature already considers rank-difference weighting and pair ranking). The scholarly contribution requires a justified integrated formulation, a controlled component ablation, and meaningful incremental evidence. No claims of new trees, novel XOR, robustness to unseen cohorts or clinical utility.

User constraint: work only within datasets in original unpublished RLR paper; no further RLMT/CRIT decision-tree expansion.
