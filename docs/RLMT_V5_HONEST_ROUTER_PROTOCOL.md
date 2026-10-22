# RLMT-v5 — a safe, interpretable outer-train-only router (2026-10-23)

STATUS: **NEW PROTOTYPE SAVED; NOT EXECUTED ON UBUNTU.**
The goal is not to guarantee improvement on known heldout labels.
It is to prevent unreliable corrections from being applied everywhere.

## The complete preceding result that motivates this NEW experiment

19 previously EXPOSED CompCancer/biomedical cohorts ×5 frozen
original seed42 folds. v4 ran 95 complete outer-training
bundles, 475 method-fold outcome vectors, and independently
audited **PASS** with exact original-paper RLR patient-level
parity. Exit code 0. Equal-dataset mean heldout logloss:

| Variant | Actual route folds/95 | Dataset wins/losses vs ORIGINAL_RLR | Mean ΔLL vs ORIGINAL_RLR | Mean ΔLL vs GLOBAL_CAL |
|---|---:|---:|---:|---:|
| ORIGINAL_RLR | 0 | 0/0 | +0.000000000 | −0.001301927 |
| GLOBAL_CAL | 0 | 11/8 | +0.001301927 | 0 |
| CONF_ROUTE | 91 | 10/9 | +0.001615592 | +0.000313665 |
| ACTIVITY_ROUTE | 89 | 11/8 | +0.001233117 | −0.000068810 |
| SOFT_GAP_ROUTE | 91 | 7/12 | **+0.050673367** | **+0.049371440** |

No v4 method demonstrated better mean logloss than ORIGINAL RLR.
The earlier 2-cohort v4 SOFT_GAP mean −0.009452910
did **NOT** replicate and must not be presented as evidence
of superiority. All 19 datasets have influenced model
development, and all future trials on them remain exposed.
SOFT_GAP is excluded here due to retrospective evidence
of unstable large negative corrections, with this
selection disclosed.

## Precise v5 hypothesis

Treat "no change to original RLR" as the default.
A second-stage ROUTER selects whether a small
relational correction has *internally* demonstrated benefit
on patients that were not used to train the proposing
model. If not: exact saved original RLR probabilities.

For every outer training fold independently:

1. Make fixed 3 inner stratified folds, seed=1403.
   All 19 previous outer test samples remain unseen.
2. For each inner heldout fold, fit the FULL frozen
   \`rlmt4_fit\` on the other 2/3 of outer training.
   Its root uses original manuscript RLR;
   feature nomination, training-only threefold original
   RLR parent OOF, β and γ are completely re-fit
   **without inner heldout labels**.
3. Obtain honest heldout probabilities from four
   candidates: original RLR (default),
   GLOBAL_CAL, CONF_ROUTE, ACTIVITY_ROUTE.
4. For each of three corrections test fixed convex
   original/correction mixes with weights
   {0.25,0.50,0.75,1}. Paired 3-fold mean gain:
      Δ=logloss(original)-logloss(mixture).
   Score = mean fold gain - 0.5×SE(fold gains)
           - 0.00025 complexity penalty.
   Accept a correction only when score>0.00025.
   Otherwise weight 0 and *exact* original RLR.
5. Refit frozen \`rlmt4_fit\` on the **entire outer
   training**. The selected class/weight is frozen
   BEFORE prediction on outer heldout. New predictor
   p_final=(1-w)p_ORIGINAL_RLR+w p_corrector.
6. Critically, w=0 is an EXACT nested original-RLR
   predictor and there is NEVER a required split.
   When minority-class support is too small for
   honest nested model fits, fall back to original.

The 0.00025 penalty, 0.50×SE penalty and four blending
weights are novel development choices after exposure of
19 cohorts. They are **not optimized against v5 results**;
changing them after looking at these test labels must
create a NEW research version and will not be
independent validation.

**Interpretability:** at most the published RLR's ten
selected signed gene-order pairs, plus one named
relational correction pair, 2 low-dimensional
logistic coefficients (β,γ) and one scalar correction
weight. No local full RLR refit on tiny leaves.
Although the term model tree is used for research
continuity, this is an abstaining/soft-repair router,
NOT a classic hard decision tree.

## Experiment and scientific gates

First: real synthetic parent/correction training,
true disjoint inner CV, RDS exact replay, artificial
positive and negative train-only router unit gates.
Second: six previously exposed outer fold positions
on chen-2002 and gordon-2002, with frozen archived
patient train/test ID parity and source SHA256 checks,
and exact original-paper RLR probability parity.
Third: independent Python audit reruns heldout
patient LL/Brier, validates original RLR parity, route
identity, whether the final decision was made before
outer heldout predictions. No overwrite of original RLR.

Only if the pilot is technically sound and
promising (without post-hoc tuning), run full fixed
19×5 audit, including all 19 datasets equally weighted,
balanced accuracy/Macro-F1, LODO, selection frequencies,
actual route usage and internal-gain predictive reliability.

**No guarantees:** strong original RLR may be Bayes-near-optimal
for some cohorts. This router can still overfit meta-selection
to small inner folds. Do NOT claim reliable improvement until
new, independently acquired patient/lab-batch/cohort validation
with the algorithm locked before any outcomes is completed.
