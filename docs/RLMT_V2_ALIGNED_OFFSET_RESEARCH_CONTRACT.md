# RLMT-v2: aligned parent estimation and conservative relational leaf offsets

**NEW STUDY DESIGN, NOT CODE OR RESULTS. 2026-10-21.**

Research branch: `feature/rlmt-aligned-probability-v2`. This study is
*separate* from published RLR, RLMT-v1, RLMT-DC-v0,
RLMT-CA-PS-v1, and historically faithful TSPDT/Weight-kTSP.

## 0. Evidence boundary and decision

The 19×5 publication cohorts are EXPOSED method-development data.
Their complete 285 archived RLR/RLMT results and 95 CA-PS
results have been inspected at the aggregate and accepted-tree levels.
Consequently **neither holdout tuning nor a claim of independent
confirmation on these 19 cohorts is scientifically defensible**.

Frozen paper `main=4079a5ac400f1d2cfded523c7000f655470925ce`.
Frozen original DC code `4873878f57f806bf5b52144756f8087b17fa540c`.
Frozen completed CA-PS 95-fit source
`e049fe20530ae7c5378db73f4e9186b5f83f789e`.
CA-PS 13 accepted trees, 26 leaves, **4 benefit / 9 harm** by
outer heldout logloss; dataset-equal delta LL vs RLR
**+0.006177806**. DC-v0 **+0.017217934**.
**No acceptance gate has established CA-PS superiority over RLR.**

Interpretation is narrower than a posthoc causal diagnosis: the
0.5-threshold classification accuracy need not follow probability
quality. Some outer fold losses are adverse despite unchanged
accuracy; underdetermined sparse local leaf probability estimates
and D/T predictor misalignment are plausible, unproven causes.

**CA-PS-v1 status = STOP AND FREEZE.** No per-dataset
threshold changes, no ex-post adjustments to g, B, local
model regularization, splits or benchmark composition.
The new research program must use a new branch and model ID.

## 1. Research hypotheses: two *separable* interventions

### H1. Algorithmic alignment of train-only confirmation and full fit

Current DC/CA candidate confirmation on independent training partition
V predicts with parent/leaf models trained on D, using a D-only
supervised relation pool and deterministic CV. If accepted, the
outer-test predictor is subsequently fit on the full T using
T's separately calculated relation pool and, for the original RLR
root, a different CV-fold implementation.

The D-trained and T-trained *models* can never be identical; this
is not the proposed invariant. The needed invariant is that both
models are instantiated by **exactly the same learning procedure**
A(S), up to using different training sets S.

Define one reusable `fit_parent_RLR(S)` with:
- original supervised top1000 pair pool computed ONLY on S;
- original deterministic K=10,m=10,d=2,alpha=0.5 model;
- same fixed CV scheme, fold assignment rule, lambda selection
  and pair encoding for S=D and S=T;
- saved pool ordering, fold IDs, lambda, coefficients, fitted
  estimator metadata; no implicit dependence on output path;
- original RLR published definition remains unmodified
  in \`main\`, and its archived scores stay a comparator.

Likewise define `fit_candidate_leaf(S, nominated_route)`
using the identical S-only reranking, inner fitting and prediction
procedure when S=D (validation candidate) and S=T
(final refit). Final T may use the historical RLR parent
code if parity has been **demonstrated** on training-only
simulations. If it cannot, use a separately identified
`RLR_ALIGNED` control rather than silently
relabeling a different parent as original RLR.

H1 ablation: retain original DC candidate score,
class support, 1SE confirmation rule, and the old
local leaf estimator, changing **only** the alignment
of full and validation learner procedures. This isolates
whether estimator inconsistency is materially related
to unstable confirmation. Alignment is a process
property, **not a guarantee of improved expected loss**.

### H2. A conservative **log-odds correction**, not a free local fit

At a candidate leaf l, let p_parent(x) be the consistently
fitted global RLR probability. Define

\`\`\`
eta_parent(x) = logit(clip(p_parent(x), epsilon, 1-epsilon))
p_new(x) = sigmoid(eta_parent(x) + delta_l(x)).
\`\`\`

The correction is trained **only from independent parent out-of-fold
(OOF) scores on leaf-training observations**, never a full
T in-sample parent score reused as a supposedly honest
offset. This is mandatory to avoid local target leakage.

Candidate implementation uses the already available
`glmnet::glmnet` binomial offset / newoffset functionality,
not a custom logistic optimizer. Penalize all correction
coefficients; no unpenalized free leaf intercept that
can force near-zero or near-one probabilities from a
tiny subgroup. The new leaf model should contain
training-only selected relational features using original
`encode_pair_matrix` and at most K=10 pairs.

**Numerical rule:** any probability clipping/offset
epsilon is a computational safeguard, not an
after-the-fact calibration tuned to previous
breast2 outcomes. Exact fixed epsilon and lambda
selection must be frozen *before* new cohort outcomes.

**Critical leakage hazard:** tuning glmnet lambda via
ordinary internal leaf CV on parent OOF offsets
can reintroduce dependencies (one held-out leaf sample
may influence another's OOF parent fit). The
implementation must either:

1. perform **nested** parent OOF reconstruction
   inside each leaf CV training fold; OR
2. avoid internal leaf-outcome-tuned lambda altogether
   with a defensible, prospectively fixed penalty and
   no data-driven external adjustment.

Pick and document ONE before implementation.
Do **not** treat the existing CA convex probability
mixture as also active in H2: to test H2, prediction
should be the offset corrected parent, not a posthoc
hybrid of two rescue devices.

If delta(x)=0 identically, the offset estimator must
give exactly the parent probability (within declared
numerical tolerance), whether or not a split is
proposed. If no split is accepted, return the **exact**
parent predictor, not a refitted alternate root.

H2 is **not identical** to CA-PS-v1's convex mixture,
which did not penalize leaf logit corrections against
parent odds.

## 2. Frozen research sequence: no mixing mechanisms

**A — scientific feasibility / implementation gate:**

- Validate original RLR parent parity (or register a
  separately named aligned parent baseline).
- Define the consistent fit/predict interface and
  replay all deterministic RLR/leaf artifacts.
- Use *synthetic only* nonlinear relational data to
  ensure accepted/rejected branches, all class-support
  rules and leakage invariance can be exercised.
- Record R warnings (not only final exit code).
- Do not run 19 old datasets to choose design constants.

**B — H1-only ablation:** Same DC split nomination and
free leaf predictor, but identical D/T fitting procedure.
Use a distinct model ID, fixed splits and *one* planned
independent test benchmark, conditional on A passing.

**C — H2 offset leaf:** Only after B is frozen, compare
H2 vs aligned parent and H1-only. Both nomination and
acceptance must be based solely on train-only discovery
and confirmation; do not adapt them based on external
test outcomes.

**D — mechanistic and external evaluation:** Decide in
advance if the primary application is
probability/risk prediction (logloss + calibration)
or class discrimination (balanced accuracy/F1).
Do not use test discrimination gains to reinterpret
bad logloss as a positive calibration result.

## 3. Development gate: what must be tested

1. **No leakage / source integrity:** candidate feature
   scores, class-support, supervised relation selection,
   OOF offsets and confirmation use no outer test labels.
2. **Discovery candidate freeze:** at most one winner,
   chosen entirely on D, with independent V acceptance.
3. **Full estimator identity:** the same A(S) helper is
   used on D and on T; no manual change from
   stratified fold IDs to implicit CV-fold generation.
4. **Counterfactual label permutation:** changing
   *only* outer-test outcomes leaves every fitted
   model, nominated split, accepted status, prediction
   and leaf weight unchanged.
5. **Original parent replay:** saved RDS, selected pairs,
   coefficient vector, fold IDs, lambda, per-sample
   probabilities; fail on mismatches.
6. **Conservative correction:** delta=0 returns parent,
   probabilities finite/in [0,1], no unpenalized
   correction offset drift; test small classes.
7. **Independent confirmation correctness:**
   full paired proper logloss, gate consistency,
   no post-V second candidate.
8. **Runtime and portability:** uses existing R
   switchBox/glmnet; works without cloud token/API
   charges; experimental results isolated from
   all older output roots.
9. **Human/clinical interpretation:** class label order,
   prevalence in TRAIN and chosen cohorts, reliability
   and clinical costs recorded; do not force
   50:50 child prevalence.
10. **Reproducible model loading:** full fit saved,
    reloading yields identical test probabilities.

## 4. Independent biomedical cohort data gate BEFORE model outcomes

Do not label the old 19 cohorts as external because
a new method name is used. New clinical/bulk-microarray
RNA expression comparisons must have:
- verifiable repository record + raw source and SHA256;
- lawful research availability and usage/license terms;
- binary clinically meaningful target and class counts;
- explicit sample identifier, patient/group
  identifier if multiple biopsies/time points,
  exclusion list, harmonized class labels;
- no patient overlap across train/test or cohort sites;
- a precommitted split/sampling/CV plan;
- all preprocessing, filtering and gene mapping
  learned within the outer train fold;
- fixed primary **unweighted** heldout logloss
  if probabilities are the primary model output;
- Brier + calibration intercept/slope and
  balanced accuracy/macro-F1 + class-specific recall;
- all datasets/splits reported, **including failures,
  empty accepted groups and adverse outliers**.

If external data are selected after inspecting
the previous 19 cohorts, that selection process
and any class-balance motivations are transparently
reported. Never choose a favorable dataset after
checking new test results.

## 5. STOP conditions and release criterion

Do not merge into published \`main\`, rename old
\`RLMT-DC\`/TSPDT algorithms, or claim discovery
of a superior classifier from this plan.

No v2 superiority claim until:
- mathematically self-consistent learner is implemented;
- synthetic gate passes with *real* two-leaf fits;
- new untouched biomedical cohort independently
  passes provenance/split-group checks;
- all frozen heldout outcomes are generated once,
  externally and without retrospective retuning;
- original RLR remains explicitly reported,
  alongside H1/H2 and failed alternatives.

The methodological contribution should focus on
interpretable **relational effects on log-odds under
proper-scoring evaluation**, not choosing a
classification threshold that hides badly
calibrated probabilities.

**Present status:** V2 contract is a design artifact.
No v2 R implementation, new dataset acquisition,
actual v2 model fit, external results, or accepted
scientific hypothesis are claimed by this commit.
