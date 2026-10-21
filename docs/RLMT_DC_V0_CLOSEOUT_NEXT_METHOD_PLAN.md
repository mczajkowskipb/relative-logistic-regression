# RLMT-DC next-method scientific design — post-v0 research brief

**Status: hypothesis/specification only, NO algorithm implemented and
NO external validation run.** Written *after* the original 19-dataset
publication benchmark was inspected, so nothing below should be
presented as independently preregistered prior to these outcomes.

**Evidence checkpoint:** [TSPDT_RECON_19x5_AUDITED_CHECKPOINT.md](TSPDT_RECON_19x5_AUDITED_CHECKPOINT.md),
which now records the full 285-fit RLR/RLMT five-fold audit,
95 TSPDT-inspired fold replay, 380-model read-only comparison,
nine DC discovery/confirmation gates, and the final 18-leaf
test-route/logloss decomposition.

## Observation to explain, without reusing test labels for design choices

- RLMT-DC-v0 accepted **9 of 95** outer-fold candidates;
  **6/9** accepted trees increased external-test logloss.
- Dataset-equal-weight overall DC-minus-RLR logloss:
  **+0.017217934**.
- **breast2/fold01**: confirmation gain +1.604980,
  estimated standard error 0.838936, gain−SE +0.766044,
  but outer-test DC-minus-RLR logloss **+1.698460**;
  both terminal leaves contribute positively.
- Four test observations in that fold had true-class
  probability below 0.01. Locally trained leaf models
  had 10 active relations per leaf from 26 and 35
  training samples respectively.
- **leukemia** had two beneficial accepted folds with
  logloss improvements −0.060396 and −0.150278.

These facts support a *stability/generalization risk*,
**not** a proved implementation fault, fixed numerical bug,
causal feature-selection failure, or validated intervention.

## Concrete source-code observation

Current `rlmt_dc_parent_and_child_on_validation()`
trains parent and child estimators on discovery partition D,
using a supervised pair pool generated on D, then estimates
parent-minus-child logloss on disjoint internal confirmation V.
`rlmt_dc_confirm()` accepts one prespecified nomination if
mean gain−standard error of per-observation gain is >0.
If accepted, `rlmt_dc_write_two_leaf_result()`
refits both leaves on **the full outer training T**,
using the common supervised T-level pair pool and
the node-local reranking behavior of `rlmt_fit_leaf()`.

The D-trained validation model is therefore **not the same
fitted estimator** evaluated on the outer test after T refit.
This is mathematically normal in model selection; nevertheless
the rule does not automatically protect against
model instability induced by moving from D to T.

## Scientific hypotheses for a separate successor method

**H1 — estimator-aligned internal validation.**
The winning relation and local-model induction
procedure should be evaluated by an inner training-only
resampling protocol reproducing the *complete final
refit pipeline*, including how supervised pair pools
are constructed, how leaf screening is repeated, and how
CV regularization parameters are selected.
Evaluation must use nested, strictly out-of-fold
predictions for each nominated candidate and its parent.
Keep the nomination independent of the confirmation
outcomes. Test whether training-only holdout paired
losses become more predictive of outer test harm when the
same full model-building algorithm is reproduced for
confirmation. This is a **hypothesis**, not a claim that
it will improve real-world logloss.

**H2 — conservative probabilistic leaf estimates.**
Test a principled hierarchical/shrunk probabilistic
alternative in which leaf models share information with,
or shrink to, the parent when node sample sizes are small.
Candidate concepts include penalized shared-plus-local
logistic coefficients, cross-fitted partial pooling,
or a train-only fitted gating weight. Any smoothing,
mixing weight and regularization must be learned **inside
each training fold**, with no access to outer test outcomes.
Do not apply a post-hoc fixed clipping level selected
from the breast2 error magnitudes.

Each hypothesis is separable. Prioritize implementing
**one** research-grade candidate at a time rather than
an unconstrained hyperparameter grid.

## Mandatory implementation gates before any new clinical data

1. Produce mathematical specification of estimator,
   loss, candidate nomination, acceptance, probability
   prediction and fixed decision threshold. Preserve
   distinction between development, internal
   confirmation, final full-training refit and outer test.
2. On synthetic generated data, test invariant
   predictions under saving/reloading, reproducibility
   of candidate ordering, identity to original RLR when
   no split is accepted, correct no-leakage of supervised
   feature ranking, and all failure/fallback branches.
3. Include synthetic high-confidence outlier examples
   and adversarial class-support patterns; test numerical
   finiteness and record positive-class ARFF factor level
   order explicitly.
4. Log per-candidate and per-fold model identity,
   train/validation sample IDs, local support, coefficient
   complexity, predicted probability range and
   acceptance evidence; never silently delete failures.
5. Benchmark against both **original published RLR**
   and the original frozen **RLMT-DC-v0** with identical
   outer folds, including negative outcomes.
6. Treat all 19 original publication datasets as exposed
   **development/retrospective** evidence.
   Choose at least one entirely new biologically
   defensible cohort and freeze its preprocessing,
   patient grouping, outer split protocol and
   primary metric **before** inspecting its outcomes.
7. Report **logloss as primary** because the hypothesized
   failure concerns probability reliability, and
   balanced accuracy, macro-F1, Brier, calibration
   diagnostics and model complexity as prespecified
   secondary descriptors.
8. Preserve all v0 and TSPDT-inspired v1 source/code SHAs,
   benchmark result directories, and previously
   published RLR paper artifacts. A method successor
   needs its own branch, version, independent result root
   and source audits.

**Decision:** Complete v0 diagnostic stage. No further
test-label-driven intervention on breast2. Next actionable
scientific work is to select H1 or H2 by mathematical
reasoning and available source evidence, then implement
a separate experimental prototype with synthetic-only
verification before touching new biomedical cohorts.

## Independent TSPDT track

TSPDT-inspired v1 (MAD100, max depth3, Wilson approximation)
is **not** original TSPDT-2011. A scientifically faithful
reconstruction track must first clarify the original
Relief-F feature selection, node-local TSP ranking
(including secondary tie handling), root-pair equivalence,
and C4.5-like pruning. Its version 2 work should remain
separate from RLMT-DC probability-stability research.
