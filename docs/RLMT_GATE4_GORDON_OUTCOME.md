# RLMT GATE-4 outcome: Gordon-2002 (2026-10-21)

**Technical audit: PASS. Evidence of RLMT-L2 improvement: NONE.**

This is the final **observed-result checkpoint** following the pre-outcome
freeze in [RLMT_GATE4_GORDON_FREEZE.md](RLMT_GATE4_GORDON_FREEZE.md).
The pre-outcome frozen configuration, publication RLR, and pilot artifacts
were not modified.

## Provenance

- GATE-4 input: `gordon-2002` from filtered historical CompCancer.
- Source integrity: 181 samples, 1626 numeric features, 150 AD / 31 MPM,
  no nonfinite features or exact duplicate rows, publication-manifest
  SHA256 `42bb0900f7a22851244c27243b626514917c0674d1a3c63dbf77bcab43361bab`.
- Model code at execution: `193fab667acfd4ab542ebc51c96cb3a3df31ddba`,
  written to the local `results/rlmt_gate4_gordon_frozen_5fold/CODE_SHA.txt`.
- Publication `main`: `4079a5ac400f1d2cfded523c7000f655470925ce`.
- Config: `configs/rlmt/gate4_gordon_frozen.json`.
- Five folds, `rep00_fold00` through `rep00_fold04`, seed 42,
  original RLR versus RLMT-L1 versus RLMT-L2.
- All outputs below were obtained on the compute server, stored in the
  ignored `results/` tree; they were validated by
  `R/new/38_rlmt_gate4_check.R` and
  `R/new/39_rlmt_gate4_screen_support.R`.

## Outcomes

| Outer fold | Test n | RLR logloss | RLMT-L2 logloss | Leaves |
| --- | ---: | ---: | ---: | ---: |
| 00 | 37 | 0.0015 | 0.0015 | 1 |
| 01 | 36 | 0.0309 | 0.0309 | 1 |
| 02 | 36 | 0.1868 | 0.1868 | 1 |
| 03 | 36 | 0.0027 | 0.0027 | 1 |
| 04 | 36 | 0.0032 | 0.0032 | 1 |

The mean fold logloss of RLR and RLMT-L2 was 0.045020 each;
mean paired difference `RLMT-L2 − RLR = 0`. Mean Brier difference
also equaled zero. RLMT-L1 matched the frozen RLR exactly in audited
probabilities/metrics. There were 0 accepted splits and hence no
two-leaf out-of-sample evaluation on Gordon.

**Do not report '0 of 100 split candidates was significant'**:
there were `0/0` lookahead evaluations; the candidate support
gate rejected the complete saved pool before screening.

## Why 0/0 lookahead: read-only reconstruction of the saved pool

Audit `R/new/39_rlmt_gate4_screen_support.R` checked **all**
1000 saved outer-training supervised global pairs in **each** of
five folds. There was no model refitting, no change to split gates,
and no use of outer-test outcomes for selecting a split.

| Check (across five folds) | Passed | Failed |
| --- | ---: | ---: |
| Both candidate children contain at least 16 training observations | 5000 | 0 |
| Both children contain at least 6 of **each** training class | 0 | 5000 |
| Both requirements simultaneously | 0 | 5000 |

The maximum observed minimum cell count over the 2 × 2
`child side × class` tables was **2** across the saved pairs
(per fold maximum: 2/2/2/2/2), versus the preregistered required **6**.
Maximum of the smaller MPM-side counts ranged from 2 to 3, depending
on fold. Every single rejected candidate passed the total-child-size
criterion and failed **only** the class-support criterion.

Although Gordon's class imbalance restricts feasible candidate splits,
the root outer-training folds contained 24–25 MPM samples, so **the
constraint of at least 6 MPM in both children is mathematically possible
in principle**. The result applies to the **saved supervised top-1000
candidate pool**, not to the whole `p(p−1)/2` space of possible pairs.
The diagnosis is therefore a combination of the fixed local-model
class-support requirement and the selected pool's splitting behavior,
not proof that Gordon has no usable relational partitions whatsoever.

Since the supervised top-ranked global pairs are selected for class
discrimination, strong class separation in the selected pool is a
plausible *mechanism* for low cross-class support in candidate
children. The audit establishes the observed class counts, not the
counterfactual performance or properties of unselected pairs.

## Scientific disposition

**GATE-4 technical: PASS. GATE-4 demonstration of RLMT improvement: NO.**

- There is no empirical benefit from L2 over the already strong root
  classifier on this dataset: RLMT-L2 remains the identical one-leaf RLR.
- This is **not** proof that an RLMT-like tree can never improve, nor
  that relaxed leaf thresholds are justified.
- GATE-3 Singh had one accepted split with poorer holdout logloss
  (0.5173 vs RLR 0.2900), associated with a very overconfident wrong
  probability in a small minority-class leaf.
- Together, these observations emphasize the tension between
  (i) searching for highly discriminative relation pairs for splits and
  (ii) requiring enough minority observations for stable probabilistic
  models in both child regions.
- Results from historical CompCancer do not constitute verification on
  independent external biomedical cohorts. Patient independence in
  processed ARFF was not established.
- `mean_gain-SE` among screen-selected candidates remains a model
  selection heuristic subject to selection optimism and reuse of
  supervised outer-training global pairs inside inner OOF.

No threshold relaxation, retrospective data exclusion, post-hoc
selection of folds, or modification of the frozen method is warranted
from GATE-4. The experiment should remain a fixed negative-result
record. Any revised tree-generation method must be separately named
and evaluated under a **new preregistered protocol** on further,
previously unexamined cases or cohorts.

## Candidate *research* questions for a separate future RLMT variant

1. Can candidate splits be derived from **residual heterogeneity and
   minimum class overlap**, rather than reusing only the strongest
   global discriminative pairs? Keep local leaf predictors and split
   proposal features conceptually distinct. Do not compute or tune
   an alternate candidate set on Gordon and call it independent.
2. Can data-dependent candidate discovery, orientation and residual
   screening be performed **within the inner training subset**, with
   candidate evaluation on strictly unseen inner validation subsets?
   The final split/leaf models should still use outer training only;
   this requires an explicitly different algorithm and potentially
   further nesting or sample splitting.
3. Can a decision rule quantify the **calibration/tail risk** of local
   conditional probabilities without tuning to the already examined
   Singh holdout? All model selection should use training-only data.

These are hypotheses and design alternatives, **not conclusions or
approved modifications of the original RLMT**.
