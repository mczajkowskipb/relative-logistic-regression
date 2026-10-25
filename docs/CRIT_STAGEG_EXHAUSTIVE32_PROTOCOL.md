# Stage G — complete 32-gene rank-pair candidate-pool intervention

Protocol frozen 2026-10-25. **SCIENTIFIC_STAGE=SYNTHETIC_DEVELOPMENT_ONLY**. No independent biomedical performance inference.

## Established mechanism (Stage F; completed independently audited)

On six XOR synthetic studies, the planted gene01>gene02 gate was absent from the discovery proposal pool **6/6**, while planted gene03>gene04 was present **6/6**. Consequently, CRIT-v9 could not nominate the exact pair and rejected all 6 trees. An oracle tree supplied both true pairs passed the original V gate 6/6 and achieved **mean heldout Δlogloss −0.2167817** relative to unchanged original RLR; a global interaction oracle was V accepted 6/6 with **mean Δlogloss −0.2877593**. Null oracle accepted 1/6 and worsened heldout logloss, underscoring the need for later statistical calibration.

The original `crit7_pool` assembles 12 rank-variable genes, original RLR pairs and up to 256 pseudorandom pairs. Its pseudorandom function is called with fixed seed `cfg$seed+617`, so the sampled pair subset is repeated across replicates sharing feature labels. The observed 0/6 omission is therefore **not six independent misses of randomized discovery**. None of the original RLR, v9 or Stage F sources may be changed.

## Isolated intervention

For a synthetic feature matrix with exactly **32 uniquely named genes**, enumerate **all C(32,2)=496 undirected canonical pairs** as potential gate and leaf relations. No outcome or heldout label enters this enumeration. Replace **only** the v9 D proposal pool via a separately named `critG_complete32_pool`. Preserve all other v9 code paths and fixed parameters:

- exact published ORIGINAL_RLR fitted on T, D and OOF;
- v7 D-only residual interaction screening;
- shortlist at most 6 candidates;
- v9 discovery-only internal CV criterion for choosing **one** candidate;
- same independent confirmation V, one test per five predeclared control arms, exactly unchanged conservative threshold and regularization;
- same accepted/refit model, exact original RLR fallback for rejected corrections;
- the same independent E test-set procedure.

The copied, namespaced `critG_fit` and `critG_predict` are used only in this isolated worktree. They are not replacements for the original files. **An exhaustive 496-pair synthetic pool is not a deployable multi-omics candidate generator:** genome-wide C(p,2) and gate×pair screening have prohibitive scaling. It is an *ablation* that isolates candidate-pool failure.

## Fixed data

1. **24 reused developmental simulations:** all four frozen v9 scenario families (null, global_pair, gate_only, xor_interaction) × 6 original seeds, 180 train and 120 E-test patients, 32 Gaussian genes, seeded Bernoulli labels. Reproduce archived original-RLR E sample IDs, labels and probabilities exactly (1e-10). These are exposed data.
2. **24 additional newly seeded independent synthetic DEVELOPMENT simulations:** same four predeclared mechanisms × 6 new seeds per mechanism, original n=180/120 and class support. The new random samples are independent from the old ones, but the **algorithm was explicitly designed after examining the old XOR data**; consequently these are further development experiments, not independent evaluation of an unselected research hypothesis.
3. In every study, report planted gene01>gene02 and gene03>gene04 coverage in the *unchanged* original `crit7_pool` and in the full 496-pair intervention; true-gate/top-six proposal count, true joint gate×relation proposal count, exact D selected gate×relation, E logloss/Brier and primary conditional-tree accept flag. All arms saved, no heldout E labels used for selection.

Total **48 × 6 = 288 model-study comparisons** and **48 × 6 × 120 = 34,560 heldout patient prediction rows**, with exact original RLR parity on reused 24 studies. Python independently audits all patient metrics, candidate provenance, thresholds and exact fallback. No biological cohort training/score repetition in Stage G; no originals overwritten.

## Scientific outcome and decision gates

- If full pool retrieves planted gate yet the top six does not, then `crit7_nominate` interaction screening / shortlist dominates the failure. The next method should improve supervised discovery on D under computational constraints.
- If the planted gate reaches top six, but the internal v9 D CV does not select it, then the D scoring heuristic is the bottleneck. Feature nomination itself uses all D labels, so reported D CV is **not unbiased nested CV**.
- If true joint selected on D and V still rejects, inspect confirmation variance and effective V sample size without adjusting thresholds on reused E outcomes.
- If V accepts and E improves, pool coverage can rescue conditional modeling in this 32-gene simulation, but still no proof for high-dimensional biomedical cohorts.
- Compare no-signal acceptances and test harms separately. Six new null replicates are too few to estimate real type-I error.
- If **GLOBAL_INTERACTION** outperforms **CONDITIONAL_TREE**, note original distinct ridge parameterizations despite overlapping representation and require dedicated head-to-head experiments; do not interpret it as automatic superiority of a tree.

## Source and outputs

- Frozen original v9 ref SHA `7b32f804b7be92158a14fd5645eb65d56e7da763`.
- Stage F evidence SHA `94d2bc049bf987d08d87637847427c00cf1acee9`.
- New isolated branch `analysis/crit-v9-stageG-exhaustive32-ablation-v1`.
- Output `results/crit_v9_stageG_exhaustive32_v1/`: `study_summary.csv`, `candidate_provenance.csv`, `test_patients.csv`, `48_studies.log`, `independent_audit.log`, `status.txt`, source SHA files.

Required runtime PASS:

`STAGEG_48_STUDIES_288_METHODS=PASS`;
`STAGEG_34560_TEST_PATIENT_VECTORS=PASS`;
`STAGEG_EXHAUSTIVE_496_RELATIONS=PASS`;
`STAGEG_34560_PATIENT_METRICS_INDEPENDENTLY_AUDITED=PASS`;
`STAGEG_24_REUSED_STUDIES_EXACT_ORIGINAL_RLR_PARITY=PASS`;
`STAGEG_BATCH_EXIT=0`.

**No tuning of candidate count, penalties or decision thresholds on observed results.**
