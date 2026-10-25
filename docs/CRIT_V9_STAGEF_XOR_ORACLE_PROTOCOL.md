# CRIT-v9 Stage F — XOR discovery vs oracle expressivity

Protocol fixed on 2026-10-25; **REUSED_SYNTHETIC_DEV_MECHANISM_ONLY**.

## Why stage F is needed

Frozen v9 passed all synthetic structural gates, 17,280 patient-level heldout synthetic comparisons, original RLR benchmark parity, and 36 original biomedical method-fold vectors. Stage-v9 results:

| Mechanism (6 independent replicates) | CONDITIONAL_TREE accepted | Other controls |
|---|---:|---|
| No signal | 0/6 | GATE_ONLY 1/6 (test Δlogloss +0.001926, average across six +0.000321) |
| Global planted rank pair | 3/6 | ADDITIVE_PAIRS 4/6 |
| Planted gate-only label | 4/6 | ADDITIVE_PAIRS 4/6 |
| Planted XOR rank interaction | **0/6** | GLOBAL_INTERACTION 1/6 |

One exposed Chen fold 02 recovered v7's tree exactly, with a D-selected rank-six candidate, v8 rejected D rank 1. Stage v7–D demonstrated V-selection optimism under a constructed artificial-null stress test. The central remaining question is why genuine planted XOR fails even though a two-leaf rank tree can represent that mechanism.

The old 24 datasets have already been inspected, so **nothing in Stage F is an independent performance evaluation**. This is a frozen mechanistic debugging study.

## Scientific hypotheses, tested separately

H1 **Candidate coverage:** Is the planted gene01>gene02 gate and gene03>gene04 rank pair in the unsupervised/random/original-RLR candidate pool?

H2 **Discovery ranking:** Do any of the six D-screened candidate trees contain the planted gate and interaction pair? Does the D crossfit ranking choose the planted gate×pair? How high is the planted joint cross-product score relative to candidate screening scores?

H3 **Confirmation power:** If the planted gate and rank pair are given **as an oracle**, can a four-coefficient conditional tree pass the *exact existing v9* one-shot D-fit V-validation threshold? Can the equivalent explicit global interaction pass? Do not change any regularization or confirmation threshold based on past Chen or synthetic experiments.

H4 **Representational capacity and refit:** Does a full-T oracle conditional model improve independent synthetic E logloss **without a V filter**? Does the V-gated oracle preserve the benefit or reject it? Distinguish inability to fit the true structure from inability to find it or from low validation power.

## Experimental design

- Original frozen source CRIT-v9 SHA `7b32f804b7be92158a14fd5645eb65d56e7da763`; use a separate worktree and new result output.
- Regenerate **exactly six XOR and six label-independent null** studies from the already committed synthetic panel, each with 180 train / 120 heldout test patients, 32 measured independent Gaussian genes, predeclared seeds. Refit the genuine original RLR and CRIT-v9 once per study. Confirm exact model heldout probability vectors and E sample IDs / labels agree with the frozen v9 archive within 1e−10.
- Recreate the original v9 D/V split and D-only parent/OOf predictions; verify D-only candidate rank replay. Record pool size, planted pair presence, whether top-six candidates include correct gate/interaction, D selected rank and planted truth flags. Explicitly distinguish "gate is present" from "gate and interaction are jointly present".
- For two **oracle-known** designs: (a) two-leaf tree with opposite coefficients for gene03>gene04 conditional on gene01>gene02; (b) the equivalent global gate×rank-pair GLM. Train corrections using original D-only RLR OOF offsets with **unchanged v9 ridge and decision thresholds**; evaluate each exactly once on untouched V. Refit accepted oracle(s) on all T using strict original-RLR T-OOF. Also predict with full-T oracle *unconditionally* as a diagnostic of representational capacity.
- Store original RLR, discovered CRIT-v9 tree, conditional oracle (V-gated/unconditional) and global-interaction oracle (V-gated/unconditional): **12 studies × 6 arms × 120 heldout patients = 8,640 predictions**. Never choose an arm based on test E. All prior synthetic and biomedical result files unchanged.
- Independent Python audit recomputes test logloss/Brier and v9 replay parity, verifies D-candidate rankings, truth flags, one-shot V acceptance flags and exact RLR fallback.

## Interpretation

- If both true relations absent from pool: candidate discovery coverage is at fault; consider a mathematically justified broader low-dimensional pair pool but control `O(p^4)` gate×relation screening.
- If present, but absent from top-six proposals: discovery score/ranking and multiple-pair-search burden are plausible failure points.
- If top-six includes truth, but D choice differs: the D-only selection heuristic may be weak or distorted by reuse of D outcomes for nomination; this is **not** unbiased nested cross-validation.
- If oracle unconditionally improves heldout E but V rejects: one-shot confirmation may have insufficient power given 120 discovery/60 confirmation patients or a mismatched baseline. These known outcomes should inform **new independent synthetic studies**, not threshold tuning on the same E labels.
- If even oracle unconditionally fails on XOR: investigate correction model, original-RLR residual offsets, penalty, feature orientation, or signal-to-noise; avoid claiming a discovery-only bottleneck.

We are not doing a p-value or familywise-error calibration in Stage F. The six null replicates cannot quantify rare false acceptance. The correct future confirmation pathway is a larger, independently generated null/power benchmark and a genuinely unexposed biomedical cohort/batch after algorithm structure is fixed.

## Output and PASS gates

`results/crit_v9_stageF_xor_oracle_v1/`:

- `D_candidate_provenance.csv` — top-six true-feature coverage and selected rank for twelve studies
- `mechanism_summary.csv` — six-arm model statistics and V oracle decision per study
- `mechanism_test_patients.csv` — all 8,640 patient-level test probabilities and labels
- `oracle_reconstruction.log`, `independent_audit.log`, `status.txt` and protected source SHA files.

Completion requires:

`STAGEF_12_PAIRED_SYNTHETIC_STUDIES=PASS`,
`STAGEF_8640_PATIENT_PREDICTIONS=PASS`,
`STAGEF_8640_INDEPENDENT_PATIENT_REPLAY=PASS`,
`STAGEF_FROZEN_V9_PATIENT_PARITY=PASS`,
`STAGEF_BATCH_EXIT=0`.

All prior algorithms unchanged. All results labeled **reused synthetic development data**, not external validation.
