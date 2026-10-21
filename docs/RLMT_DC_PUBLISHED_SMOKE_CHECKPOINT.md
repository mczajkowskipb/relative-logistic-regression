# RLMT-DC — publication-data first-fold pilot: observed checkpoint

**2026-10-21**. The following is based strictly on the supplied
compute-server terminal transcript; it is **not** a completed independent
audit until the corrected read-only scripts pass.

## Completed computation (observed)

- 19/19 historical publication datasets (CompCancer11 and legacy8);
  first paper-defined stratified fold `rep00_fold00`, seed 42.
- 57/57 method fits marked `OK` in the run manifest; **0 ERROR**.
  Models: unchanged RLR K10/m10/d2/pool1000; original RLMT-v1 L2;
  experimental RLMT-DC-v0.
- Detached runner exited with `RUNNER_EXIT=0`,
  `PUBLISHED_DATA_PILOT=PASS` in run log.
- **At least one accepted RLMT-DC split**, on `leukemia`
  (legacy8). Transcript logloss:
  `RLR = 0.10620000`, `RLMT-DC = 0.04578766`,
  `realized_leaves = 2`, `accepted_split = TRUE`.
  These are **single-fold exploratory logloss values**, not a
  claim of superiority, external validation or clinical value.
- Values for the remaining datasets, counts of all accepted
  splits and complete paired contrasts must be obtained from
  `run_manifest.csv` and independently audited before
  cross-dataset interpretation. The terminal transcript showed
  only its tail.

**Immutable execution code SHA** recorded by the original launcher:
`814fee71fe6306db1e62b2ca47b1038df82d77a8`.

Output root:
`results/rlmt_dc_published_data_smoke_v1`.

## Initial audit: observed FAIL, identified audit-only bug

The original `R/new/45_rlmt_dc_published_audit.R` stopped after
all 57 recorded fits were complete with:

```
'length = 57' in coercion to 'logical(1)'
```

The cause was two `which(...)` row filters using `&&` on
57-row logical vectors rather than elementwise `&`. That operation
raises the observed error in R. It occurred only in the
post-run audit; model fitting had already finished.

The fixed read-only auditor preserves the **original code SHA** above,
and no model parameters, predictions or original experiment outputs
are replaced. The separate
`R/new/46_rlmt_dc_published_accepted_replay.R` independently
reconstructs every actually accepted two-leaf model's test
probabilities from the original ARFF, saved tree routing,
coefficients, pair matrix and `cvfit.rds`; it also verifies
the confirmation criterion where relevant.

## Next release gate (pending actual server execution)

1. Run corrected `R/new/45_rlmt_dc_published_audit.R` on the
   existing complete 57-fit output. Require
   `PUBLISHED_DATA_AUDIT=PASS`.
2. Run the separate nonfitting accepted-leaf replay checker.
   Require `PUBLISHED_ACCEPTED_REPLAY=PASS`.
3. Analyze **all 19 datasets** including losses, accepted-split counts,
   class-support failures and no-split cases. No exclusion or
   threshold retuning based on these exploratory test outcomes.
4. Only after the two audit checks pass prepare a *separately
   registered* all-fold exploration. No repeat of the first-fold
   models is required to fix the audit.

The original RLR publication branch `main`
(`4079a5ac400f1d2cfded523c7000f655470925ce`)
and the frozen RLMT-v1 prototype branch are untouched.


## Completed independent audit — observed PASS (2026-10-21)

The operator ran the corrected audit scripts against the **original**
57-fit artifact tree (execution SHA `814fee71...`). Both completed:

- `PUBLISHED_DATA_AUDIT=PASS`; 57/57 model records, 0 failures.
- `PUBLISHED_ACCEPTED_REPLAY=PASS`; all 4 actually accepted
  two-leaf trees were reconstructed from original ARFF feature rows,
  `test_routes.csv`, saved `cvfit.rds`, and the stored relation
  coefficients. Replayed logloss matched within the locked tolerance.

Complete first-fold outcome relative to RLR:

| RLMT variant | Wins | Losses | Exact ties | Accepted two-leaf trees |
| --- | ---: | ---: | ---: | --- |
| RLMT-v1 | 0 | 1 | 18 | Singh |
| RLMT-DC-v0 | 2 | 1 | 16 | Gordon, breast2, leukemia |

Only changed RLMT-DC holdout logloss values:

| Dataset | Frozen RLR logloss | RLMT-DC logloss | Paired difference DC-RLR |
| --- | ---: | ---: | ---: |
| gordon-2002 | 0.00152436 | 0.00057627 | -0.000948088 |
| breast2 | 1.01415402 | 1.14046901 | +0.126314990 |
| leukemia | 0.10618332 | 0.04578766 | -0.060395658 |

RLMT-v1 accepted Singh and was worse: 0.51731647 versus
RLR 0.29002105 (difference +0.227295421). All remaining
per-dataset/variant first-fold differences were zero.

**No claim of superiority:** this is one shared first-fold split per
dataset and all datasets belong to the already exposed publication /
RLMT development corpus. The result is a descriptive development
checkpoint, not independent scientific validation. Follow-up work
must preserve all 19 datasets and both gains and harms without changing
candidate thresholds based on these holdout results.
