# RLMT-DC-v0 synthetic DEV gate ledger

This file tracks **observed** test results, not clinical/scientific performance.
GATE-3/4 original RLMT runs and published RLR must not be changed.

## DEV v1 — observed PASS

Execution SHA: `c1c29cd9a7b3f12da66c5e4a35e05893dcdfc5ff`.

- `RLMT_DC_PRIMITIVES=PASS`.
- `RLMT_DC_SYNTHETIC_L1_PARITY=PASS`.
- `RLMT_DC_SYNTHETIC_ROOT_PARITY=PASS`.
- `RLMT_DC_SYNTHETIC_PIPELINE=PASS`.
- Natural synthetic algorithm decision: `ACCEPTED_SPLIT=FALSE`.
- No end-to-end two-leaf replay at this point.

## DEV v2 — observed FAIL, **not** a successful gate

Execution SHA: `e6cd9aed77385f9efd8c9b598cbbdaa304dbfd76`.

- Parse, primitive regressions, strict discovery-OOF validation-label
  canary and RLR/L1 parity passed.
- Natural synthetic decision rejected the split and correct single-leaf
  recovery passed.
- The **deliberately forced** accepted two-leaf functional fixture
  completed training, saved artifacts, then failed the *combined*
  `"forced split structure/routing/confirmation"` assertion.
- Because multiple conditions were joined into one boolean assertion,
  **the exact failing predicate cannot be inferred from this log**.
- There is a concrete possible mechanism: the new JSON writer
  inherited `jsonlite::write_json`'s limited default numeric precision,
  whereas confirmation lower-bound gain was checked against the
  in-memory double with tolerance `1e-9`.
  This hypothesis requires re-execution; it is not an observed cause.

## DEV v3 — pending server execution

Changes in `R/new/42_rlmt_dc_run.R` and
`R/new/43_rlmt_dc_synthetic_smoke.R`:

1. Preserve confirmation gains and logloss at `digits=12` in
   **new RLMT-DC JSON artifacts only**, without changing selection,
   inner folds, class support or the frozen RLR implementation.
2. Replace the combined acceptance assertion by individual checks
   for tree flag, realized leaves, method identity, numerical precision,
   prediction sample IDs, route sample IDs and exact leaf routes.
3. Record stored vs original confirmation gain and its difference
   as `FORCED_CONFIRMATION_PRECISION` in synthetic output.
4. Run under separate `results/rlmt_dc_v0_devcheck_v3` to
   preserve the v1 and failed v2 evidence.
5. Require full serialized leaf-model prediction replay at
   tolerance `1e-12` as before.

**Passing DEV v3 may be declared only after the server reports**
`RLMT_DC_SYNTHETIC_ACCEPTED_BRANCH=PASS` and
`RLMT_DC_DEV_GATE=PASS`; a successful parser or primitive test alone
is insufficient.

The forced acceptance fixture is deliberately **not** a model-selection
or scientific performance result. A later, fully *naturally accepted*
synthetic experiment would give additional integration evidence but is
not a prerequisite for recording coverage of the forced persistence route.

No independent omics cohort has been tested with RLMT-DC; do not claim
predictive superiority.
