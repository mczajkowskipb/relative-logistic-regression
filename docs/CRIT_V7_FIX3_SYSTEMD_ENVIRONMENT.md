# CRIT-v7 fix3 — systemd-run invalid environment block

## Observed Ubuntu failure (2026-10-25)

Operator launched fix2 from a new isolated worktree at
`a4e2a3da24b72e1a61e3c4c404073bd4f2286791`.
Source and worktree checks passed, user linger was enabled.
The launcher then failed **before service creation**:

```text
Failed to start transient service unit: Invalid environment block.
```

## Root cause

Fix2 passed `--setenv=PATH` and `--setenv=HOME`, and optional
R environment variables without values, to `systemd-run`.
Each `--setenv` argument must provide a valid `NAME=VALUE` assignment.
This is an infrastructure failure; it says nothing about the synthetic
or biomedical performance of CRIT.

## Isolated technical fix

- Launcher now passes `--setenv=PATH=<caller PATH>` and
  `--setenv=HOME=<caller HOME>`, with documented fallbacks.
- Optional `R_LIBS_USER`, `R_LIBS_SITE`, `LD_LIBRARY_PATH`
  values are passed as complete `NAME=VALUE` arguments.
- Refuses any embedded NUL, carriage return, or newline in values.
- Use branch `feature/crit-v7-systemd-env-fix3` and unique service `crit7-fix3.service`.
- Update guarded batch's identity check to that new pinned branch.
- Use a new worktree named `relative-logistic-regression-crit-v7-fix3`.
- Preserve failed fix1 and fix2 worktrees and logs.
- Do not change any original RLR, CRIT-v7 scientific algorithm,
  synthetic generator, biomedical pilot, or independent auditor source.

## Runtime status

The source-level fix is committed to GitHub. No successful
`systemd-run` launch, synthetic R test, biomedical fit, or independent
audited result is claimed until executed on the operator's Ubuntu host.

Scientific stage remains `EXPOSED_RETROSPECTIVE_DEV_ONLY`.
