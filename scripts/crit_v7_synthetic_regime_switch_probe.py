#!/usr/bin/env python3
"""CRIT-v7: deterministic SYNTHETIC representational falsification probe.

This is NOT original RLR, not yet a clinical model, not evidence of
unseen biomedical improvement. Baseline is a deliberately strong
ADDITIVE relational logistic model with ALL 45 pairs of 10 genes.
A literal two-leaf model adds a DIFFERENT branch intercept and slope
for ONE extra relation (genuine conditional interaction).

Two-phase discovery/confirmation is confined to TRAIN patients;
an independent outer TEST set provides the only reported performance.
Parameters and seed range are fixed before the holdout 8100..8199.
All gene-pair relation data are rank-order comparisons only.

Reproduce: OPENBLAS_NUM_THREADS=1 python3 \
    scripts/crit_v7_synthetic_regime_switch_probe.py
Requires: numpy, scipy, scikit-learn.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss


N_GENES = 10
N_TRAIN = 300
N_TEST = 1400
BASELINE_L2_C = 0.23
LEAF_RIDGE_PENALTY = 0.08
CANDIDATES_TO_CONFIRM = 8
DISCOVERY_FRACTION = 0.68
CONFIRMATION_SE_MULTIPLIER = 1.50
CONFIRMATION_GAIN_PENALTY = 0.015
FIRST_HOLDOUT_SEED = 8100
REPETITIONS = 100
PAIR_LIST = tuple(combinations(range(N_GENES), 2))
PAIR_INDEX = {pair: i for i, pair in enumerate(PAIR_LIST)}
DISJOINT_PAIR_MASK = np.asarray([
    [not (set(a) & set(b)) for b in PAIR_LIST]
    for a in PAIR_LIST
], dtype=bool)
np.fill_diagonal(DISJOINT_PAIR_MASK, False)


def rank_relations(x: np.ndarray) -> np.ndarray:
    """Pairwise within-patient expression order, no absolute units."""
    return np.column_stack([
        np.where(x[:, a] > x[:, b], 1.0, -1.0)
        for a, b in PAIR_LIST
    ])


def generate(n: int, seed: int, scenario: str) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    genes = rng.standard_normal((n, N_GENES))
    features = rank_relations(genes)
    gate = features[:, PAIR_INDEX[(0, 1)]]
    relation = features[:, PAIR_INDEX[(2, 3)]]
    if scenario == "switch":
        # The slope of the relation REVERSES between two gate branches.
        eta = 2.0 * gate * relation
    elif scenario == "global":
        # Strong stable additive signal: no split is truly required.
        eta = 1.4 * relation
    elif scenario == "noise":
        # Negative control: labels have no conditional signal at all.
        eta = np.zeros(n)
    else:
        raise ValueError(f"Unknown frozen synthetic scenario {scenario!r}")
    labels = rng.binomial(1, expit(eta)).astype(int)
    return features, labels


def fit_additive_baseline(x: np.ndarray, y: np.ndarray) -> LogisticRegression:
    clf = LogisticRegression(C=BASELINE_L2_C, max_iter=300, solver="lbfgs")
    clf.fit(x, y)
    return clf


def base_predict(clf: LogisticRegression, x: np.ndarray) -> np.ndarray:
    return clf.predict_proba(x)[:, 1]


def branch_design(gate: np.ndarray, relation: np.ndarray) -> np.ndarray:
    """Two leaf-specific intercepts and slopes: [1L,1R,r*1L,r*1R]."""
    left = (gate < 0).astype(float)
    right = (gate > 0).astype(float)
    return np.column_stack((left, right, left * relation, right * relation))


def logit(p: np.ndarray) -> np.ndarray:
    t = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(t / (1 - t))


def fit_leaf_offsets(
    y: np.ndarray, baseline_probability: np.ndarray, x: np.ndarray
) -> np.ndarray:
    fixed_offset = logit(baseline_probability)

    def loss_and_grad(t: np.ndarray) -> tuple[float, np.ndarray]:
        eta = fixed_offset + x @ t
        loss = (np.mean(np.logaddexp(0, eta) - y * eta)
                + 0.5 * LEAF_RIDGE_PENALTY * np.sum(t * t))
        grad = (x.T @ (expit(eta) - y)) / len(y) + LEAF_RIDGE_PENALTY * t
        return float(loss), grad

    result = minimize(
        loss_and_grad, np.zeros(x.shape[1]), method="BFGS", jac=True,
        options={"gtol": 1e-6, "maxiter": 100},
    )
    if not np.all(np.isfinite(result.x)):
        raise RuntimeError("Nonfinite two-leaf offset coefficients")
    return result.x


def individual_logloss(y: np.ndarray, p: np.ndarray) -> np.ndarray:
    q = np.clip(p, 1e-12, 1 - 1e-12)
    return -(y * np.log(q) + (1 - y) * np.log1p(-q))


def select_tree_on_outer_train(
    x: np.ndarray, y: np.ndarray, seed: int
) -> tuple[tuple[LogisticRegression, int | None, int | None,
                 np.ndarray | None], float, tuple[int, int]]:
    # Completely separate discovery D and confirmation V.
    order = np.random.default_rng(seed).permutation(len(y))
    n_d = int(DISCOVERY_FRACTION * len(y))
    d, v = order[:n_d], order[n_d:]

    root_d = fit_additive_baseline(x[d], y[d])
    pd = base_predict(root_d, x[d])
    pv = base_predict(root_d, x[v])

    # Detect COEFFICIENT INSTABILITY, not a scalar mean residual shift:
    # score(g,r) ∝ |sum_i (y_i - p_root_i) g_i r_i|.
    # Screen ONLY on D. The outer test and V labels are unseen here.
    score = np.abs(x[d].T @ (x[d] * (y[d] - pd)[:, None]))
    score[~DISJOINT_PAIR_MASK] = -np.inf
    shortlist = np.argsort(score.ravel())[-CANDIDATES_TO_CONFIRM:][::-1]

    best_gain = -np.inf
    best_lower_gain = -np.inf
    selected = None
    for flat_ix in shortlist:
        gate_j, model_j = np.unravel_index(flat_ix, score.shape)
        fd = branch_design(x[d, gate_j], x[d, model_j])
        theta = fit_leaf_offsets(y[d], pd, fd)
        fv = branch_design(x[v, gate_j], x[v, model_j])
        corrected = expit(logit(pv) + fv @ theta)

        # Paired per-patient logloss: nested train-only 1.5 SE gate.
        delta = individual_logloss(y[v], pv) - individual_logloss(
            y[v], corrected
        )
        mean_gain = float(np.mean(delta))
        conservative_gain = (mean_gain -
                             CONFIRMATION_SE_MULTIPLIER *
                             float(np.std(delta, ddof=1)) /
                             np.sqrt(len(v)))
        if conservative_gain > best_lower_gain:
            best_lower_gain = conservative_gain
            best_gain = mean_gain
            selected = (int(gate_j), int(model_j))

    assert selected is not None
    accepted = best_lower_gain > CONFIRMATION_GAIN_PENALTY
    # Full outer TRAIN refit is permitted, while accepted gate/pair
    # and the 2-leaf architecture are already frozen.
    full_root = fit_additive_baseline(x, y)
    if not accepted:
        return (full_root, None, None, None), best_gain, selected
    gate_j, model_j = selected
    corrected_design = branch_design(x[:, gate_j], x[:, model_j])
    leaf_theta = fit_leaf_offsets(y, base_predict(full_root, x),
                                  corrected_design)
    return (full_root, gate_j, model_j, leaf_theta), best_gain, selected


def one_independent_trial(seed: int, scenario: str) -> dict:
    xtr, ytr = generate(N_TRAIN, seed, scenario)
    xte, yte = generate(N_TEST, 100000 + seed, scenario)
    model, inner_gain, selected = select_tree_on_outer_train(
        xtr, ytr, seed + 1000
    )
    baseline_model, gate_j, model_j, leaf_theta = model
    root_prob = base_predict(baseline_model, xte)
    if gate_j is None:
        predicted = root_prob
    else:
        predicted = expit(logit(root_prob) +
                          branch_design(xte[:, gate_j],
                                        xte[:, model_j]) @ leaf_theta)
    root_ll = log_loss(yte, root_prob)
    tree_ll = log_loss(yte, predicted)
    return dict(seed=seed, scenario=scenario,
                activated=gate_j is not None, selected=selected,
                tree_minus_additive_ll=tree_ll - root_ll,
                root_ll=root_ll, tree_ll=tree_ll,
                best_train_validation_gain=inner_gain)


def main() -> None:
    target = {PAIR_INDEX[(0, 1)], PAIR_INDEX[(2, 3)]}
    for scenario in ("switch", "global", "noise"):
        records = [one_independent_trial(FIRST_HOLDOUT_SEED + i, scenario)
                   for i in range(REPETITIONS)]
        delta = np.array([r["tree_minus_additive_ll"] for r in records])
        activated = sum(r["activated"] for r in records)
        wins = int(np.sum(delta < -1e-10))
        losses = int(np.sum(delta > 1e-10))
        exact_recovered = sum(set(r["selected"]) == target for r in records)
        print(f"{scenario}: N={REPETITIONS} activated={activated} "
              f"W/L/T={wins}/{losses}/{REPETITIONS - wins - losses} "
              f"mean_dLL_tree_vs_additive={delta.mean():+.9f} "
              f"mean_additive_LL={np.mean([r['root_ll'] for r in records]):.9f} "
              f"exact_planted_pair_set_found={exact_recovered}/{REPETITIONS}")
    print("CRIT_V7_SYNTHETIC_ONLY=TRUE")
    print("BASELINE=ADDITIVE_ALL_45_PAIR_LOGISTIC_NOT_ORIGINAL_RLR")
    print("SCIENTIFIC_STAGE=SIMULATION_REPRESENTATIONAL_FALSIFICATION")
    print("NO_BIOMEDICAL_DATA_USED=TRUE")


if __name__ == "__main__":
    main()
