"""Tree-side scalar features used in the V.03 spectral-gap relations.

Each function is a single scalar identity, kept separate so callers can mix
and match (e.g. in a DataFrame .apply over a sweep grid).
"""
from __future__ import annotations

import numpy as np


def imbalance_eta(n1: int, n2: int) -> float:
    """η = max(n1, n2) / min(n1, n2) ∈ [1, ∞)."""
    return float(max(n1, n2)) / float(min(n1, n2))


def n_min(n1: int, n2: int) -> int:
    """n_min = min(n1, n2) — the smaller clan size."""
    return int(min(n1, n2))


def structural_margin_rho(
    S_in_max: float, S_out_max: float, S_out_min: float
) -> float:
    """ρ = (S_in^max − S_out^max) − (S_out^max − S_out^min).

    Under the molecular clock (single S_out value, so S_out_min = S_out_max),
    this collapses to ρ = S_in − S_out.
    """
    return (S_in_max - S_out_max) - (S_out_max - S_out_min)


def estimate_features_from_M(M: np.ndarray, v_pop: np.ndarray) -> dict:
    """Infer (η, ρ, S_in/out_max/min, n_min margin) from a real similarity matrix.

    Splits taxa into two clans by the sign of ``v_pop`` (the reference Fiedler /
    population partition), then reads off the empirical extrema of the
    intra-clan vs cross-clan entries of ``M`` to instantiate the non-balanced
    CBM bound parameters.
    """
    pos = np.where(v_pop > 0)[0]
    neg = np.where(v_pop <= 0)[0]
    n1, n2 = int(len(pos)), int(len(neg))
    eta = imbalance_eta(n1, n2)

    M_in_pos = M[np.ix_(pos, pos)]
    M_in_neg = M[np.ix_(neg, neg)]
    iu_pos = np.triu_indices(len(pos), k=1)
    iu_neg = np.triu_indices(len(neg), k=1)
    in_vals = np.concatenate([M_in_pos[iu_pos], M_in_neg[iu_neg]])
    out_vals = M[np.ix_(pos, neg)].ravel()

    s_in_max = float(in_vals.max())
    s_out_max = float(out_vals.max())
    s_out_min = float(out_vals.min())
    rho = structural_margin_rho(s_in_max, s_out_max, s_out_min)
    margin = rho - eta * s_out_max
    return dict(
        n1=n1, n2=n2, eta=eta,
        S_in_max=s_in_max, S_out_max=s_out_max, S_out_min=s_out_min,
        rho=rho, margin=margin,
    )
