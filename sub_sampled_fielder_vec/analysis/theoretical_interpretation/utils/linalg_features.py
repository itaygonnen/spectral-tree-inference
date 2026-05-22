"""Linear-algebra-side scalar features for the rank-2 subspace U = [v⁽⁰⁾, v⁽¹⁾].

Topology-agnostic: every helper takes plain dense arrays and returns scalars.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np


def coherence_mu(U: np.ndarray) -> float:
    """μ(U) = (n / r) · max_i ‖e_iᵀ U‖² for an n×r orthonormal basis U."""
    n, r = U.shape
    row_norm_sq = np.einsum("ij,ij->i", U, U)
    return float((n / r) * row_norm_sq.max())


def compute_top_eigenpairs(L: np.ndarray, k: int = 3) -> Tuple[np.ndarray, np.ndarray]:
    """Return the k smallest eigenvalues and eigenvectors of symmetric L, ascending."""
    eigvals, eigvecs = np.linalg.eigh(L)
    return eigvals[:k], eigvecs[:, :k]


def spectral_gap(eigvals: np.ndarray) -> float:
    """Δλ = λ₃ − λ₂ from an ascending eigenvalue array."""
    return float(eigvals[2] - eigvals[1])


def compute_lemma04_row(
    S: np.ndarray,
    n1: int,
    n2: int,
    S_in: float,
    S_out: float,
    alpha: float,
) -> dict:
    """One ``(η, m)`` row of the Lemma 0.4 verification dataframe.

    Computes the HBM-corrected structural margin
    ``ρ = S_in * α^{D_max} − S_out`` together with the V.03 identity predictions
    (``μ_pred = (1+η)/2``, ``λ_2_pred = (n1+n2)·S_out``) and the Lemma-0.4 lower
    bound ``gap_lb = n_min·ρ − η·n_min·S_out``. Returns a flat dict suitable for
    ``pd.DataFrame(rows)`` — including ``gap_slack = gap_emp − gap_lb`` which
    should stay ≥ 0 under the corrected definition.
    """
    from .tree_features import hbm_d_max, hbm_s_in_min, imbalance_eta, n_min as _n_min

    n = n1 + n2
    eta = imbalance_eta(n1, n2)
    nmin = _n_min(n1, n2)
    d_max = hbm_d_max(max(n1, n2))
    s_in_min = hbm_s_in_min(S_in, alpha, max(n1, n2))
    rho = s_in_min - S_out

    L = np.diag(S.sum(axis=1)) - S
    eigvals, eigvecs = compute_top_eigenpairs(L, k=3)
    U = eigvecs[:, :2]

    mu_emp = coherence_mu(U)
    mu_pred = (1.0 + eta) / 2.0
    lambda2_emp = float(eigvals[1])
    lambda2_pred = n * S_out
    gap_emp = spectral_gap(eigvals)
    gap_lb = nmin * rho - eta * nmin * S_out

    return {
        "eta": eta, "m": n, "n1": n1, "n2": n2, "n_min": nmin,
        "D_max": d_max, "S_in_min": s_in_min, "rho": rho,
        "mu_emp": mu_emp,           "mu_pred": mu_pred,
        "lambda2_emp": lambda2_emp, "lambda2_pred": lambda2_pred,
        "gap_emp": gap_emp,         "gap_lb": gap_lb,
        "mu_residual":      mu_emp - mu_pred,
        "lambda2_residual": lambda2_emp - lambda2_pred,
        "gap_slack":        gap_emp - gap_lb,
    }
