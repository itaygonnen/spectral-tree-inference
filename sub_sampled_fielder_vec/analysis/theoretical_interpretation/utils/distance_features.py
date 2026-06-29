"""Distance-side scalar features for the M1–M12 extension (see distance_extension.tex).

Companion to ``linalg_features.py`` (similarity / Laplacian side). Each helper
takes a plain dense array and returns either a scalar feature or a small dict
suitable for ``pd.DataFrame(rows)``.

The distance-side recovery operator is
    B = H D H,   H = I − (1/m) 1·1ᵀ,
and its leading eigenvector ``u^(1)`` replaces the Fiedler vector ``v^(2)`` of
``L = Deg(S) − S`` as the sign target.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import scipy.linalg


def tree_diameter(D: np.ndarray) -> float:
    """r(T) = max pairwise distance — the tree-diameter contamination factor in M7."""
    return float(np.asarray(D).max())


def beta0_d0_from_partition(
    S: np.ndarray,
    A_idx: np.ndarray,
    B_idx: np.ndarray,
    tr_q_norm: float = 1.0,
) -> Tuple[float, float]:
    """Cross-clan mean similarity ``β₀`` and matching distance scale ``d₀``.

    ``β₀ = mean_{i∈A, j∈B} S_ij`` and ``d₀ = −log β₀ / (4·‖tr Q‖)`` per M2 of
    distance_extension.tex. ``tr_q_norm`` defaults to 1 (absorbs the rate
    normalisation into a common constant when only the *ratio* d₀²/β₀² is used,
    as in A1 and Sim 1).
    """
    A_idx = np.asarray(A_idx)
    B_idx = np.asarray(B_idx)
    beta_0 = float(S[np.ix_(A_idx, B_idx)].mean())
    d_0 = -np.log(max(beta_0, 1e-12)) / (4.0 * float(tr_q_norm))
    return beta_0, d_0


def laplacian_features(S: np.ndarray) -> dict:
    """Spectral features of ``L = Deg(S) − S``.

    Returns: ``L``, top-3 ascending eigenvalues ``eigvals_L``, the Fiedler
    eigenvector ``v2`` (= eigvecs[:,1]), ``lambda2 = eigvals[1]``, and the
    spectral gap ``gap_L = eigvals[2] − eigvals[1]`` (the gap that controls
    Davis–Kahan stability of v2).
    """
    S = np.asarray(S)
    deg = S.sum(axis=1)
    L = np.diag(deg) - S
    eigvals, eigvecs = scipy.linalg.eigh(L, subset_by_index=[0, 2])
    return {
        "L": L,
        "eigvals_L": eigvals,
        "v2": np.asarray(eigvecs[:, 1], dtype=np.float64),
        "lambda2": float(eigvals[1]),
        "gap_L": float(eigvals[2] - eigvals[1]),
    }


def griffing_features(D: np.ndarray) -> dict:
    """Spectral features of ``B = H D H`` (Griffing's double-centered distance operator).

    ``B`` is negative-semi-definite for raw phylogenetic distance ``D`` (a
    classical MDS fact); the structural eigenvector is therefore the one with
    *largest absolute* eigenvalue — i.e. the most-negative one. The "top gap"
    of M6 is the gap between the dominant magnitude and the next, on the
    same end of the spectrum.

    Returns: ``B``, sorted-by-magnitude eigenvalues ``eigvals_B`` (descending
    |·|), the structural eigenvector ``u1`` (largest |eigvalue|), ``lambda1``
    and the *magnitude-gap* ``gap_B = |λ₁| − |λ₂|``.
    """
    D = np.asarray(D)
    row = D.mean(axis=1, keepdims=True)
    col = D.mean(axis=0, keepdims=True)
    grand = float(D.mean())
    B = D - row - col + grand
    B = 0.5 * (B + B.T)
    eigvals, eigvecs = scipy.linalg.eigh(B)
    order = np.argsort(np.abs(eigvals))[::-1]
    eigvals_by_mag = eigvals[order]
    eigvecs_by_mag = eigvecs[:, order]
    return {
        "B": B,
        "eigvals_B": eigvals_by_mag,
        "u1": np.asarray(eigvecs_by_mag[:, 0], dtype=np.float64),
        "lambda1": float(eigvals_by_mag[0]),
        "gap_B": float(abs(eigvals_by_mag[0]) - abs(eigvals_by_mag[1])),
    }


def distance_from_similarity(S: np.ndarray, tr_q_norm: float = 1.0) -> np.ndarray:
    """Tree-additive distance ``D†_ij = −log S_ij / (4·‖tr Q‖)`` (M1).

    Uses ``log(max(S, 1e-12))`` to keep the operator well-defined when entries
    of S are tiny (sub-sampled or numerically zero). The diagonal stays at 0
    after the negation because log(1) = 0 — but S has zero diagonal by
    construction in this codebase, so we explicitly zero the diagonal of D.
    """
    S = np.asarray(S)
    D = -np.log(np.clip(S, 1e-12, None)) / (4.0 * float(tr_q_norm))
    np.fill_diagonal(D, 0.0)
    return D
