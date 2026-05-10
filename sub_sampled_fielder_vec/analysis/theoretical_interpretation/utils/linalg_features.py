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
