"""Topology-agnostic spectral operations on a similarity matrix.

`compute_fiedler_of_S` builds the unnormalised Laplacian L = D - S and returns
the eigenvector of its second-smallest eigenvalue (the Fiedler vector). For an
S with constant row sums (e.g. balanced binary), the eigenvectors of S and L
coincide with reversed ordering, so this is equivalent to the second-largest
eigenvector of S used in the thesis.

`subsample_S` applies a symmetric Bernoulli(p) mask and rescales kept entries
by 1/p, leaving the diagonal untouched (matrix-completion convention).
"""
from __future__ import annotations

import numpy as np

from src.core.fiedler_computer import FiedlerVectorComputer
from src.core.utils import compute_laplacian


def compute_fiedler_of_S(S: np.ndarray, sampling_prob: float | None = None) -> np.ndarray:
    """Return the Fiedler vector of L = D - S with the project's sign convention."""
    laplacian = compute_laplacian(S)
    return FiedlerVectorComputer().compute(laplacian, sampling_prob=sampling_prob)


def subsample_S(S: np.ndarray, p: float, seed: int) -> np.ndarray:
    """Symmetric Bernoulli(p) sub-sample of S with 1/p rescaling on kept entries.

    The diagonal of S is preserved as-is (Sii = 1 in the population model and is
    always observed). Off-diagonal entries are kept independently with
    probability p across the upper triangle and mirrored to keep S_hat
    symmetric, then rescaled by 1/p.
    """
    if not 0.0 < p <= 1.0:
        raise ValueError(f"p must be in (0, 1], got {p}")
    n = S.shape[0]
    rng = np.random.default_rng(seed)

    upper_mask = rng.random((n, n), dtype=np.float32) < p
    upper_mask = np.triu(upper_mask, k=1)
    mask = upper_mask | upper_mask.T

    S_hat = np.zeros_like(S)
    S_hat[mask] = S[mask] / np.asarray(p, dtype=S.dtype)
    np.fill_diagonal(S_hat, np.diag(S))
    return S_hat
