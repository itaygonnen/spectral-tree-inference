"""Griffing's distance-matrix partitioning helpers.

The SNJ paper (Jaffe & Kluger, Section 4.1 + Appendix D) references the
classical distance-based partitioning of Griffing (Chapter 4 of Griffing,
2012). Given the matrix of pairwise distances ``D ∈ R^{m×m}`` between
terminal nodes, the method:

  1. Double-mean-centers: ``B = (I − 11ᵀ/m) D (I − 11ᵀ/m)``.
  2. Takes the eigenvector ``v`` with largest |eigenvalue| of ``B``.
  3. Partitions terminal nodes by ``sign(v)``.

Why this is interesting under sub-sampling: if the IPW estimator gives
``D̂ ≈ D/p`` (deterministic systematic bias under mean-imputation), then
``B̂ ≈ B/p`` and the leading eigenvector is **invariant under positive
scalar multiplication**. The *sign pattern* survives the bias even when
``||D − D̂||_2`` is enormous — which is exactly the regime where the
"more samples → better" effect should become visible (the stochastic
component of D̂ shrinks with n).
"""
from __future__ import annotations

import numpy as np
import scipy.linalg


def griffing_centered(D: np.ndarray) -> np.ndarray:
    """Return ``B = (I − 11ᵀ/m) D (I − 11ᵀ/m)`` for a (possibly non-symmetric) D.

    Implemented in ``O(m²)`` via row/col/grand means rather than the explicit
    matrix product (which would be ``O(m³)``). For symmetric D the result is
    symmetric.
    """
    row = D.mean(axis=1, keepdims=True)
    col = D.mean(axis=0, keepdims=True)
    grand = float(D.mean())
    return D - row - col + grand


def griffing_leading_eigvec(D: np.ndarray) -> np.ndarray:
    """Return the eigenvector of ``B = J D J`` with the largest |eigenvalue|.

    B is symmetric for symmetric D. Uses ``scipy.linalg.eigh`` on the
    symmetrized B; the eigenvalue of largest magnitude corresponds to the
    dominant phylogenetic axis (the first coordinate of classical MDS),
    which Griffing partitions by sign.
    """
    B = griffing_centered(D)
    B = 0.5 * (B + B.T)
    eigvals, eigvecs = scipy.linalg.eigh(B)
    idx = int(np.argmax(np.abs(eigvals)))
    return np.asarray(eigvecs[:, idx], dtype=np.float64)
