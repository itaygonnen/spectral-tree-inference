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
from scipy.sparse.linalg import ArpackError, ArpackNoConvergence, eigsh

from .logging import log_info

# Solver for the one eigenpair this module actually needs.
#   "eigh_full" -- LAPACK on the whole spectrum, then argmax|lambda|. The
#                  historical path; kept as the DEFAULT so every existing caller
#                  and every cached sweep key stays byte-identical.
#   "lm_k1"     -- ARPACK for the single largest-magnitude eigenpair. Same vector
#                  (verified |cos| = 1 to 8 decimals and identical sign
#                  partitions on the Kingman pool at n = 500..6000) at a fraction
#                  of the cost, because eigh_full is O(m^3) to use one column:
#                  14.94 s vs 0.12 s at m = 6000. That gap is what makes a
#                  25-point p-grid x 10 reps sweep over the pool feasible at all.
SOLVERS = ("eigh_full", "lm_k1")
DEFAULT_SOLVER = "eigh_full"


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


def griffing_leading_eigvec(
    D: np.ndarray, solver: str = DEFAULT_SOLVER,
) -> np.ndarray:
    """Return the eigenvector of ``B = J D J`` with the largest |eigenvalue|.

    B is symmetric for symmetric D, and negative semi-definite for a raw
    (non-squared) distance matrix, so the eigenvalue of largest *magnitude* is
    the most negative one -- that is the dominant phylogenetic axis (the first
    coordinate of classical MDS), which Griffing partitions by sign. Picking the
    largest *algebraic* eigenvalue instead lands on the degenerate ~0 of the
    all-ones null space and gives a random partition.

    ``solver`` selects how that single eigenpair is obtained; see :data:`SOLVERS`.
    Both branches return the same vector up to sign, and callers here always
    sign-align against a reference, so they are interchangeable. ``"lm_k1"``
    falls back to the dense path when ARPACK cannot help or cannot converge.
    """
    if solver not in SOLVERS:
        raise ValueError(f"solver must be one of {SOLVERS}, got {solver!r}")
    B = griffing_centered(D)
    B = 0.5 * (B + B.T)

    # ARPACK needs k < m-1, and below a few dozen rows the dense solve is
    # already faster than setting up an iterative one.
    if solver == "lm_k1" and B.shape[0] > 3:
        try:
            _vals, vecs = eigsh(B, k=1, which="LM")
            return np.asarray(vecs[:, 0], dtype=np.float64)
        except (ArpackNoConvergence, ArpackError) as exc:
            log_info("griffing",
                     f"eigsh(k=1, LM) failed on m={B.shape[0]} ({exc}); "
                     f"falling back to full eigh", force=True)

    eigvals, eigvecs = scipy.linalg.eigh(B)
    idx = int(np.argmax(np.abs(eigvals)))
    return np.asarray(eigvecs[:, idx], dtype=np.float64)
