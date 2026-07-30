"""Flat constant block model (CBM): a non-balanced two-clan similarity matrix.

Topology-specific module for the features ↔ linear-algebra notebook. The clans
have arbitrary sizes (n1, n2) so the imbalance ratio η = max(n1, n2) / min
sweeps freely over [1, n−1]. Within-clan entries equal S_in, cross-clan entries
equal S_out, the diagonal is zero (matrix-completion convention).
"""
from __future__ import annotations

import numpy as np


def build_flat_cbm_S(
    n1: int, n2: int, S_in: float, S_out: float, dtype=np.float32
) -> np.ndarray:
    """Return the (n1+n2)×(n1+n2) flat-CBM similarity matrix with zero diagonal."""
    n = n1 + n2
    S = np.full((n, n), S_out, dtype=dtype)
    S[:n1, :n1] = S_in
    S[n1:, n1:] = S_in
    np.fill_diagonal(S, 0.0)
    return S


def _within_clan_dist_matrix(size: int) -> np.ndarray:
    """Tree-distance matrix for a clan of ``size`` leaves under greedy near-equal binary split.

    The clan is recursively split into halves of sizes ``size//2`` and
    ``size - size//2`` until each leaf is alone. ``dist[i, j]`` is the depth (in
    merges from the leaves) of the smallest subtree containing both ``i`` and
    ``j``. Sibling pairs in any subtree get ``dist = 1``; cousins ``2``; pairs
    on opposite halves of the clan root ``ceil(log2(size))``.
    """
    if size <= 1:
        return np.zeros((size, size), dtype=np.int32)
    dist = np.zeros((size, size), dtype=np.int32)

    def recurse(start: int, end: int) -> int:
        sz = end - start
        if sz <= 1:
            return 0
        mid = start + sz // 2
        left_depth = recurse(start, mid)
        right_depth = recurse(mid, end)
        d = max(left_depth, right_depth) + 1
        dist[start:mid, mid:end] = d
        dist[mid:end, start:mid] = d
        return d

    recurse(0, size)
    return dist


def build_decay_cbm_S(
    n1: int,
    n2: int,
    S_in: float,
    S_out: float,
    alpha: float = 1.0,
    dtype=np.float32,
) -> np.ndarray:
    """Two-clan hierarchical CBM with geometric decay along **tree distance**.

    Within each clan, leaves are organised by recursive greedy near-equal
    binary split (sizes ``n_k//2`` and ``n_k - n_k//2`` at every level). For
    ``i != j`` in the same clan,

        ``S_ij = S_in * alpha ** dist(i, j)``,

    where ``dist(i, j)`` is the merge-depth of their smallest common subtree
    (siblings = 1, cousins = 2, ..., opposite halves of the clan root =
    ``ceil(log2(n_k))``). Cross-clan entries are constant ``S_out`` and the
    diagonal is zero (matrix-completion convention).

    ``alpha = 1.0`` collapses every within-clan entry to ``S_in`` and recovers
    :func:`build_flat_cbm_S` exactly. ``alpha < 1`` produces the nested-block
    pattern of a hierarchical CBM, with maximum off-diagonal within-clan entry
    ``S_in * alpha`` (siblings) — matching the ``S_in_max_eff = S_in * alpha``
    convention used by the V.03 identities.
    """
    n = n1 + n2
    S = np.full((n, n), S_out, dtype=dtype)
    for offset, size in [(0, n1), (n1, n2)]:
        if size <= 1:
            continue
        dist = _within_clan_dist_matrix(size)
        block = (S_in * np.power(float(alpha), dist)).astype(dtype)
        S[offset:offset + size, offset:offset + size] = block
    np.fill_diagonal(S, 0.0)
    return S


def flat_cbm_population_fiedler(n1: int, n2: int) -> np.ndarray:
    """Closed-form unit-norm Fiedler vector for the flat-CBM Laplacian.

    Constant ``a`` inside C₁ and ``b`` inside C₂, with n1·a + n2·b = 0
    (orthogonality to v⁽⁰⁾) and n1·a² + n2·b² = 1 (unit norm).
    """
    n = n1 + n2
    a = np.sqrt(n2 / (n * n1))
    b = -np.sqrt(n1 / (n * n2))
    v = np.empty(n, dtype=np.float64)
    v[:n1] = a
    v[n1:] = b
    return v
