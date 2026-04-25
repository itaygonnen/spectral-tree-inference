"""Topology-specific code for the symmetric balanced binary tree.

Anything that depends on the balanced-binary structure lives here. Everything
else in this package (Fiedler computation, sub-sampling, caching, sweep,
plotting) is topology-agnostic so that other tree models (e.g. caterpillar,
sequence-based) can drop in their own builder module.
"""
from __future__ import annotations

import numpy as np


def build_balanced_binary_S(n: int, alpha: float, dtype=np.float32) -> np.ndarray:
    """Build the population similarity matrix S_ij = alpha ** d(i, j).

    For a symmetric balanced binary tree with n = 2**k leaves, the graph
    distance between leaves i and j is d(i, j) = 2 * (k - lca_level(i, j)).
    With leaves indexed 0..n-1 by position, the LCA depth-from-root is
    k - bit_length(i XOR j), so d(i, j) = 2 * bit_length(i XOR j) for i != j
    and 0 on the diagonal. This is computed in vectorised form via XOR.
    """
    if n < 2 or (n & (n - 1)) != 0:
        raise ValueError(f"n must be a power of 2 and >= 2, got {n}")
    if not 0.0 < alpha < 1.0:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")

    idx = np.arange(n, dtype=np.int64)
    xor = idx[:, None] ^ idx[None, :]
    bit_length = np.zeros_like(xor, dtype=np.int64)
    nz = xor > 0
    bit_length[nz] = np.floor(np.log2(xor[nz])).astype(np.int64) + 1
    d = 2 * bit_length
    return (np.asarray(alpha, dtype=dtype) ** d.astype(dtype)).astype(dtype, copy=False)


def balanced_binary_population_fiedler(n: int) -> np.ndarray:
    """Closed-form population Fiedler vector v^(1) = (1/sqrt(n)) * (+1..., -1...).

    Used to verify the optimal-coherence assumption (|v_i| = 1/sqrt(n)) before
    running a sub-sampling sweep.
    """
    if n < 2 or (n & (n - 1)) != 0:
        raise ValueError(f"n must be a power of 2 and >= 2, got {n}")
    half = n // 2
    v = np.empty(n, dtype=np.float64)
    v[:half] = 1.0
    v[half:] = -1.0
    v /= np.sqrt(n)
    return v
