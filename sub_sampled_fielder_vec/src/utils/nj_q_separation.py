"""Fig. 3 analog for classical NJ: distribution of the Q-criterion

    Q(i, j) = (r - 2) * D̂(i, j) - sum_k D̂(k, i) - sum_k D̂(k, j)

at the *initial* NJ step (all singletons, r = m), for adjacent vs.
non-adjacent leaf pairs. NJ merges argmin_{i≠j} Q(i,j), so a clean
separation between the Q distributions for cherries vs. non-cherries at
step 0 is the analog of σ₂ separation in the SNJ paper's Fig. 3.

Adjacency is defined on the *true* tree: a pair is adjacent iff the two
leaves form a cherry (share an immediate parent with no other children).
We reuse ``cherry_mask_from_tree`` from ``snj_sigma2_separation`` so the
adjacency set matches the SNJ experiment exactly.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np

from .snj_sigma2_separation import cherry_mask_from_tree


def _pairwise_q_singletons(D_hat: np.ndarray) -> np.ndarray:
    """Return an m×m symmetric matrix Q where

        Q[i, j] = (m - 2) * D̂[i, j] - sum_k D̂[k, i] - sum_k D̂[k, j].

    Computed vectorised over all (i, j); diagonal is zeroed since i ≠ j.
    """
    D = np.asarray(D_hat, dtype=np.float64)
    m = D.shape[0]
    row_sum = D.sum(axis=1)
    Q = (m - 2) * D - row_sum[:, None] - row_sum[None, :]
    np.fill_diagonal(Q, 0.0)
    Q = 0.5 * (Q + Q.T)  # numerical symmetrise
    return Q


def compute_q_separation(D_hat: np.ndarray, tree, taxa_metadata
                         ) -> Tuple[np.ndarray, np.ndarray]:
    """Return (adj, nonadj) — arrays of Q values for adjacent vs.
    non-adjacent leaf pairs at step 0 of NJ.

    Each array has length equal to the count of (i, j) pairs in the
    respective category, taking the upper triangle only (i < j).
    """
    m = D_hat.shape[0]
    Q = _pairwise_q_singletons(D_hat)
    cherry = cherry_mask_from_tree(tree, taxa_metadata)

    iu = np.triu_indices(m, k=1)
    q_upper = Q[iu]
    cherry_upper = cherry[iu]

    adj = q_upper[cherry_upper]
    nonadj = q_upper[~cherry_upper]
    return adj, nonadj
