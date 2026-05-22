"""Fig. 3 analog: distribution of Λ(i,j) = σ₂(R̂^{A_i ∪ A_j}) for adjacent vs.
non-adjacent leaf pairs at the *initial* SNJ step (all singletons).

At step 0, A_i = {i} and A_j = {j}, so the cross-block matrix is
``R̂[[i, j], :]`` with columns i and j removed — a 2×(m-2) slice. σ₂ of a
2×k matrix has a closed form via the 2×2 Gram matrix M M^T, which lets us
vectorise over all pairs (i, j) using only one m×m matrix-multiply.

Adjacency is defined on the *true* tree: a pair is adjacent iff the two
leaves share an immediate parent with no other leaves attached (i.e., they
form a cherry).
"""
from __future__ import annotations

from typing import Tuple

import numpy as np


def _pairwise_sigma2_singletons(R_hat: np.ndarray) -> np.ndarray:
    """Return an m×m symmetric matrix L where L[i, j] = σ₂(R̂[[i,j], ~{i,j}]).

    For a 2×k matrix M, the singular values are sqrt of eigenvalues of
    M M^T (a 2×2 matrix), so σ₂² = (tr − sqrt(tr² − 4·det)) / 2. We assemble
    the 2×2 Gram entries for every pair (i, j) from the full m×m product
    ``G = R̂ R̂^T``, after subtracting off the contributions of the two
    excluded columns i and j.
    """
    R = np.asarray(R_hat, dtype=np.float64)
    m = R.shape[0]

    # G[i, j] = <R[i,:], R[j,:]> (all m columns, includes columns i, j).
    G = R @ R.T

    # Subtract off contributions of the excluded columns i and j to get
    # <R[i, ~{i,j}], R[i, ~{i,j}]>, etc.
    diag_R_sq_i = R ** 2  # used for diag terms
    # a[i, j] = G[i, i] - R[i, i]^2 - R[i, j]^2
    a = G.diagonal()[:, None] - (R.diagonal() ** 2)[:, None] - R ** 2
    # b[i, j] = G[j, j] - R[j, i]^2 - R[j, j]^2
    b = G.diagonal()[None, :] - (R.T) ** 2 - (R.diagonal() ** 2)[None, :]
    # c[i, j] = G[i, j] - R[i, i]*R[j, i] - R[i, j]*R[j, j]
    c = G - R.diagonal()[:, None] * R.T - R * R.diagonal()[None, :]

    trace = a + b
    det = a * b - c ** 2
    disc = np.maximum(trace ** 2 - 4.0 * det, 0.0)  # numerical floor at 0
    sigma2_sq = 0.5 * (trace - np.sqrt(disc))
    sigma2_sq = np.maximum(sigma2_sq, 0.0)
    sigma2 = np.sqrt(sigma2_sq)

    # Diagonal isn't a valid pair (i == j); blank it out.
    np.fill_diagonal(sigma2, 0.0)
    # Symmetrise (formula is symmetric up to FP noise).
    sigma2 = 0.5 * (sigma2 + sigma2.T)

    # Suppress used-but-unread warning from linters; this is kept for clarity.
    _ = diag_R_sq_i
    return sigma2


def cherry_mask_from_tree(tree, taxa_metadata) -> np.ndarray:
    """Return an m×m boolean mask where True iff (i, j) is a cherry on `tree`.

    A cherry is a pair of leaves whose immediate parent has no other leaf
    children attached (and ideally no internal-node siblings, which would
    mean the parent has more than two children). We accept the standard
    "exactly two leaf children, no other children" definition.
    """
    m = len(taxa_metadata)
    mask = np.zeros((m, m), dtype=bool)

    parent_to_leaves: dict = {}
    parent_to_total_children: dict = {}
    for leaf in tree.leaf_node_iter():
        parent = leaf.parent_node
        if parent is None:
            continue
        parent_to_leaves.setdefault(id(parent), []).append(leaf)
        parent_to_total_children.setdefault(id(parent), parent)

    for parent_id, leaves in parent_to_leaves.items():
        if len(leaves) != 2:
            continue
        parent = parent_to_total_children[parent_id]
        if len(parent.child_nodes()) != 2:
            # Parent has additional non-leaf children; not a strict cherry.
            continue
        idx_i = taxa_metadata[leaves[0].taxon]
        idx_j = taxa_metadata[leaves[1].taxon]
        mask[idx_i, idx_j] = True
        mask[idx_j, idx_i] = True

    return mask


def compute_sigma2_separation(R_hat: np.ndarray, tree, taxa_metadata
                              ) -> Tuple[np.ndarray, np.ndarray]:
    """Return (adj, nonadj) — arrays of σ₂ values for adjacent vs. non-adjacent
    leaf pairs at step 0 of SNJ.

    Each array has length equal to the count of (i, j) pairs in the
    respective category, taking the upper triangle only (i < j).
    """
    m = R_hat.shape[0]
    sigma2 = _pairwise_sigma2_singletons(R_hat)
    cherry = cherry_mask_from_tree(tree, taxa_metadata)

    iu = np.triu_indices(m, k=1)
    sig_upper = sigma2[iu]
    cherry_upper = cherry[iu]

    adj = sig_upper[cherry_upper]
    nonadj = sig_upper[~cherry_upper]
    return adj, nonadj
