"""Recursive Fiedler-vector bipartitioning (STDR partition phase only, no merge).

Replicates the top-down splitting in
``spectraltree/spectral_tree_reconstruction.py::deep_spectral_tree_reconstruction``
(lines 614-694) without the bottom-up ``margeTreesLeftRight`` phase.
"""
from __future__ import annotations

from typing import FrozenSet, List, Tuple

import numpy as np
import scipy.linalg

from spectraltree.spectral_tree_reconstruction import partition_taxa


def recursive_split(
    similarity: np.ndarray,
    threshold: int,
    num_gaps: int = 0,
    min_split: int = 1,
) -> Tuple[np.ndarray, List[FrozenSet[int]]]:
    """Run the recursive Fiedler split once; return cluster ids and bipartitions.

    At each internal split on subset ``I`` we emit the bipartition
    ``frozenset(I_L)`` of the *full* taxon set, canonicalised to the smaller
    side (tiebreak: side containing the smallest taxon index). Recursion
    stops when ``|I| <= threshold``; that ``I`` becomes one leaf cluster
    (no bipartition emitted for it).
    """
    n = similarity.shape[0]
    cluster_ids = np.full(n, -1, dtype=np.int64)
    bipartitions: List[FrozenSet[int]] = []
    stack: list[np.ndarray] = [np.arange(n, dtype=np.int64)]
    next_id = 0

    while stack:
        I = stack.pop()
        if I.size <= threshold:
            cluster_ids[I] = next_id
            next_id += 1
            continue

        sub = similarity[np.ix_(I, I)]
        L = np.diag(sub.sum(axis=0)) - sub
        _, V = scipy.linalg.eigh(L, subset_by_index=(0, 1))
        fiedler = V[:, 1]

        try:
            mask = partition_taxa(fiedler, sub, num_gaps, min_split)
        except Exception:
            cluster_ids[I] = next_id
            next_id += 1
            continue

        left, right = I[mask], I[~mask]
        if left.size == 0 or right.size == 0:
            cluster_ids[I] = next_id
            next_id += 1
            continue

        # Canonicalise the bipartition by smaller side; tiebreak by min index.
        a, b = frozenset(left.tolist()), frozenset(right.tolist())
        if len(a) < len(b) or (len(a) == len(b) and min(a) < min(b)):
            bipartitions.append(a)
        else:
            bipartitions.append(b)

        stack.append(left)
        stack.append(right)

    return cluster_ids, bipartitions


def recursive_partition(
    similarity: np.ndarray,
    threshold: int,
    num_gaps: int = 0,
    min_split: int = 1,
) -> np.ndarray:
    """Recursively bipartition taxa using Fiedler vectors; return leaf cluster ids."""
    cluster_ids, _ = recursive_split(similarity, threshold, num_gaps, min_split)
    return cluster_ids


def recursive_bipartitions(
    similarity: np.ndarray,
    threshold: int,
    num_gaps: int = 0,
    min_split: int = 1,
) -> List[FrozenSet[int]]:
    """Return the set of bipartitions produced by the recursive Fiedler split.

    Each entry is the smaller side of an internal split (canonicalised by the
    smallest taxon index on ties), expressed in the *full* taxon index space.
    """
    _, bipartitions = recursive_split(similarity, threshold, num_gaps, min_split)
    return bipartitions
