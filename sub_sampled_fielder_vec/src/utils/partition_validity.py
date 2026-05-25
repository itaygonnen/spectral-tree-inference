"""Partition-validity check against a ground-truth dendropy tree.

A partition (boolean mask of length n over the tree's leaves) is "valid"
iff it corresponds to a single edge cut of the unrooted tree — i.e. there
exists one edge whose removal separates exactly the True-side leaves from
the False-side leaves.

Wraps ``spectraltree.utils.check_is_bipartition``. Handles the TaxaMetadata
construction in leaf-iteration order (which is the same order used by
``simulate_sequences`` / ``_get_cached_similarity_matrix`` when building M).
"""
from __future__ import annotations

import numpy as np

from spectraltree.utils import TaxaMetadata, check_is_bipartition


def check_partition_valid_in_tree(tree, partition_mask: np.ndarray) -> bool:
    """True iff ``partition_mask`` is a single-edge bipartition of ``tree``."""
    if partition_mask is None:
        return False
    tree.encode_bipartitions()
    leaf_taxa = [leaf.taxon for leaf in tree.leaf_nodes()]
    meta = TaxaMetadata(tree.taxon_namespace, leaf_taxa)
    return bool(check_is_bipartition(tree, partition_mask, meta))
