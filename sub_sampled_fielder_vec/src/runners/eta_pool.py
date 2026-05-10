"""Single-attempt builder for the eta-binned matrix pool.

One call to ``attempt_one`` produces a fresh (M, v_pop) pair for one Kingman
tree at a given seed. The script that owns the pool decides whether to keep
or discard the result based on the realized eta.

Kept deliberately small — no I/O, no filtering — so it stays useful as a
building block for future pool variants (different tree models, different
sequence models, etc.).
"""
from __future__ import annotations

from typing import Tuple

import numpy as np

from ..models.tree_models import get_tree_factory
from ..models.sequence_models import get_sequence_factory
from ..core.utils import generate_sequences
from ..utils.random_entries import (
    _get_cached_similarity_matrix,
    compute_fiedler_from_similarity,
)


def attempt_one(
    seed: int,
    n: int,
    mu: float,
    pop_size: float,
    seq_len: int,
    tree_model: str = "kingman",
    seq_model: str = "JC69",
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate one tree, simulate sequences, return (M, v_pop).

    Args:
        seed: RNG seed; sets ``np.random.seed`` for reproducibility.
        n: number of taxa.
        mu: mutation rate (passed to ``simulate_sequences`` via ``generate_sequences``).
        pop_size: Kingman pop_size parameter.
        seq_len: simulated sequence length.
        tree_model: registered tree-model name (default ``"kingman"``).
        seq_model: registered sequence-model name (default ``"JC69"``).

    Returns:
        (M, v_pop): n x n similarity matrix and its reference Fiedler vector.
    """
    np.random.seed(int(seed))

    tree_factory = get_tree_factory(
        tree_model, {"num_taxa": int(n), "pop_size": float(pop_size)}
    )
    seq_factory = get_sequence_factory(seq_model, {"mutation_rate": float(mu)})

    tree = tree_factory()
    seq_obj = seq_factory()
    observations = generate_sequences(int(n), int(seq_len), float(mu), tree, seq_obj)

    M = _get_cached_similarity_matrix(observations)
    v_pop = compute_fiedler_from_similarity(M)
    return M, v_pop
