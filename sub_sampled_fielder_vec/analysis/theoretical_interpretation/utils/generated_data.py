"""Generated-tree data loader for the simulation distance-vs-similarity notebook.

Generated analogue of the real-data notebook's ``load_S`` + ``fasta_to_distance``:
build a known tree (Kingman coalescent or birth-death), simulate Jukes-Cantor
sequences on it, and return the JC similarity ``S``, the JC distance ``D``, the
ground-truth dendropy tree, and the taxon labels. Everything downstream (validity
gate, screens, recovery sweeps) is data-source agnostic and reuses these outputs
exactly as the real-data pipeline does.

``D = similarity2distance(S) = -log(clip(S))`` reproduces ``fasta_to_distance``
(see ``src/runners/real_data_bpart.py``).
"""
from __future__ import annotations

import random
from typing import Callable, Dict, List, Tuple

import numpy as np
import spectraltree

# --- generator registry: tag -> (m, random.Random) -> dendropy tree ----------
GENERATORS: Dict[str, Callable[[int, random.Random], object]] = {
    "kingman": lambda m, rng: spectraltree.unrooted_pure_kingman_tree(m, rng=rng),
    "bd": lambda m, rng: spectraltree.unrooted_birth_death_tree(m, rng=rng),
}

# per-generator seed offset so kingman_000 and bd_000 draw different streams
_SEED_OFFSET: Dict[str, int] = {"kingman": 0, "bd": 1_000_000}


def build_ids(n_per_gen: int) -> List[str]:
    """Parallel of the real-data ``SCREEN_IDS`` list: ``['kingman_000', ..., 'bd_000', ...]``."""
    return [f"{tag}_{i:03d}" for tag in GENERATORS for i in range(n_per_gen)]


def parse_id(gid: str) -> Tuple[str, int]:
    """Split ``'kingman_007'`` -> ``('kingman', 7)``."""
    tag, idx = gid.rsplit("_", 1)
    return tag, int(idx)


def make_generated(
    gid: str, m: int, seq_len: int, p_no_change: float = 0.9
) -> Tuple[np.ndarray, List[str], object, np.ndarray]:
    """Generate one tree + sequences and return ``(S, labels, tree, D)``.

    Deterministic in ``gid``: re-calling with the same id rebuilds the identical
    tree and alignment, so the screen can cache lightweight results and the sweep
    cell can rebuild ``S`` on demand without caching the (large) matrix.
    """
    tag, idx = parse_id(gid)
    if tag not in GENERATORS:
        raise KeyError(f"unknown generator tag {tag!r}; known: {list(GENERATORS)}")
    seed = _SEED_OFFSET[tag] + idx
    tree_rng = random.Random(seed)
    seq_rng = np.random.default_rng(seed)

    tree = GENERATORS[tag](m, tree_rng)

    jc = spectraltree.Jukes_Cantor()
    mr = jc.p2t(p_no_change)
    obs, meta = spectraltree.simulate_sequences(
        seq_len, tree_model=tree, seq_model=jc, mutation_rate=mr,
        rng=seq_rng, alphabet="DNA",
    )
    labels = [str(t.label) for t in list(meta)]

    S = spectraltree.JC_similarity_matrix(obs)
    D = spectraltree.similarities.similarity2distance(S)
    np.fill_diagonal(D, 0.0)

    tree.encode_bipartitions()
    return S, labels, tree, D
