"""Turn a tree id into a simulated tree plus its similarity and distance matrices.

Generated analogue of the real-data path's ``load_S`` + ``fasta_to_distance``:
build a known tree (Kingman coalescent, birth-death, lopsided, balanced binary),
simulate Jukes-Cantor sequences on it, and return the JC similarity ``S``, the
taxon labels, the ground-truth dendropy tree and the JC distance ``D``.
Everything downstream (validity gate, screens, recovery sweeps) is data-source
agnostic and consumes these four exactly as the real-data pipeline does.

``D = similarity2distance(S) = -log(clip(S))`` reproduces ``fasta_to_distance``
(see ``src/runners/real_data_bpart.py``).

The id grammar lives in ``src/utils/tree_ids.py``. This module used to live at
``analysis/utils/generated_data.py``, which still re-exports it so the
simulation notebook keeps working; new code should import from here.
"""
from __future__ import annotations

import random
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import spectraltree

from ..utils.tree_ids import DEFAULT_TAGS, parse_spec

# --- tree builders: tag -> (m, rng, params) -> dendropy tree ------------------
# The params dict is model-specific and always optional; every builder falls back
# to the same defaults ``src/models/tree_models.py`` uses. Unlike that registry,
# these take an explicit ``rng``, which is what makes a tree reproducible from
# its id.
TREE_BUILDERS: Dict[str, Callable[[int, random.Random, dict], object]] = {
    "kingman": lambda m, rng, p: spectraltree.unrooted_pure_kingman_tree(
        m, pop_size=float(p.get("pop_size", 1.0)), rng=rng),
    "bd": lambda m, rng, p: spectraltree.unrooted_birth_death_tree(
        m, birth_rate=float(p.get("birth_rate", 0.5)),
        death_rate=float(p.get("death_rate", 0.0)), rng=rng),
    "lopsided": lambda m, rng, p: spectraltree.lopsided_tree(
        m, edge_length=float(p.get("edge_length", 1.0))),
    "balanced_binary": lambda m, rng, p: spectraltree.balanced_binary(
        m, edge_length=float(p.get("edge_length", 1.0))),
}

# Legacy alias, and the pair ``build_ids(n)`` emits when given no model list.
GENERATORS: Dict[str, Callable] = {k: TREE_BUILDERS[k] for k in DEFAULT_TAGS}

# per-generator seed offset so kingman_000 and bd_000 draw different streams
_SEED_OFFSET: Dict[str, int] = {
    "kingman": 0, "bd": 1_000_000,
    "lopsided": 2_000_000, "balanced_binary": 3_000_000,
}


def make_generated(
    gid: str,
    m: Optional[int] = None,
    seq_len: int = 1000,
    p_no_change: float = 0.9,
    *,
    params: Optional[dict] = None,
    mutation_rate: Optional[float] = None,
) -> Tuple[np.ndarray, List[str], object, np.ndarray]:
    """Generate one tree + alignment and return ``(S, labels, tree, D)``.

    Deterministic in ``gid``: re-calling with the same id rebuilds the identical
    tree and alignment, so the screen can cache lightweight results and the sweep
    can rebuild ``S`` on demand without caching the (large) matrix.

    ``m`` is only consulted when the id carries no ``_n`` field; the id always
    wins. ``mutation_rate`` defaults to ``p2t(p_no_change)``, so the default call
    is numerically identical to what it was before the knob existed.

    Ids naming an eta bin are refused: such a tree *passed* a rejection test when
    the pool was built and cannot be rebuilt from a seed (dendropy's Kingman draw
    uses its own global RNG). Those are served by
    ``src/runners/eta_pool_bridge.py``, which reads the stored tree.
    """
    spec = parse_spec(gid)
    if spec.eta is not None:
        raise ValueError(
            f"{gid!r} names an eta bin; eta-binned trees are served from the "
            "eta pool (src.runners.eta_pool_bridge.load_from_pool), not "
            "regenerated from a seed")
    if spec.model not in TREE_BUILDERS:
        raise KeyError(f"unknown generator tag {spec.model!r}; "
                       f"known: {list(TREE_BUILDERS)}")
    n = spec.n_taxa if spec.n_taxa is not None else m
    if n is None:
        raise ValueError(f"{gid!r} carries no taxon count and none was passed")

    seed = _SEED_OFFSET[spec.model] + spec.index
    tree_rng = random.Random(seed)
    seq_rng = np.random.default_rng(seed)

    tree = TREE_BUILDERS[spec.model](int(n), tree_rng, params or {})

    jc = spectraltree.Jukes_Cantor()
    mr = jc.p2t(p_no_change) if mutation_rate is None else float(mutation_rate)
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


__all__ = ["TREE_BUILDERS", "GENERATORS", "make_generated"]
