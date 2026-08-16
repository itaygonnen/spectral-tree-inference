"""Serve eta-binned trees to the benchmark from the eta pool.

The benchmark's plain generated trees are rebuilt from a seed on demand. An
eta-binned tree cannot be: it is a tree that *passed* a rejection test at
pool-build time (``|eta - target| <= eta_tol`` plus a valid-bipartition check),
and dendropy's Kingman draw uses its own global RNG, so the seed alone does not
reproduce it. The pool therefore stores the tree, and this module reads it back.

One pool sample supplies everything a benchmark loader must return::

    S      = M                       the stored similarity matrix
    labels = tree leaf order         M's rows are in tree leaf-iteration order
    tree   = parsed stored newick    ground truth for the validity screen
    D      = similarity2distance(S)  same transform the real-data path uses

That row-order claim is the one load-bearing assumption here; it is the same one
``src/utils/partition_validity.py`` documents, and it is checked directly by
:func:`verify_pool_alignment`.

``scripts/build_eta_pool.py`` stays the canonical pool builder — :func:`ensure_cells`
shells out to it rather than reimplementing the accept/reject loop, so there is
only ever one definition of what lands in a bin. Note that build_eta_pool
rewrites ``samples_per_bin`` in the pool manifest to whatever it was invoked
with; bin counts are always recounted from disk, so nothing is lost.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cache_io import CACHE_ROOT
from ..utils.eta_pool_cache import (
    list_completed_samples, load_pool_entry, param_key, sample_dir,
)

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_BUILDER = _PROJECT_ROOT / "scripts" / "build_eta_pool.py"

# What ``build_eta_pool.py --tree-model`` accepts. A model outside this set can
# still be run by the benchmark, just not in eta-pooled mode.
POOLABLE_MODELS = ("kingman", "kingman_mean", "lopsided", "birth_death",
                   "balanced_binary")


def pool_key(model: str, n_taxa: int, seq_len: int, mu: float,
             pop_size: float, seq_model: str = "JC69") -> str:
    """The pool's ``param_key`` for one (model, size) configuration."""
    return param_key(n=int(n_taxa), seq_len=int(seq_len), mu=float(mu),
                     tree_model=str(model), pop_size=float(pop_size),
                     seq_model=str(seq_model))


def available_samples(model: str, n_taxa: int, eta: int, *, seq_len: int,
                      mu: float, pop_size: float, seq_model: str = "JC69",
                      cache_root: Path = CACHE_ROOT) -> List[int]:
    """Completed sample indices in one bin that carry a stored tree.

    Samples predating the tree-persistence change are skipped: without the tree
    there is no ground truth, so they would be screened but never make the
    both-valid cohort.
    """
    key = pool_key(model, n_taxa, seq_len, mu, pop_size, seq_model)
    out = []
    for idx in list_completed_samples(cache_root, key, int(eta)):
        d = sample_dir(cache_root, key, int(eta), idx)
        if (d / "tree.txt").exists():
            out.append(idx)
    return out


def load_from_pool(model: str, n_taxa: int, eta: int, slot: int, *,
                   seq_len: int, mu: float, pop_size: float,
                   seq_model: str = "JC69", cache_root: Path = CACHE_ROOT):
    """Return ``(S, labels, tree, D)`` for the ``slot``-th sample of one bin.

    ``slot`` is a position in the sorted list of available samples, not a pool
    sample index — that is what makes a tree id stable while the pool grows:
    new samples get higher indices and land after the existing slots.
    """
    import dendropy
    import spectraltree

    idxs = available_samples(model, n_taxa, eta, seq_len=seq_len, mu=mu,
                             pop_size=pop_size, seq_model=seq_model,
                             cache_root=cache_root)
    if slot >= len(idxs):
        return None

    key = pool_key(model, n_taxa, seq_len, mu, pop_size, seq_model)
    entry = load_pool_entry(cache_root, key, int(eta), idxs[slot])
    if entry is None:
        return None
    M, _v_pop, _partition, tree_newick, _meta = entry
    if tree_newick is None:
        return None

    tree = dendropy.Tree.get(data=tree_newick, schema="newick",
                             preserve_underscores=True)
    tree.encode_bipartitions()
    labels = [str(leaf.taxon.label) for leaf in tree.leaf_nodes()]

    S = np.asarray(M, float)
    D = spectraltree.similarities.similarity2distance(S)
    np.fill_diagonal(D, 0.0)
    return S, labels, tree, D


def verify_pool_alignment(model: str, n_taxa: int, eta: int, slot: int = 0,
                          **kw) -> Optional[bool]:
    """Check the row-order assumption on one sample.

    The stored partition is in M-row order; if that really is tree leaf order,
    it must be a single-edge bipartition of the re-parsed tree. ``None`` when
    there is nothing on disk to check.
    """
    from ..utils.partition_validity import check_partition_valid_in_tree

    cache_root = kw.pop("cache_root", CACHE_ROOT)
    idxs = available_samples(model, n_taxa, eta, cache_root=cache_root, **kw)
    if slot >= len(idxs):
        return None
    key = pool_key(model, n_taxa, kw["seq_len"], kw["mu"], kw["pop_size"],
                   kw.get("seq_model", "JC69"))
    entry = load_pool_entry(cache_root, key, int(eta), idxs[slot])
    if entry is None or entry[2] is None or entry[3] is None:
        return None
    loaded = load_from_pool(model, n_taxa, eta, slot, cache_root=cache_root, **kw)
    if loaded is None:
        return None
    return bool(check_partition_valid_in_tree(loaded[2], entry[2]))


def counts_by_cell(models: Sequence[str], n_values: Sequence[int],
                   etas: Sequence[int], *, seq_len: int, mu: float,
                   pop_size: float, seq_model: str = "JC69",
                   cache_root: Path = CACHE_ROOT) -> Dict[Tuple[str, int, int], int]:
    """How many usable samples each (model, n, eta) cell already has."""
    return {
        (mo, int(n), int(e)): len(available_samples(
            mo, int(n), int(e), seq_len=seq_len, mu=mu, pop_size=pop_size,
            seq_model=seq_model, cache_root=cache_root))
        for mo in models for n in n_values for e in etas
    }


def ensure_cells(models: Sequence[str], n_values: Sequence[int],
                 etas: Sequence[int], per_cell: int, *, seq_len: int, mu: float,
                 pop_size: float, seq_model: str = "JC69",
                 max_attempts: int = 2000, eta_tol: float = 1.0,
                 cache_root: Path = CACHE_ROOT) -> Dict[Tuple[str, int], int]:
    """Top the pool up to ``per_cell`` samples per bin, one builder run per size.

    Returns the builder's exit code per (model, n). Bins that are already full
    cost one process start and nothing else — the builder recounts from disk and
    returns immediately. Output is streamed, not captured: a rejection-sampling
    run for a high eta at a large n can take a long time and its progress log is
    the only signal that it is working.
    """
    codes: Dict[Tuple[str, int], int] = {}
    for model in models:
        for n in n_values:
            have = counts_by_cell([model], [n], etas, seq_len=seq_len, mu=mu,
                                  pop_size=pop_size, seq_model=seq_model,
                                  cache_root=cache_root)
            short = sorted(int(e) for (_m, _n, e), c in have.items()
                           if c < per_cell)
            if not short:
                codes[(model, int(n))] = 0
                continue
            cmd = [sys.executable, str(_BUILDER),
                   "--n", str(int(n)),
                   "--mu", repr(float(mu)),
                   "--pop-size", repr(float(pop_size)),
                   "--seq-len", str(int(seq_len)),
                   "--tree-model", str(model),
                   "--seq-model", str(seq_model),
                   "--eta-tol", repr(float(eta_tol)),
                   "--samples-per-bin", str(int(per_cell)),
                   "--max-attempts", str(int(max_attempts)),
                   "--eta-targets", *[str(e) for e in short]]
            print(f"\n  building pool: {model} n={n} bins={short} "
                  f"(target {per_cell}/bin, ≤{max_attempts} attempts)")
            codes[(model, int(n))] = subprocess.call(cmd, cwd=str(_PROJECT_ROOT))
    return codes


__all__ = ["POOLABLE_MODELS", "pool_key", "available_samples", "load_from_pool",
           "verify_pool_alignment", "counts_by_cell", "ensure_cells"]
