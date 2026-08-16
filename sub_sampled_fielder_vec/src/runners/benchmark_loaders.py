"""Data loaders for the operator-comparison benchmark.

Both loaders implement the notebooks' ``_load_all`` contract::

    loader(tree_id) -> (S, labels, tree, D) | None

with ``D`` permuted into ``labels`` order and ``tree`` optionally ``None`` (no
ground-truth topology available → the validity screen is skipped for that tree).

They are plain dataclasses with ``__call__`` rather than closures because
``multiprocessing`` uses the *spawn* start method on macOS: workers must be able
to unpickle the loader. ``GeneratedLoader`` is deterministic in ``tree_id``, so
workers regenerate matrices instead of receiving them over the pipe.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from ..utils.dataset_scan import FASTA_EXTS, NEWICK_EXTS
from .real_data_bpart import fasta_to_distance

LoadResult = Optional[Tuple[np.ndarray, List[str], object, np.ndarray]]


def _find(directory: Optional[Path], stem: str, exts) -> Optional[Path]:
    if directory is None:
        return None
    for ext in exts:
        p = directory / f"{stem}{ext}"
        if p.exists():
            return p
    return None


@dataclass(frozen=True)
class RealLoader:
    """FASTA (+ optional Newick) loader — ports the real notebook's ``_load_all``."""

    fasta_dir: Path
    newick_dir: Optional[Path] = None

    def group_of(self, tree_id: str) -> Optional[str]:  # no cohort grouping
        return None

    def __call__(self, tree_id: str) -> LoadResult:
        import dendropy
        import spectraltree

        fpath = _find(Path(self.fasta_dir), tree_id, FASTA_EXTS)
        if fpath is None:
            return None

        cm = dendropy.DnaCharacterMatrix.get(path=str(fpath), schema="fasta")
        obs, meta = spectraltree.charmatrix2array(cm)
        labels = [str(t.label) for t in list(meta)]
        S = spectraltree.JC_similarity_matrix(obs)

        tree = None
        npath = _find(Path(self.newick_dir) if self.newick_dir else None,
                      tree_id, NEWICK_EXTS)
        if npath is not None:
            tree = dendropy.Tree.get(path=str(npath), schema="newick",
                                     preserve_underscores=True,
                                     taxon_namespace=cm.taxon_namespace)
            tree.encode_bipartitions()

        # fasta_to_distance orders rows by taxon_namespace, charmatrix2array by
        # `meta` — realign D onto `labels` before anything downstream sees it.
        D, names = fasta_to_distance(str(fpath))
        nidx = {n: i for i, n in enumerate(names)}
        perm = np.array([nidx[l] for l in labels])
        D = D[perm][:, perm]
        return S, labels, tree, D


@dataclass(frozen=True)
class GeneratedLoader:
    """Simulated-tree loader — no files on disk, everything keyed by the tree id.

    The id carries the model, the taxon count and (optionally) the eta bin, so
    one loader serves a run that spans several models, several sizes and several
    eta bins at once. ``n_taxa`` is only the fallback for legacy ids that predate
    the ``_n`` field.

    Two routes, chosen per id:

    * no eta field — regenerate the tree and alignment from the id's seed.
    * an eta field — read the tree back from the eta pool, which is the only
      place a rejection-sampled tree exists (see ``eta_pool_bridge``).
    """

    seq_len: int
    n_taxa: Optional[int] = None
    mutation_rate: Optional[float] = None
    params: Optional[dict] = None
    seq_model: str = "JC69"

    @property
    def pop_size(self) -> float:
        return float((self.params or {}).get("pop_size", 1.0))

    def group_of(self, tree_id: str) -> Optional[str]:
        """Cohort-balancing key — one group per (model, size, eta bin)."""
        from ..utils.tree_ids import parse_spec
        try:
            return parse_spec(tree_id).cell
        except ValueError:
            return tree_id.rsplit("_", 1)[0]

    def __call__(self, tree_id: str) -> LoadResult:
        from ..models.generated_trees import make_generated
        from ..utils.tree_ids import MODEL_OF_TAG, parse_spec

        spec = parse_spec(tree_id)
        if spec.eta is None:
            return make_generated(tree_id, self.n_taxa, self.seq_len,
                                  params=self.params,
                                  mutation_rate=self.mutation_rate)

        from .eta_pool_bridge import load_from_pool

        if spec.n_taxa is None:
            raise ValueError(f"{tree_id!r} names an eta bin but no taxon count; "
                             "the pool is keyed by size")
        if self.mutation_rate is None:
            raise ValueError("eta-pooled ids need an explicit mutation_rate — it "
                             "is part of the pool's cache key")
        return load_from_pool(
            MODEL_OF_TAG.get(spec.model, spec.model), spec.n_taxa, spec.eta,
            spec.index, seq_len=self.seq_len, mu=float(self.mutation_rate),
            pop_size=self.pop_size, seq_model=self.seq_model)
