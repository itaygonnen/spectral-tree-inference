"""The generated-tree id grammar.

A tree id is a *recipe*, not a file on disk: it names the model, optionally the
taxon count and an eta bin, and carries an index that seeds the RNG::

    kingman_007                 legacy — taxon count comes from the caller
    kingman_n1000_007           1000 taxa
    kingman_n1000_eta05_003     1000 taxa, eta bin 5

Pure string handling — no numpy, no spectraltree, no cache. The generator that
turns an id into a tree lives in ``src/models/generated_trees.py``; the eta-binned
ids are served by ``src/runners/eta_pool_bridge.py``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

# The two default generators, in this order: ``build_ids(n)`` with no explicit
# model list emits exactly these, which is what the supporting notebook's
# SCREEN_IDS expects. Do not reorder.
DEFAULT_TAGS: Tuple[str, ...] = ("kingman", "bd")

# Tag as it appears in an id -> model name in ``src/models/tree_models.py`` and in
# the eta pool's cache key. Only 'bd' differs; kept so ids already on disk stay
# readable.
MODEL_OF_TAG = {
    "kingman": "kingman", "bd": "birth_death",
    "lopsided": "lopsided", "balanced_binary": "balanced_binary",
}

# model matched lazily, so 'kingman_mean_n500_003' keeps its underscored tag
_ID_RE = re.compile(
    r"^(?P<model>.+?)(?:_n(?P<n>\d+))?(?:_eta(?P<eta>\d+))?_(?P<idx>\d+)$")


@dataclass(frozen=True)
class GenSpec:
    """One parsed tree id."""

    model: str
    index: int
    n_taxa: Optional[int] = None
    eta: Optional[int] = None

    @property
    def cell(self) -> str:
        """Grouping key for cohort balancing — one cell per (model, n, eta)."""
        parts = [self.model]
        if self.n_taxa is not None:
            parts.append(f"n{self.n_taxa}")
        if self.eta is not None:
            parts.append(f"eta{self.eta:02d}")
        return "|".join(parts)


def make_tree_id(model: str, index: int, n_taxa: Optional[int] = None,
                 eta: Optional[int] = None) -> str:
    """``('kingman', 3, 1000, 5)`` -> ``'kingman_n1000_eta05_003'``."""
    mid = model
    if n_taxa is not None:
        mid += f"_n{int(n_taxa)}"
    if eta is not None:
        mid += f"_eta{int(eta):02d}"
    return f"{mid}_{int(index):03d}"


def parse_spec(gid: str) -> GenSpec:
    """Parse a tree id into its fields; raises ``ValueError`` on a malformed id."""
    m = _ID_RE.match(gid)
    if m is None:
        raise ValueError(f"malformed tree id {gid!r}")
    return GenSpec(
        model=m.group("model"),
        index=int(m.group("idx")),
        n_taxa=int(m.group("n")) if m.group("n") else None,
        eta=int(m.group("eta")) if m.group("eta") else None,
    )


def parse_id(gid: str) -> Tuple[str, int]:
    """Legacy two-field view: ``'kingman_007'`` -> ``('kingman', 7)``."""
    spec = parse_spec(gid)
    return spec.model, spec.index


def build_ids(
    n_per_gen: Optional[int] = None,
    *,
    models: Optional[Sequence[str]] = None,
    n_values: Optional[Sequence[int]] = None,
    etas: Optional[Sequence[int]] = None,
    per_cell: Optional[int] = None,
) -> List[str]:
    """Tree ids for one run.

    Legacy form — ``build_ids(300)`` -> ``['kingman_000', ..., 'bd_299']``, the
    two default generators, no size or eta in the id.

    Grid form — every combination of ``models`` x ``n_values`` x ``etas`` is one
    *cell*, and each cell gets ``per_cell`` ids. Ids come out **index-major**, so
    any prefix of the list is balanced across cells: an interrupted screen has
    covered every cell evenly rather than finishing the first one.
    """
    grid = any(x is not None for x in (models, n_values, etas, per_cell))
    if not grid:
        if n_per_gen is None:
            raise TypeError("build_ids needs n_per_gen, or the keyword grid form")
        return [f"{tag}_{i:03d}"
                for tag in DEFAULT_TAGS for i in range(int(n_per_gen))]

    count = per_cell if per_cell is not None else n_per_gen
    if count is None:
        raise TypeError("the grid form of build_ids needs per_cell")
    cells = [(mo, n, e)
             for mo in list(models or DEFAULT_TAGS)
             for n in list(n_values or [None])
             for e in list(etas or [None])]
    return [make_tree_id(mo, i, n, e)
            for i in range(int(count)) for mo, n, e in cells]


__all__ = ["DEFAULT_TAGS", "MODEL_OF_TAG", "GenSpec", "make_tree_id",
           "parse_spec", "parse_id", "build_ids"]
