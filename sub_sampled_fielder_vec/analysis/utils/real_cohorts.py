"""Real benchmark cohorts: any ``<name>/{fasta,newick}`` directory pair.

A cohort holds aligned FASTA files and the true tree for each of them, matched by stem
(``random_tree_0007.fasta`` <-> ``random_tree_0007.nwk``). Cohorts live in the repo-level
``data/cohorts/`` so both projects in this repo -- and the cluster copy -- point at one
place; ``$STR_DATA_DIR`` overrides it (a scratch filesystem, a shared mount). Two are on
disk today, ``1000 taxa`` and ``6000 taxa``, and anything copied in the same shape is
discovered without a code change.

One load per tree gives every consumer both operators' inputs:

    S = JC_similarity_matrix(obs)     -> the similarity route,  L(S)
    D = fasta_to_distance(path)       -> the distance route,    B = H D H

``D`` is permuted into ``labels`` order so a single leaf-index map serves both.

Cost scales as O(m^2 L) for the two matrix builds and O(m^3) for every eigensolve; at
m=6000 that is ~1 min to load a tree and ~9 s per sub-sampled Fiedler vector, so cache
anything computed on top of a load.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]          # sub_sampled_fielder_vec
_REPO = _ROOT.parent                                 # spectral-tree-inference
for _p in (str(_ROOT), str(_REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

def search_roots() -> List[Path]:
    """Where cohorts are looked for, most specific first.

    ``$STR_DATA_DIR`` wins so a cluster can keep the data on scratch without editing code;
    the package-local path is the pre-2026-09 location, kept so an old checkout still works.
    """
    roots = []
    env = os.environ.get("STR_DATA_DIR")
    if env:
        roots.append(Path(env).expanduser())
    roots.append(_REPO / "data" / "cohorts")
    roots.append(_ROOT / "data" / "real_datasets" / "Datasets")
    return roots


DATASETS_ROOT = _REPO / "data" / "cohorts"


@dataclass(frozen=True)
class Cohort:
    """One real dataset: ``<root>/<name>/{fasta,newick}``."""

    name: str
    root: Optional[Path] = None

    @property
    def dir(self) -> Path:
        if self.root is not None:
            return Path(self.root) / self.name
        for r in search_roots():
            if (r / self.name / "fasta").is_dir():
                return r / self.name
        return DATASETS_ROOT / self.name

    @property
    def fasta_dir(self) -> Path:
        return self.dir / "fasta"

    @property
    def newick_dir(self) -> Path:
        return self.dir / "newick"

    def ids(self, limit: Optional[int] = None) -> List[str]:
        """Tree ids with both an alignment and a true tree on disk, sorted."""
        ids = sorted(p.stem for p in self.fasta_dir.glob("*.fasta")
                     if (self.newick_dir / f"{p.stem}.nwk").exists())
        return ids[:limit] if limit else ids

    def shape(self) -> Tuple[int, int]:
        """``(m, L)`` read from the first alignment's header block, without parsing it."""
        first = next(iter(sorted(self.fasta_dir.glob("*.fasta"))), None)
        if first is None:
            return (0, 0)
        m, seq_len, cur = 0, 0, 0
        with open(first) as fh:
            for line in fh:
                if line.startswith(">"):
                    m += 1
                    if m == 2:
                        seq_len = cur
                elif m == 1:
                    cur += len(line.strip())
        return m, seq_len or cur

    def load_all(self, tree_id: str):
        """Return ``(S, labels, tree, D)`` for one tree, or ``None`` if files are missing.

        ``D`` is aligned to ``labels`` order. ``tree`` is read into the alignment's taxon
        namespace with bipartitions encoded, ready for the validity gate.
        """
        import dendropy
        import spectraltree
        from src.runners.real_data_bpart import fasta_to_distance

        fpath = self.fasta_dir / f"{tree_id}.fasta"
        tpath = self.newick_dir / f"{tree_id}.nwk"
        if not (fpath.exists() and tpath.exists()):
            return None

        cm = dendropy.DnaCharacterMatrix.get(path=str(fpath), schema="fasta")
        obs, meta = spectraltree.charmatrix2array(cm)
        labels = [str(t.label) for t in list(meta)]
        S = spectraltree.JC_similarity_matrix(obs)

        D, names = fasta_to_distance(str(fpath))
        nidx = {n: i for i, n in enumerate(names)}
        perm = np.array([nidx[l] for l in labels])
        D = D[perm][:, perm]

        tree = dendropy.Tree.get(path=str(tpath), schema="newick",
                                 preserve_underscores=True,
                                 taxon_namespace=cm.taxon_namespace)
        tree.encode_bipartitions()
        return S, labels, tree, D


def list_cohorts() -> List[Cohort]:
    """Every dataset dir with both a ``fasta/`` and a ``newick/``, first root to define it."""
    out, seen = [], set()
    for root in search_roots():
        if not root.is_dir():
            continue
        for d in sorted(root.iterdir()):
            if d.name in seen:
                continue
            if (d / "fasta").is_dir() and (d / "newick").is_dir():
                out.append(Cohort(d.name, root))
                seen.add(d.name)
    return out


def get_cohort(name: str) -> Cohort:
    for c in list_cohorts():
        if c.name == name:
            return c
    raise FileNotFoundError(
        f"no cohort {name!r} in {[str(r) for r in search_roots()]} "
        f"(have: {[x.name for x in list_cohorts()]})")


CACHE_ROOT = (_ROOT / "analysis" / "notebooks_cache" / "distance_vs_similarity_real")

# The m=6000 caches were written before cohorts were a parameter; keep their paths so the
# 100-tree screen already on disk is not orphaned by the rename.
_LEGACY = {"6000 taxa": ("eta_screen_6000.npz", "sweep6000")}


def _slug(name: str) -> str:
    return name.replace(" ", "_").lower()


def screen_cache_path(name: str) -> Path:
    legacy = _LEGACY.get(name)
    return CACHE_ROOT / (legacy[0] if legacy else f"eta_screen_{_slug(name)}.npz")


def sweep_cache_dir(name: str) -> Path:
    legacy = _LEGACY.get(name)
    return CACHE_ROOT / (legacy[1] if legacy else f"sweep_{_slug(name)}")
