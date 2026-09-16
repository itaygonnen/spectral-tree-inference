"""Discovering real benchmark datasets: any ``<name>/{fasta,newick}`` directory pair.

A dataset holds aligned FASTA files and the true tree for each of them, matched by stem
(``random_tree_0007.fasta`` <-> ``random_tree_0007.nwk``). Cohorts live in the repo-level
``data/datasets/`` so both projects in this repo -- and the cluster copy -- point at one
place; ``$STR_DATA_DIR`` overrides it (a scratch filesystem, a shared mount). Two are on
disk today, ``1000 taxa`` and ``6000 taxa``, and anything copied in the same shape is
discovered without a code change.

This module is the *discovery* half of a real data source -- where the files are, which
ids have both halves, whether the alignments are actually aligned. Loading one tree is
:class:`src.runners.benchmark_loaders.RealLoader`'s job, and running anything on it is
:mod:`src.runners.experiment_run`'s. :meth:`Dataset.source` joins the three:

    get_dataset("6000 taxa").source()  ->  Source(name, loader, ids, m)

which is the only thing the screen, the sweep and the exporter ever see. They do not
know a FASTA file exists, which is why the same runner also drives simulated trees.

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
    """Where datasets are looked for, most specific first.

    ``$STR_DATA_DIR`` wins so a cluster can keep the data on scratch without editing code;
    the package-local path is the pre-2026-09 location, kept so an old checkout still works.
    """
    roots = []
    env = os.environ.get("STR_DATA_DIR")
    if env:
        roots.append(Path(env).expanduser())
    roots.append(_REPO / "data" / "tree_sets")     # canonical
    roots.append(_REPO / "data" / "cohorts")       # what the first installs used
    roots.append(_ROOT / "data" / "real_datasets" / "Datasets")   # pre-2026-09 layout
    return roots


DATASETS_ROOT = _REPO / "data" / "tree_sets"


@dataclass(frozen=True)
class Dataset:
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

    def probe(self) -> dict:
        """Scan the first alignment: ``{m, L_min, L_max, path}``, without parsing it.

        ``L_min != L_max`` means the file is NOT aligned -- the records have different
        lengths, which no similarity or distance matrix can be built from. Reporting only
        the first record's length (what this used to do) hid exactly that.
        """
        first = next(iter(sorted(self.fasta_dir.glob("*.fasta"))), None)
        if first is None:
            return dict(m=0, L_min=0, L_max=0, path=None)
        lengths: List[int] = []
        cur = 0
        with open(first) as fh:
            for line in fh:
                if line.startswith(">"):
                    if lengths or cur:
                        lengths.append(cur)
                    cur = 0
                else:
                    cur += len(line.strip())
        if cur:
            lengths.append(cur)
        return dict(m=len(lengths), L_min=min(lengths) if lengths else 0,
                    L_max=max(lengths) if lengths else 0, path=first)

    def shape(self) -> Tuple[int, int]:
        """``(m, L)`` of the first alignment; ``L`` is its longest record."""
        pr = self.probe()
        return pr["m"], pr["L_max"]

    def is_ragged(self) -> bool:
        """True when the first alignment's records differ in length (not an alignment)."""
        pr = self.probe()
        return pr["L_min"] != pr["L_max"]

    def loader(self):
        """The callable that turns one of this dataset's ids into its matrices."""
        from src.runners.benchmark_loaders import RealLoader
        return RealLoader(self.fasta_dir, self.newick_dir)

    def source(self, limit: Optional[int] = None):
        """This dataset as a :class:`src.runners.experiment_run.Source`."""
        from src.runners.experiment_run import Source
        m, seq_len = self.shape()
        return Source(name=self.name, loader=self.loader(),
                      ids=self.ids(limit), m=m, seq_len=seq_len)

    def load_all(self, tree_id: str):
        """``(S, labels, tree, D)`` for one tree, or ``None`` if files are missing.

        A thin pass-through to :class:`RealLoader`, kept because notebooks call it.
        """
        return self.loader()(tree_id)


def list_datasets() -> List[Dataset]:
    """Every dataset dir with both a ``fasta/`` and a ``newick/``, first root to define it."""
    out, seen = [], set()
    for root in search_roots():
        if not root.is_dir():
            continue
        for d in sorted(root.iterdir()):
            if d.name in seen:
                continue
            if (d / "fasta").is_dir() and (d / "newick").is_dir():
                out.append(Dataset(d.name, root))
                seen.add(d.name)
    return out


def _key(name: str) -> str:
    """Compare dataset names ignoring case, spaces and underscores."""
    return "".join(str(name).lower().split()).replace("_", "").replace("-", "")


def get_dataset(name: str) -> Dataset:
    """Find a dataset by name, forgiving case and spacing ("6000 taxa" == "6000_Taxa")."""
    datasets = list_datasets()
    for c in datasets:
        if c.name == name:
            return c
    for c in datasets:                       # then the forgiving match
        if _key(c.name) == _key(name):
            return c
    raise FileNotFoundError(
        f"no dataset {name!r}.\n" + describe_search())


def describe_search() -> str:
    """Where datasets were looked for, what was found, and why a folder was rejected.

    A dataset is a directory with BOTH a fasta/ and a newick/ subdirectory; the usual
    failure is data unpacked one level too deep, or the two halves side by side under
    different names.
    """
    lines = []
    for root in search_roots():
        if not root.is_dir():
            lines.append(f"  {root}  (does not exist)")
            continue
        entries = sorted(d for d in root.iterdir() if d.is_dir())
        if not entries:
            lines.append(f"  {root}  (empty)")
            continue
        lines.append(f"  {root}")
        for d in entries:
            has_f, has_n = (d / "fasta").is_dir(), (d / "newick").is_dir()
            if has_f and has_n:
                lines.append(f"    ok       {d.name!r}")
            else:
                missing = ", ".join(x for x, ok in (("fasta/", has_f),
                                                    ("newick/", has_n)) if not ok)
                inner = ", ".join(sorted(x.name for x in d.iterdir() if x.is_dir())[:6])
                lines.append(f"    skipped  {d.name!r}: no {missing}"
                             + (f" (contains: {inner})" if inner else ""))
    return ("a dataset is <root>/<name>/fasta/*.fasta beside <name>/newick/*.nwk\n"
            "looked in:\n" + "\n".join(lines))


# Where results and resumable state live is not a property of the data source, so it
# lives in src/utils/run_paths.py now. Re-exported because notebooks and scripts import
# these names from here.
from src.utils.run_paths import (cache_dir, new_run_dir, results_root,   # noqa: E402,F401
                                 screen_cache_path, sweep_cache_dir)
