"""Discover real-data tree datasets (FASTA alignments + optional Newick trees).

The benchmark only *needs* the FASTA files: ``S`` and ``D`` are both derived from
the alignment, and every downstream metric (recovery curve, scale plot) is
measured against the full-matrix reference. Newick trees are optional — they feed
the validity screen (`is this bipartition a real single-edge clan of the true
tree?`) and nothing else.

Directory layout in the wild is not tidy. On the reference machine the complete
600-alignment set and the complete 600-tree set live under *different* roots::

    data/Datasets/Fasta 1000 taxa/extracted/fasta/   600 .fasta   (no sibling newick/)
    sub_sampled_fielder_vec/data/real_datasets/Datasets/1000 taxa/newick/   600 .nwk

so pairing is done on a *normalized size label* across all roots rather than on
directory siblinghood.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

FASTA_EXTS = (".fasta", ".fa", ".fna")
NEWICK_EXTS = (".nwk", ".newick", ".tree")

# Directory names that describe packaging, not a dataset — skipped when walking up
# from a fasta/newick dir to find the size label.
_STRUCTURAL = {"extracted", "datasets", "dataset", "data", "real_datasets",
               "fasta", "newick", "trees", "alignments"}

_MAX_DEPTH = 5


def _normalize_label(name: str) -> str:
    """``'Fasta 1000 taxa'`` and ``'1000 taxa'`` both -> ``'1000 taxa'``."""
    s = re.sub(r"^(fasta|newick|alignments?|trees?)\b[\s_-]*", "", name.strip(),
               flags=re.IGNORECASE)
    return re.sub(r"[\s_]+", " ", s).strip().lower() or name.strip().lower()


def _label_for(d: Path, root: Path) -> str:
    """Nearest non-structural ancestor of ``d``, normalized."""
    for parent in d.parents:
        if parent == root.parent:
            break
        if parent.name.lower() not in _STRUCTURAL:
            return _normalize_label(parent.name)
        if parent == root:
            break
    return _normalize_label(d.parent.name)


def _stems(d: Path, exts: Sequence[str]) -> List[str]:
    return sorted({p.stem for p in d.iterdir()
                   if p.is_file() and p.suffix.lower() in exts})


@dataclass
class DatasetEntry:
    """One runnable (or not-yet-runnable) dataset, keyed by its size label."""

    label: str
    fasta_dir: Optional[Path] = None
    newick_dir: Optional[Path] = None
    fasta_stems: List[str] = field(default_factory=list)
    newick_stems: List[str] = field(default_factory=list)

    @property
    def tree_ids(self) -> List[str]:
        """Ids the pipeline can actually run: fasta stems, ∩ newick when present."""
        if self.newick_stems:
            return sorted(set(self.fasta_stems) & set(self.newick_stems))
        return list(self.fasta_stems)

    @property
    def runnable(self) -> bool:
        return bool(self.tree_ids)

    @property
    def status(self) -> str:
        if not self.fasta_stems:
            return "no fasta — not runnable"
        if not self.newick_stems:
            return "no newick — validity screen will be skipped"
        missing = len(set(self.fasta_stems) - set(self.newick_stems))
        if missing:
            return f"partial — {missing} fasta without a matching newick"
        return "complete"

    def describe(self) -> str:
        return (f"{self.label:<16} — {len(self.fasta_stems):>5} fasta, "
                f"{len(self.newick_stems):>5} newick   [{self.status}]")


def scan_roots(roots: Sequence[Path]) -> List[DatasetEntry]:
    """Find every fasta/newick directory under ``roots`` and pair them by label.

    Directories holding zero matching files are ignored. When two directories map
    to the same label (e.g. a partial download plus the full extract), the one
    with more files wins.
    """
    entries: Dict[str, DatasetEntry] = {}

    for root in roots:
        root = Path(root)
        if not root.is_dir():
            continue
        for d in _walk_dirs(root):
            kind = d.name.lower()
            if kind not in ("fasta", "newick"):
                continue
            exts = FASTA_EXTS if kind == "fasta" else NEWICK_EXTS
            try:
                stems = _stems(d, exts)
            except OSError:
                continue
            if not stems:
                continue
            label = _label_for(d, root)
            e = entries.setdefault(label, DatasetEntry(label=label))
            if kind == "fasta":
                if len(stems) > len(e.fasta_stems):
                    e.fasta_dir, e.fasta_stems = d, stems
            else:
                if len(stems) > len(e.newick_stems):
                    e.newick_dir, e.newick_stems = d, stems

    return sorted(entries.values(),
                  key=lambda e: (not e.runnable, -len(e.tree_ids), e.label))


def _walk_dirs(root: Path):
    """Yield directories under ``root`` down to ``_MAX_DEPTH``, skipping dotdirs."""
    stack = [(root, 0)]
    while stack:
        d, depth = stack.pop()
        if depth > _MAX_DEPTH:
            continue
        try:
            children = [c for c in d.iterdir()
                        if c.is_dir() and not c.name.startswith(".")]
        except OSError:
            continue
        for c in children:
            yield c
            stack.append((c, depth + 1))


def entry_from_paths(fasta_dir: Path, newick_dir: Optional[Path] = None,
                     label: str = "custom") -> DatasetEntry:
    """Build an entry from explicitly supplied paths (the 'type a path' branch)."""
    fasta_dir = Path(fasta_dir)
    if not fasta_dir.is_dir():
        raise NotADirectoryError(f"fasta dir not found: {fasta_dir}")
    e = DatasetEntry(label=label, fasta_dir=fasta_dir,
                     fasta_stems=_stems(fasta_dir, FASTA_EXTS))
    if newick_dir:
        newick_dir = Path(newick_dir)
        if not newick_dir.is_dir():
            raise NotADirectoryError(f"newick dir not found: {newick_dir}")
        e.newick_dir = newick_dir
        e.newick_stems = _stems(newick_dir, NEWICK_EXTS)
    return e
