"""Does a dataset satisfy the contract ``RealLoader`` reads it under?

Five ways a directory of alignments fails the runners, checked in one pass: a FASTA that
is not rectangular, a gap character where none is expected, a stem with no partner, a
taxon label that differs between the alignment and the tree, and an alphabet that makes
``charmatrix2array`` report k != 4.

The label check is the one that bites. ``RealLoader`` parses the newick into the
alignment's taxon namespace (``benchmark_loaders.py:77``), so a header that does not match
a leaf label exactly is a hard failure rather than a warning.

Generic on purpose -- it takes a directory, not a generator -- so any dataset can be held
to it, not only the simulated ones.
"""
from __future__ import annotations

from pathlib import Path
from typing import List

ACGT = set("ACGT")

def _scan_fasta(path: Path) -> tuple:
    """``(headers, lengths, alphabet)`` without parsing into a matrix."""
    headers, lengths, alphabet = [], [], set()
    cur = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].strip())
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                s = line.strip()
                cur += len(s)
                alphabet.update(s)
    if cur:
        lengths.append(cur)
    return headers, lengths, alphabet


def check_dataset(dataset: Path, expect_gapless: bool = True) -> List[str]:
    """Every way this dataset could fail ``RealLoader``, as a list of messages."""
    import dendropy

    fails: List[str] = []
    fasta_dir, newick_dir = dataset / "fasta", dataset / "newick"
    if not fasta_dir.is_dir() or not newick_dir.is_dir():
        return [f"{dataset} has no fasta/ and newick/ pair"]

    fa_stems = {p.stem for p in fasta_dir.glob("*.fasta")}
    nwk_stems = {p.stem for p in newick_dir.glob("*.nwk")}
    for stem in sorted(fa_stems - nwk_stems):
        fails.append(f"{stem}: fasta without a newick")
    for stem in sorted(nwk_stems - fa_stems):
        fails.append(f"{stem}: newick without a fasta")
    paired = sorted(fa_stems & nwk_stems)
    if not paired:
        return fails + [f"{dataset}: no paired stems"]

    for stem in paired:
        headers, lengths, alphabet = _scan_fasta(fasta_dir / f"{stem}.fasta")
        if len(set(lengths)) > 1:
            fails.append(f"{stem}: not aligned, lengths "
                         f"{min(lengths)}..{max(lengths)}")
        if expect_gapless and not alphabet <= ACGT:
            fails.append(f"{stem}: alphabet {sorted(alphabet)} is not a subset of ACGT "
                         f"-- gaps make JC_similarity_matrix use k>4")
        tree = dendropy.Tree.get(path=str(newick_dir / f"{stem}.nwk"),
                                 schema="newick", preserve_underscores=True)
        leaves = {str(t.taxon.label) for t in tree.leaf_node_iter() if t.taxon}
        if leaves != set(headers):
            missing, extra = leaves - set(headers), set(headers) - leaves
            fails.append(f"{stem}: labels differ -- {len(missing)} only in the newick, "
                         f"{len(extra)} only in the fasta")

    # the one expensive check, on the first alignment only: what the screen will see
    import numpy as np
    import spectraltree
    first = fasta_dir / f"{paired[0]}.fasta"
    cm = dendropy.DnaCharacterMatrix.get(path=str(first), schema="fasta")
    obs, _ = spectraltree.charmatrix2array(cm)
    k = len(np.unique(obs))
    if expect_gapless and k != 4:
        fails.append(f"{paired[0]}: charmatrix2array gives k={k}, not 4 -- "
                     f"JC_similarity_matrix would apply a {k}-state correction")
    print(f"  {len(paired)} paired stems, first alignment {obs.shape}, k={k}")
    return fails

