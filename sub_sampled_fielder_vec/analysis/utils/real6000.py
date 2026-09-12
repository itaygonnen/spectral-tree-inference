"""The m=6000 dataset, as module-level constants.

Thin view over :mod:`analysis.utils.real_cohorts` for the notebook and for callers that
only ever want this one dataset (100 alignments, 6000 taxa x 5000 sites, JC, no indels,
beside 3000 true trees). Anything that should work on *any* real dataset -- the cluster
runner, the interactive launcher -- takes a ``Dataset`` instead.
"""
from __future__ import annotations

from .real_cohorts import Dataset, get_dataset

DATASET_NAME = "6000 taxa"
DATASET: Dataset = Dataset(DATASET_NAME)

SIX_DIR = DATASET.dir
FASTA_DIR = DATASET.fasta_dir
NEWICK_DIR = DATASET.newick_dir

M_TAXA, SEQ_LEN = 6000, 5000
SCREEN_IDS = DATASET.ids() if DATASET.fasta_dir.is_dir() else []


def available_ids(limit=None):
    return DATASET.ids(limit)


def load_all(tree_id: str):
    """``(S, labels, tree, D)`` for one tree of this dataset, or None if files are missing."""
    return DATASET.load_all(tree_id)


__all__ = ["DATASET", "DATASET_NAME", "SIX_DIR", "FASTA_DIR", "NEWICK_DIR",
           "M_TAXA", "SEQ_LEN", "SCREEN_IDS", "available_ids", "load_all", "get_dataset"]
