"""The m=6000 cohort, as module-level constants.

Thin view over :mod:`analysis.utils.real_cohorts` for the notebook and for callers that
only ever want this one dataset (100 alignments, 6000 taxa x 5000 sites, JC, no indels,
beside 3000 true trees). Anything that should work on *any* real cohort -- the cluster
runner, the interactive launcher -- takes a ``Cohort`` instead.
"""
from __future__ import annotations

from .real_cohorts import Cohort, get_cohort

COHORT_NAME = "6000 taxa"
COHORT: Cohort = Cohort(COHORT_NAME)

SIX_DIR = COHORT.dir
FASTA_DIR = COHORT.fasta_dir
NEWICK_DIR = COHORT.newick_dir

M_TAXA, SEQ_LEN = 6000, 5000
SCREEN_IDS = COHORT.ids() if COHORT.fasta_dir.is_dir() else []


def available_ids(limit=None):
    return COHORT.ids(limit)


def load_all(tree_id: str):
    """``(S, labels, tree, D)`` for one tree of this cohort, or None if files are missing."""
    return COHORT.load_all(tree_id)


__all__ = ["COHORT", "COHORT_NAME", "SIX_DIR", "FASTA_DIR", "NEWICK_DIR",
           "M_TAXA", "SEQ_LEN", "SCREEN_IDS", "available_ids", "load_all", "get_cohort"]
