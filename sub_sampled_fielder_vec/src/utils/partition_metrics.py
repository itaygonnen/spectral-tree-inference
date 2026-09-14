"""Scoring a bipartition: imbalance, and agreement with a reference.

These four functions were written six and four times respectively across the repo --
``screening``, both eta screens, the real recovery sweep, the real-data benchmark and
``tree_features`` each had their own eta; the NMI/ARI/agreement trio existed in
``real_data_bpart``, ``bpart_sweep_cache``, ``p_sweep_inner`` and the sweep. Identical
arithmetic every time, so a change to one of them (a different tie rule, a different
orientation convention) would have silently split the numbers in two.

Every function here is orientation-invariant: a bipartition and its complement describe
the same split, so ``A|B`` and ``B|A`` must score identically.

Note for callers on the published-figure path (``p_sweep_inner``, ``bpart_sweep_cache``):
those keep their private copies on purpose. They back Figures 2-5 and the gain from
de-duplicating them does not justify touching them.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np


def eta(partition: np.ndarray) -> float:
    """Imbalance of a boolean bipartition: larger clan / smaller clan, in [1, m].

    A one-taxon split of m taxa gives m-1, not infinity: the denominator floors at 1 so a
    degenerate split is a large number rather than a crash in the middle of a sweep.
    """
    part = np.asarray(partition)
    n1 = int(np.sum(part.astype(bool)))
    n2 = int(part.size) - n1
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


def eta_from_sizes(n1: int, n2: int) -> float:
    """``eta`` for a split already known by its two clan sizes."""
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


def agreement(ref: np.ndarray, pred: np.ndarray) -> float:
    """Percent of taxa on the same side, taking the better of the two orientations."""
    a = np.asarray(ref).astype(bool)
    b = np.asarray(pred).astype(bool)
    return 100.0 * float(max((a == b).mean(), (a != b).mean()))


def nmi(ref: np.ndarray, pred: np.ndarray) -> float:
    """Normalized mutual information; label-permutation invariant, so orientation-free."""
    from sklearn.metrics import normalized_mutual_info_score
    return float(normalized_mutual_info_score(np.asarray(ref).astype(int),
                                              np.asarray(pred).astype(int)))


def ari(ref: np.ndarray, pred: np.ndarray) -> float:
    """Adjusted Rand index; 1 = identical split, ~0 = no better than chance."""
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(np.asarray(ref).astype(int),
                                     np.asarray(pred).astype(int)))


def score(ref: np.ndarray, pred: Optional[np.ndarray]) -> Dict[str, float]:
    """``{"nmi", "ari", "agreement"}`` of one partition against a reference.

    ``pred=None`` (a sweep step where every replicate failed) yields NaNs, so a caller
    can keep one value per p without a special case.
    """
    if pred is None:
        return dict(nmi=float("nan"), ari=float("nan"), agreement=float("nan"))
    return dict(nmi=nmi(ref, pred), ari=ari(ref, pred), agreement=agreement(ref, pred))
