"""Topology-agnostic sign-recovery metric and discrete threshold lookup."""
from __future__ import annotations

from typing import Optional

import numpy as np

from src.utils.metrics import compute_sign_agreement
from analysis.comparison.phase_transition_utils import find_discrete_threshold


def compute_recovery(v_full: np.ndarray, v_hat: np.ndarray) -> float:
    """Sign-agnostic agreement fraction in [0.5, 1] between v_full and v_hat.

    Eigenvectors are only defined up to a global sign, so we report the better
    of the two cluster-label orientations: ``max(s, 1 - s)``. This is the
    natural cluster-recovery floor (random partition labelling gives 0.5),
    independent of the magnitude-weighted dot product, which can disagree with
    the majority sign count on noisy sub-sampled Fiedler vectors.

    The upstream ``_apply_sign_convention`` flips by the first nonzero entry
    of v_hat, which is unreliable at low ``p`` — that's what produced the
    sub-50% scores in figure_3.
    """
    s = compute_sign_agreement(v_full, v_hat) / 100.0
    return max(s, 1.0 - s)


def compute_ari(v_full: np.ndarray, v_hat: np.ndarray) -> float:
    """Adjusted Rand index between the sign-bipartitions of ``v_full`` and ``v_hat``.

    Mirrors the contract of :func:`compute_recovery` — same inputs, scalar
    output — so it slots straight into the ``metrics`` dict consumed by
    :func:`utils.sweep.run_sweep`. ARI is chance-corrected: a random
    bipartition relative to the reference scores 0, a perfect (sign-flipped
    or not) match scores 1. Zero-sign entries (which can appear in noisy
    sub-sampled Fiedlers) get assigned to a third label so they neither
    inflate nor deflate the score.
    """
    from sklearn.metrics import adjusted_rand_score
    a = np.sign(v_full).astype(int)
    b = np.sign(v_hat).astype(int)
    return float(adjusted_rand_score(a, b))


def compute_nmi(v_full: np.ndarray, v_hat: np.ndarray) -> float:
    """Normalized mutual information between the sign-bipartitions of ``v_full`` and ``v_hat``.

    Same contract as :func:`compute_ari` / :func:`compute_recovery` — same
    inputs, scalar output — so it slots straight into the ``metrics`` dict
    consumed by :func:`utils.sweep.run_sweep`. NMI is label-permutation
    invariant (robust to a global sign flip): a perfect match scores 1,
    independent partitions score ~0. Zero-sign entries (which can appear in
    noisy sub-sampled Fiedlers) get assigned to a third label, matching
    :func:`compute_ari`.
    """
    from sklearn.metrics import normalized_mutual_info_score
    a = np.sign(v_full).astype(int)
    b = np.sign(v_hat).astype(int)
    return float(normalized_mutual_info_score(a, b))


def find_threshold_p_star(
    p_vals: np.ndarray,
    agreements: np.ndarray,
    threshold: float = 0.95,
) -> Optional[float]:
    """Smallest p in p_vals (sorted ascending) with agreement >= threshold.

    Pure wrapper around `find_discrete_threshold` — no sigmoid fitting. Both
    `agreements` and `threshold` are in fraction units in [0, 1].
    """
    p_arr = np.asarray(p_vals, dtype=float)
    a_arr = np.asarray(agreements, dtype=float)
    order = np.argsort(p_arr)
    return find_discrete_threshold(p_arr[order], a_arr[order], threshold=threshold)
