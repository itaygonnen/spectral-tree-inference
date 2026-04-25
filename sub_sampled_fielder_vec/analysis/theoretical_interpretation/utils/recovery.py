"""Topology-agnostic sign-recovery metric and discrete threshold lookup."""
from __future__ import annotations

from typing import Optional

import numpy as np

from src.utils.metrics import compute_sign_agreement
from analysis.comparison.phase_transition_utils import find_discrete_threshold


def compute_recovery(v_full: np.ndarray, v_hat: np.ndarray) -> float:
    """Sign-agreement fraction in [0, 1] between v_full and v_hat.

    The Fiedler sign convention is enforced upstream by FiedlerVectorComputer,
    so callers can assume v_full and v_hat are already aligned.
    """
    return compute_sign_agreement(v_full, v_hat) / 100.0


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
