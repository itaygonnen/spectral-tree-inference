"""Utilities for detecting phase transition thresholds in partition agreement curves.

Also the home of the three p*-detection / power-law helpers the paper layer needs.
They previously lived in ``analysis/comparison/phase_transition_utils.py``, which
made the paper-figure code (``analysis/theoretical_interpretation/utils/``) depend
on a superseded analysis package. Nothing under ``src/`` may import from
``analysis/``, so they moved here and the dependency now points the right way.
"""

from __future__ import annotations

from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

ThresholdPoint = Optional[Tuple[float, float]]


def compute_phase_transition_thresholds(
    series: Mapping[str, Mapping[str, Sequence[float]]],
    *,
    target: float = 50.0,
    plateau_tolerance: float = 1.0,
    first_rise_threshold: float = 55.0,
    second_rise_threshold: float = 60.0,
) -> Dict[str, ThresholdPoint]:
    """
    Identify the last ``target``% point before a sustained rise in agreement.

    We scan each curve from left (low p) to right. Whenever a point sits within
    ``plateau_tolerance`` of ``target`` it becomes the current candidate
    threshold. A candidate is confirmed only when the data shows a sustained
    increase: the first point after the plateau must be >= ``first_rise_threshold``
    and the following point must be >= ``second_rise_threshold``. This avoids
    treating brief spikes as true phase transitions.

    Args:
        series: Mapping of label -> {"x": [...], "mean": [...]} sequences.
        target: Expected plateau level (defaults to 50%).
        plateau_tolerance: Allowed deviation from the plateau value.
        first_rise_threshold: Minimum agreement for the first point after the
            plateau.
        second_rise_threshold: Minimum agreement for the second point after the
            plateau; ensures the rise persists.

    Returns:
        Dict mapping each label to a tuple ``(p_value, y_value)`` or ``None`` if
        no transition is detected.
    """

    def on_plateau(value: float) -> bool:
        return value is not None and abs(value - target) <= plateau_tolerance

    thresholds: Dict[str, ThresholdPoint] = {}

    for label, curve in series.items():
        x_vals = list(curve.get("x") or [])
        y_vals = list(curve.get("mean") or [])

        last_plateau: ThresholdPoint = None
        threshold_point: ThresholdPoint = None
        post_plateau_values: list[float] = []

        for y_idx, (p_val, y_val) in enumerate(zip(x_vals, y_vals)):
            if on_plateau(y_val):
                last_plateau = (p_val, y_val)
                post_plateau_values = []
                continue

            if last_plateau is None:
                continue

            post_plateau_values.append(y_val)

            if len(post_plateau_values) == 1:
                if y_val < first_rise_threshold:
                    # First point after plateau failed the guard; wait for a new plateau.
                    last_plateau = None
                    post_plateau_values = []
                continue

            if len(post_plateau_values) == 2:
                first_ok = post_plateau_values[0] >= first_rise_threshold
                second_ok = post_plateau_values[1] >= second_rise_threshold
                if first_ok and second_ok:
                    threshold_point = last_plateau
                    break

                # Guard failed; abandon this plateau candidate until another plateau appears.
                last_plateau = None
                post_plateau_values = []

        if threshold_point is None and last_plateau is not None:
            threshold_point = last_plateau

        thresholds[label] = threshold_point

    return thresholds


def find_discrete_threshold(
    p_vals: np.ndarray, agreements: np.ndarray, threshold: float = 100.0
) -> Optional[float]:
    """First grid ``p`` whose metric reaches ``threshold``; ``None`` if never reached.

    This is the p-hat-star read-off used throughout the paper figures. It is
    deliberately *discrete* -- the smallest value on the sampled p-grid, with no
    interpolation -- so the reported p-hat-star is always a p that was actually
    measured.

    Callers pass the metric on its own scale, so ``threshold`` must match it: 0.90
    for NMI in [0, 1] (the paper's read-off), 90.0 for agreement in percent. The
    paper uses 0.90 rather than 0.95 because at 0.95 the read-off lands on the flat
    top of the transition, where one noisy sample moves p-hat-star a whole grid step.
    """
    for p, agreement in zip(p_vals, agreements):
        if agreement >= threshold:
            return p
    return None


def fit_power_law(n_vals: np.ndarray, p_vals: np.ndarray) -> Tuple[float, float, str]:
    """Least-squares fit of ``p* = A * n^alpha`` in log-log space.

    NaN ``p_vals`` (sizes where p-hat-star was never reached) are dropped; fewer than
    two surviving points returns ``(nan, nan, "Insufficient data")``.

    This is a diagnostic, not one of the paper's metrics. In particular it is NOT
    how the constant ``C`` in ``p* = C log n / n`` is estimated -- that uses the
    median of the per-point ratios, because a log-log slope fit is dominated by
    the largest p*, i.e. by the smallest trees.

    Returns ``(alpha, A, equation_string)``.
    """
    mask = ~np.isnan(p_vals)
    n_clean = n_vals[mask]
    p_clean = p_vals[mask]

    if len(n_clean) < 2:
        return np.nan, np.nan, "Insufficient data"

    log_n = np.log(n_clean)
    log_p = np.log(p_clean)

    coeffs = np.polyfit(log_n, log_p, deg=1)
    alpha = coeffs[0]
    A = np.exp(coeffs[1])

    return alpha, A, f"p* = {A:.2e} × n^{alpha:.3f}"


def evaluate_power_law(n_vals: np.ndarray, alpha: float, A: float) -> np.ndarray:
    """Evaluate the :func:`fit_power_law` result at ``n_vals``: ``A * n^alpha``."""
    return A * (n_vals ** alpha)


__all__ = [
    "compute_phase_transition_thresholds",
    "find_discrete_threshold",
    "fit_power_law",
    "evaluate_power_law",
]

