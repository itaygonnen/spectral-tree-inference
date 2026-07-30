"""Phase-transition-next-to-scale figure for the eta-pool sweep, one file per eta.

:func:`plot_pstar_pair` emits ONE figure for ONE eta -- left panel the
recovery-vs-p phase transition (curves by n, dashed recovery-threshold line),
right panel the p*(n) scale panel. Nothing in the image names the eta: the
figures are meant to be tiled 2x2 by LaTeX ``subfigure``s whose subcaptions
carry ``$\\eta \\approx k$``, so the paneling lives in the paper, not here.

The scale panel plots the empirical p* markers against a ``C log n / n``
reference: the shape is the theorem's (eq:main_rate at fixed eta), the constant
is measured by :func:`median_ratio_C`. Neither panel carries a legend or an
in-axes annotation: the n colour order, the marker/reference roles and the
fitted C values are stated once, in the LaTeX caption for the tiled figure.

:func:`cbm_sufficient_p` still evaluates eq:main_rate with its literal constant
8(1+eta)^3 (S_out^max)^2 log m / (m*margin^2) from pool-metadata medians, but it
is not drawn: rebuilt per size, it swings by an order of magnitude between
neighbouring n because the median margin does (squared, in the denominator), and
it is nan wherever that margin is non-positive -- at eta>=10 that is most sizes.
The notebook prints it as a diagnostic.

Paper conventions (v9): no axes titles and no suptitle -- the LaTeX caption
carries that text. Type is sized for a subfigure at 0.49\\textwidth.

Reuses the ``sweep_plots`` helpers (``_agg``, ``compute_pstar``, ``ETA_COLORS``,
``_finish``, ``_label_n_axis``) rather than duplicating them.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .sweep_plots import (
    _agg, compute_pstar, ETA_COLORS, _finish, _label_n_axis,
)


# Type sizes are for the printed figure, not the PNG: each file is 7in wide and
# LaTeX renders it in a 0.49\textwidth subfigure (~3.2in), a 0.46x scale.
FS_LABEL = 13
FS_TICK = 10

# Sparse ticks: at subfigure size a full tick set turns into a grey smear, and
# these axes are read for shape, not for values. N_TICKS suits the generated
# regime (n = 500..8000); the synthesized sweep spans m = 90..12288 and passes
# decade ticks instead, via the ``n_ticks`` argument.
NMI_TICKS = [0.2, 0.4, 0.8, 1.0]
N_TICKS = [500, 2000, 8000]


def cbm_sufficient_p(eta: float, s_out_max: float, margin: float, m: int) -> float:
    """CBM sufficient sampling rate, Theorem thm:main-sim / eq:main_rate.

    ``p > 8*(1+eta)^3 * s_out_max^2 * log(m) / (m * margin^2)``. Returns ``nan``
    when ``margin <= 0`` (the certificate is vacuous, \\Cref{cor:infeasible}) or
    ``m <= 1``.
    """
    if margin <= 0.0 or m <= 1:
        return float("nan")
    return 8.0 * (1.0 + eta) ** 3 * s_out_max ** 2 * np.log(m) / (m * margin ** 2)


def _draw_recovery_panel(
    ax: Axes,
    results_for_method: Dict[Tuple[int, int], np.ndarray],
    p_values: np.ndarray,
    eta: int,
    ns: List[int],
    *,
    agg: str = "mean",
    threshold: float = 0.90,
    metric_label: str = "NMI",
) -> None:
    """Left panel: metric vs p, one curve per n, dashed threshold line.

    No legend: at subfigure size a seven-entry key eats a third of the panel.
    The curves run dark-to-light (viridis) in increasing n, which the LaTeX
    caption states once for the whole tiled figure.
    """
    cmap = plt.get_cmap("viridis")
    n_to_color = {n: cmap(i / max(1, len(ns) - 1)) for i, n in enumerate(ns)}
    for n in ns:
        arr = results_for_method.get((n, eta))
        if arr is None or arr.shape[0] == 0:
            continue
        m = _agg(arr, agg)
        ax.plot(p_values, m, "-o", ms=4, lw=1.5, color=n_to_color[n],
                label=fr"$n$={n} ($k$={arr.shape[0]})")
    ax.axhline(threshold, color="red", ls="--", lw=1.2,
               label=f"{metric_label} = {threshold}")
    ax.set_xscale("log"); ax.set_ylim(0, 1.03)
    ax.set_yticks(NMI_TICKS)
    ax.set_xlabel(r"$p$", fontsize=FS_LABEL)
    ax.set_ylabel(metric_label, fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", alpha=0.3)


def median_ratio_C(xs: np.ndarray, ys: np.ndarray) -> Tuple[float, float]:
    """Constant in ``p* = C log n / n`` as the median of the per-point ratios.

    Each measured ``(n_i, p*_i)`` implies its own constant
    ``r_i = p*_i / (log n_i / n_i)``; ``C = median_i r_i``. Robust, weights every
    size equally on log-log axes, and is the estimator the synthesized figure
    already uses (``C = median_m [p*/theory_scale]``), so the two figures'
    constants are comparable.

    Preferred over a linear-scale least-squares fit, which minimises absolute
    error and is therefore set almost entirely by the largest ``p*`` -- i.e. by
    the smallest trees.

    Returns ``(C, spread)`` with ``spread = max_i r_i / min_i r_i``, the factor
    by which a single constant fails to describe the points.
    """
    g = np.log(xs) / xs
    r = ys / g
    return float(np.median(r)), float(r.max() / r.min())


def _draw_scale_panel(
    ax: Axes,
    pstar_rows: List[Tuple[int, float]],
    ng: np.ndarray,
    eta: int,
    *,
    theory_row: Optional[np.ndarray] = None,
    n_ticks: Optional[List[int]] = None,
    size_label: str = r"$n$ (taxa)",
) -> None:
    """Right panel: p*(n) log-log against the ``C log n / n`` reference.

    ``C`` comes from :func:`median_ratio_C`; its value (and the per-point spread
    the same helper returns) is reported by the notebook and quoted in the LaTeX
    caption, not drawn in the axes. ``theory_row`` (eq:main_rate with its literal
    constant) is accepted and not drawn; see the module docstring.
    """
    color = ETA_COLORS.get(eta, "black")
    if pstar_rows:
        xs = np.array([r[0] for r in pstar_rows], float)
        ys = np.array([r[1] for r in pstar_rows], float)
        ax.plot(xs, ys, marker="o", ls="none", ms=8, color=color)
        if len(pstar_rows) >= 2:
            C, _spread = median_ratio_C(xs, ys)
            ax.plot(ng, C * np.log(ng) / ng, ls=":", lw=2.0, color=color, alpha=0.9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ticks = N_TICKS if n_ticks is None else n_ticks
    _label_n_axis(ax, np.array([t for t in ticks if ng.min() <= t <= ng.max()], float))
    ax.set_xlabel(size_label, fontsize=FS_LABEL)
    ax.set_ylabel(r"$p^\star$", fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    ax.grid(True, which="both", alpha=0.3)


def plot_pstar_pair(
    results_for_method: Dict[Tuple[int, int], np.ndarray],
    p_values: np.ndarray,
    eta: int,
    ns: List[int],
    *,
    theory_row: Optional[np.ndarray] = None,
    metric_label: str = "NMI",
    agg: str = "mean",
    threshold: float = 0.90,
    n_ticks: Optional[List[int]] = None,
    size_label: str = r"$n$ (taxa)",
    savepath: Optional[Path] = None,
) -> Figure:
    """One eta, two panels: recovery vs p (left) and p*(n) scale (right).

    Emit one of these per eta and tile them 2x2 with LaTeX ``subfigure``s; the
    eta belongs in the subcaption, so it appears nowhere in the image. Neither
    panel carries a legend -- the n colour order and the marker/reference roles
    are stated once in the figure caption. ``n_ticks`` overrides the default
    x-axis ticks (:data:`N_TICKS`) for sweeps on a different size range, and
    ``size_label`` its label (the synthesized sweep varies m, not n).
    """
    ng = np.array(sorted(ns), float)
    pstar = compute_pstar(results_for_method, p_values, [eta], ns,
                          agg=agg, threshold=threshold)

    # 7.0 x 4.3 renders at ~3.2 x 2.0in in a 0.49\textwidth subfigure: tall
    # enough that the phase transition and the log-log decay both read.
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 4.3))
    _draw_recovery_panel(axes[0], results_for_method, p_values, eta, ns,
                         agg=agg, threshold=threshold, metric_label=metric_label)
    _draw_scale_panel(axes[1], pstar.get(eta, []), ng, eta, theory_row=theory_row,
                      n_ticks=n_ticks, size_label=size_label)
    fig.tight_layout()
    _finish(fig, savepath)
    return fig
