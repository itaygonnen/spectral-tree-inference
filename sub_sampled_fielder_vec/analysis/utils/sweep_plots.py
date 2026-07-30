"""Reusable matplotlib for the eta-pool sweep figures.

Serves the APPENDIX figures of ``paper/fig03_pstar_gen_kingman.ipynb``:
``plot_nmi_grid`` -> Fig 8 (fig:recovery_grid) and ``plot_pstar_vs_n`` -> Fig 9
(fig:operator_sensitivity). The paper's MAIN per-eta panels (Fig 3) come from
:mod:`sweep_plots_two_panel` instead -- see the note on constants below.

Extracted from the finalized scratchpad scripts (cmp_fig1 / cmp_fig2 / lsym_fig2)
and the matching cells in the eta_pool_sweep_comparison / eta_pool_sweep_per_n_lsym
notebooks. Behavior is preserved exactly: log axes, free power-law least-squares
fits (label ``n^{b:.2f}``), the red dashed 0.95 threshold line, and the shared ETA
color map.

Two estimators of ``C`` coexist in this package -- deliberately, but they are NOT
interchangeable:

- :func:`_theory_ref` here fits ``C`` in ``p* = C log n / n`` by least squares.
- :func:`sweep_plots_two_panel.median_ratio_C` uses the median of the per-point
  ratios, which is the estimator the paper's constants are quoted from: an LS fit
  on a linear scale is set almost entirely by the largest ``p*``, i.e. by the
  smallest trees.

So a ``C`` printed by this module is NOT comparable with one from the two-panel
module, and the free exponent ``b`` in the ``n^{b:.2f}`` labels is scaffolding for
the operator comparison, not one of the paper's metrics. Neither curve has been
changed here -- a plotted curve is a paper claim -- but the discrepancy is on record
for the author in ``docs/overleafs/v9/open-items/19-cleanup.md``.

Data shapes
-----------
``results``        : dict ``method_name -> {(n, eta): ndarray[reps x len(p_values)]}``
``methods``        : spec list ``[(name, color, linestyle, marker), ...]``
``p_values``       : 1-D array-like of p values (log x-axis)
``eta_targets``    : list of eta values (panel / curve identity)
``ns``             : list of taxa counts n

Each public plotting function builds its own figure and ends in :func:`_finish`,
which saves when ``savepath`` is given and always calls ``plt.show()`` (so the
``savepath=None`` behaviour is exactly what the notebook cells produced before).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter

from src.utils.threshold_utils import find_discrete_threshold

# eta -> color (shared across all sweep figures)
ETA_COLORS = {1: "#4b5563", 5: "#1d4ed8", 10: "#ea580c", 15: "#b91c1c"}


def _finish(fig: Figure, savepath: Optional[Path] = None) -> None:
    """Shared figure tail: save if asked, then always show.

    ``savepath=None`` reduces to the historical ``plt.show()`` tail, so a caller
    that passes nothing behaves exactly as the notebook cells did before.
    """
    if savepath is not None:
        savepath = Path(savepath)
        savepath.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, dpi=200, bbox_inches="tight")
        print(f"saved {savepath}")
    plt.show()


def _label_n_axis(ax, ng):
    """Force integer n labels at the actual sampled n on a log x-axis.

    The default LogLocator only labels decade ticks (e.g. just ``10^3`` for a
    500-8000 span), so set explicit ticks at ng with a plain integer formatter.
    """
    ax.set_xticks(ng)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.ticklabel_format(axis="x", style="plain")
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(45)
        lbl.set_fontsize(7)


def _agg(arr: np.ndarray, agg: str) -> np.ndarray:
    """Collapse the reps axis (axis 0) by mean or median."""
    if agg == "median":
        return np.median(arr, axis=0)
    return arr.mean(0)


def compute_pstar(results_for_method, p_values, eta_targets, ns,
                  agg="mean", threshold=0.95):
    """p* per (eta, n) for a single operator's results.

    Parameters
    ----------
    results_for_method : dict ``{(n, eta): ndarray[reps x len(p_values)]}``
    agg : {"mean", "median"} -- reps aggregation before find_discrete_threshold.
    threshold : metric threshold defining p* (default 0.95).

    Returns
    -------
    dict ``eta -> list[(n, p*)]`` (grid-quantized p* via find_discrete_threshold;
    only finite thresholds are kept).
    """
    parr = np.asarray(p_values, dtype=float)
    out = {}
    for eta in eta_targets:
        rows = []
        for n in ns:
            arr = results_for_method.get((n, eta))
            if arr is None or arr.shape[0] == 0:
                continue
            t = find_discrete_threshold(parr, _agg(arr, agg), threshold=threshold)
            if t is not None and np.isfinite(t):
                rows.append((n, float(t)))
        out[eta] = rows
    return out


def _pstar_all(results, p_values, eta_targets, ns, methods, agg, threshold):
    """p* dict for every operator: name -> {eta -> [(n, p*)]}."""
    return {name: compute_pstar(results[name], p_values, eta_targets, ns,
                                agg=agg, threshold=threshold)
            for name, *_ in methods}


def _theory_ref(xs, ys, ng):
    """1-parameter LS fit of the theoretical rate ``p* = C * log(n)/n``.

    Returns (C, ng, C*log(ng)/ng) for overlaying as a reference curve. The single
    free constant C is chosen to minimize squared error against the measured p*
    points, so the curve tests the *shape* log n / n rather than a free exponent.
    """
    g = np.log(xs) / xs
    C = float(np.sum(ys * g) / np.sum(g * g))
    return C, ng, C * np.log(ng) / ng


def plot_pstar_vs_n(results, p_values, eta_targets, ns, methods,
                    agg="mean", threshold=0.95, savepath: Optional[Path] = None):
    """Figure 1 -- p* vs n, one panel per eta, one line per operator.

    Dashed = free power-law least-squares fit ``log p* = a + b*log n`` (=>
    ``p* = A * n^b``); slope b is the measured scaling exponent. Works for 1+
    operators. Ports cmp_fig1.py.
    """
    ps = _pstar_all(results, p_values, eta_targets, ns, methods, agg, threshold)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
    axf = axes.flatten()
    ng = np.array(sorted(ns), float)
    for ax, eta in zip(axf, eta_targets):
        for name, col, ls, mk in methods:
            rows = ps[name].get(eta, [])
            if not rows:
                continue
            xs = np.array([r[0] for r in rows], float)
            ys = np.array([r[1] for r in rows], float)
            ax.plot(xs, ys, marker=mk, ls="none", ms=6, color=col, label=name)
            if len(rows) >= 2:  # free power-law best-fit (log-log least squares)
                b, a = np.polyfit(np.log(xs), np.log(ys), 1)
                ax.plot(ng, np.exp(a) * ng ** b, ls="--", lw=1.3, color=col,
                        alpha=0.8, label=fr"  fit $n^{{{b:.2f}}}$")
                C, xr, yr = _theory_ref(xs, ys, ng)  # theoretical C log n / n
                ax.plot(xr, yr, ls=":", lw=1.3, color=col, alpha=0.8,
                        label=fr"  ref ${C:.2f}\log n/n$")
        ax.set_xscale("log"); ax.set_yscale("log")
        _label_n_axis(ax, ng)
        ax.set_title(fr"$\eta \approx {eta}$")
        ax.set_xlabel(r"$n$ (taxa)")
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    for row in axes:
        row[0].set_ylabel(r"$p^\star$ (NMI $\geq$ 0.95)")
    fig.suptitle("Figure 1 — $p^\\star$ vs $n$ by operator, per $\\eta$  "
                 "(dashed = free power-law fit $p^\\star\\propto n^{b}$)")
    fig.tight_layout(rect=(0, 0, 1, 0.96)); _finish(fig, savepath)
    return fig


def plot_nmi_grid(results, p_values, eta_targets, ns, methods,
                  agg="mean", threshold=0.95, metric_label="NMI", summary=True,
                  savepath: Optional[Path] = None):
    """Figure 2 -- metric vs p grid (rows: n, cols: eta), optional p* summary row.

    Top block: one row per n, metric-vs-p for each operator. When ``summary`` is
    True, a bottom row adds the p* vs n scaling panel per eta (markers = p* per
    operator; dashed = free power-law fit p* ~ n^b). Set ``summary=False`` when a
    standalone p* vs n figure (``plot_pstar_vs_n``) is already shown -- the summary
    row is identical to it. Figure-level operator legend sits below the suptitle.
    Ports cmp_fig2.py.
    """
    ps = _pstar_all(results, p_values, eta_targets, ns, methods, agg, threshold) if summary else {}
    ng = np.array(sorted(ns), float)
    nr, nc = len(ns), len(eta_targets)
    SUMMARY_H = 2.6  # summary row height relative to a metric-vs-p row
    hr = [1.0] * nr + ([SUMMARY_H] if summary else [])
    fig, axes = plt.subplots(nr + (1 if summary else 0), nc,
                             figsize=(3.6 * nc, 3.0 * (nr + (SUMMARY_H if summary else 0))),
                             gridspec_kw={"height_ratios": hr}, squeeze=False)
    # --- top block: metric vs p, one row per n ---
    for i, n in enumerate(ns):
        for j, eta in enumerate(eta_targets):
            ax = axes[i][j]
            for name, col, ls, mk in methods:
                arr = results[name].get((n, eta))
                if arr is None or arr.shape[0] == 0:
                    continue
                ax.plot(p_values, _agg(arr, agg), ls=ls, marker=mk, ms=3,
                        lw=1.3, color=col)
            ax.axhline(threshold, color="red", ls=":", lw=0.8)
            ax.set_xscale("log"); ax.set_ylim(0, 1.03)
            ax.grid(True, which="both", alpha=0.25)
            if i == 0:
                ax.set_title(fr"$\eta \approx {eta}$")
            if j == 0:
                ax.set_ylabel(f"n={n}\n{metric_label}")
            if i == nr - 1:
                ax.set_xlabel(r"$p$")
    # --- summary bottom row: p* vs n, one panel per eta (optional) ---
    if summary:
        for j, eta in enumerate(eta_targets):
            ax = axes[nr][j]
            for name, col, ls, mk in methods:
                rows = ps[name].get(eta, [])
                if not rows:
                    continue
                xs = np.array([r[0] for r in rows], float)
                ys = np.array([r[1] for r in rows], float)
                ax.plot(xs, ys, marker=mk, ls="none", ms=6, color=col)
                if len(rows) >= 2:
                    b, a = np.polyfit(np.log(xs), np.log(ys), 1)
                    ax.plot(ng, np.exp(a) * ng ** b, ls="--", lw=1.3, color=col,
                            alpha=0.8, label=fr"$n^{{{b:.2f}}}$")
                    C, xr, yr = _theory_ref(xs, ys, ng)  # theoretical C log n / n
                    ax.plot(xr, yr, ls=":", lw=1.3, color=col, alpha=0.8,
                            label=fr"${C:.2f}\log n/n$")
            ax.set_xscale("log"); ax.set_yscale("log")
            _label_n_axis(ax, ng)
            ax.grid(True, which="both", alpha=0.3)
            ax.set_xlabel(r"$n$ (taxa)")
            ax.legend(fontsize=6, loc="lower left", title="fit slope")
            if j == 0:
                ax.set_ylabel("SUMMARY\n" + r"$p^\star$ (NMI $\geq$ 0.95)")
    # --- single operator legend, centered below the title, above the plots ---
    legend_handles = [Line2D([0], [0], color=col, ls=ls, marker=mk, ms=6, lw=1.6,
                             label=name)
                      for name, col, ls, mk in methods]
    _sub = ("; bottom row = $p^\\star$ vs $n$ summary per $\\eta$" if summary else "")
    fig.suptitle(f"Figure 2 — {metric_label} vs $p$ (rows: $n$, cols: $\\eta$)" + _sub,
                 y=0.995)
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.965), ncol=len(methods), fontsize=9,
               frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); _finish(fig, savepath)
    return fig


def plot_curves_per_n(results, p_values, eta_targets, ns, methods,
                      agg="mean", metric_label="NMI", threshold=0.95,
                      savepath: Optional[Path] = None):
    """metric vs p, one panel per n, curves by eta (per_n / lsym Figure-3 style).

    Curve identity is the eta color (ETA_COLORS); operators are distinguished by
    their (linestyle, marker). Ports the lsym Figure-3 cell.
    """
    cmap = plt.get_cmap("viridis")
    eta_to_color = {e: ETA_COLORS.get(e, cmap(i / max(1, len(eta_targets) - 1)))
                    for i, e in enumerate(eta_targets)}
    ncols = 2
    nrows = (len(ns) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             sharey=True, squeeze=False)
    axf = axes.flatten()
    fig.suptitle(fr"{metric_label} vs $p$ — panel per $n$, curves by $\eta$")
    for ax, n in zip(axf, ns):
        for eta in eta_targets:
            for name, col, ls, mk in methods:
                arr = results[name].get((n, eta))
                if arr is None or arr.shape[0] == 0:
                    continue
                m = _agg(arr, agg)
                lbl = (fr"$\eta$={eta} ($k$={arr.shape[0]})"
                       if len(methods) == 1 else fr"$\eta$={eta} [{name}]")
                ax.plot(p_values, m, ls=ls, marker=mk, ms=4, lw=1.5,
                        color=eta_to_color[eta], label=lbl)
        ax.axhline(threshold, color="red", ls="--", lw=1,
                   label=f"{metric_label} = {threshold}")
        ax.set_xscale("log"); ax.set_title(fr"$n$ = {n}"); ax.set_xlabel(r"$p$")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    for ax in axf[len(ns):]:
        ax.set_visible(False)
    for row in axes:
        row[0].set_ylabel(metric_label)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); _finish(fig, savepath)
    return fig


def plot_curves_per_eta(results, p_values, eta_targets, ns, methods,
                        agg="mean", metric_label="NMI", threshold=0.95,
                        savepath: Optional[Path] = None):
    """metric vs p, one panel per eta, curves by n (lsym Figure-1 style).

    Curve identity is the n color (viridis ramp); operators are distinguished by
    their (linestyle, marker). Ports the lsym Figure-1 cell.
    """
    cmap = plt.get_cmap("viridis")
    n_to_color = {n: cmap(i / max(1, len(ns) - 1)) for i, n in enumerate(ns)}
    ncols = 2
    nrows = (len(eta_targets) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             sharey=True, squeeze=False)
    axf = axes.flatten()
    fig.suptitle(fr"{metric_label} vs $p$ — panel per $\eta$, curves by $n$")
    for ax, eta in zip(axf, eta_targets):
        for n in ns:
            for name, col, ls, mk in methods:
                arr = results[name].get((n, eta))
                if arr is None or arr.shape[0] == 0:
                    continue
                m = _agg(arr, agg)
                lbl = (fr"$n$={n} ($k$={arr.shape[0]})"
                       if len(methods) == 1 else fr"$n$={n} [{name}]")
                ax.plot(p_values, m, ls=ls, marker=mk, ms=4, lw=1.5,
                        color=n_to_color[n], label=lbl)
        ax.axhline(threshold, color="red", ls="--", lw=1,
                   label=f"{metric_label} = {threshold}")
        ax.set_xscale("log"); ax.set_title(fr"$\eta \approx {eta}$")
        ax.set_xlabel(r"$p$")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, loc="lower right")
    for ax in axf[len(eta_targets):]:
        ax.set_visible(False)
    for row in axes:
        row[0].set_ylabel(metric_label)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); _finish(fig, savepath)
    return fig
