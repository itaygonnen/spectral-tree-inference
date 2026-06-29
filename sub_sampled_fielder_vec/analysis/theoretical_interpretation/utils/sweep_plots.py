"""Reusable matplotlib for the eta-pool sweep figures.

Extracted from the finalized scratchpad scripts (cmp_fig1 / cmp_fig2 / lsym_fig2)
and the matching cells in the eta_pool_sweep_comparison / eta_pool_sweep_per_n_lsym
notebooks. Behavior is preserved exactly: log axes, free power-law least-squares
fits (label ``n^{b:.2f}``), the red dashed 0.95 threshold line, and the shared ETA
color map.

Data shapes
-----------
``results``        : dict ``method_name -> {(n, eta): ndarray[reps x len(p_values)]}``
``methods``        : spec list ``[(name, color, linestyle, marker), ...]``
``p_values``       : 1-D array-like of p values (log x-axis)
``eta_targets``    : list of eta values (panel / curve identity)
``ns``             : list of taxa counts n

Each public plotting function builds its own figure and calls ``plt.show()``
(matching the notebook cells).
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from analysis.comparison.phase_transition_utils import find_discrete_threshold

# eta -> color (shared across all sweep figures)
ETA_COLORS = {1: "#4b5563", 5: "#1d4ed8", 10: "#ea580c", 15: "#b91c1c"}


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


def plot_pstar_vs_n(results, p_values, eta_targets, ns, methods,
                    agg="mean", threshold=0.95):
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
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(fr"$\eta \approx {eta}$")
        ax.set_xlabel(r"$n$ (taxa)")
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    for row in axes:
        row[0].set_ylabel(r"$p^\star$ (NMI $\geq$ 0.95)")
    fig.suptitle("Figure 1 — $p^\\star$ vs $n$ by operator, per $\\eta$  "
                 "(dashed = free power-law fit $p^\\star\\propto n^{b}$)")
    fig.tight_layout(rect=(0, 0, 1, 0.96)); plt.show()
    return fig


def plot_nmi_grid(results, p_values, eta_targets, ns, methods,
                  agg="mean", threshold=0.95, metric_label="NMI"):
    """Figure 2 -- metric vs p grid (rows: n, cols: eta) + p* summary bottom row.

    Top block: one row per n, metric-vs-p for each operator. Bottom row: each eta
    column collapsed over all n into the scaling panel (markers = p* per operator;
    dashed = free power-law fit p* ~ n^b). Figure-level operator legend is placed
    below the suptitle; bottom-row panels carry per-panel slope-only legends.
    Ports cmp_fig2.py.
    """
    ps = _pstar_all(results, p_values, eta_targets, ns, methods, agg, threshold)
    ng = np.array(sorted(ns), float)
    nr, nc = len(ns), len(eta_targets)
    fig, axes = plt.subplots(nr + 1, nc, figsize=(3.6 * nc, 3.0 * (nr + 1)),
                             squeeze=False)
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
    # --- summary bottom row: p* vs n, one panel per eta ---
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
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_xlabel(r"$n$ (taxa)")
        ax.legend(fontsize=6, loc="lower left", title="fit slope")
        if j == 0:
            ax.set_ylabel("SUMMARY\n" + r"$p^\star$ (NMI $\geq$ 0.95)")
    # --- single operator legend, centered below the title, above the plots ---
    legend_handles = [Line2D([0], [0], color=col, ls=ls, marker=mk, ms=6, lw=1.6,
                             label=name)
                      for name, col, ls, mk in methods]
    fig.suptitle("Figure 2 — NMI vs $p$ (rows: $n$, cols: $\\eta$);  "
                 "bottom row = $p^\\star$ vs $n$ summary per $\\eta$", y=0.995)
    fig.legend(handles=legend_handles, loc="upper center",
               bbox_to_anchor=(0.5, 0.965), ncol=len(methods), fontsize=9,
               frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.95)); plt.show()
    return fig


def plot_curves_per_n(results, p_values, eta_targets, ns, methods,
                      agg="mean", metric_label="NMI", threshold=0.95):
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
    fig.tight_layout(rect=(0, 0, 1, 0.96)); plt.show()
    return fig


def plot_curves_per_eta(results, p_values, eta_targets, ns, methods,
                        agg="mean", metric_label="NMI", threshold=0.95):
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
    fig.tight_layout(rect=(0, 0, 1, 0.96)); plt.show()
    return fig
