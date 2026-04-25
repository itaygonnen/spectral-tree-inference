"""Topology-agnostic two-panel recovery figure and asymptotic constant.

The right-panel constant follows the thesis annotation:
``C = median_n[ p_star * n / ln(n) ]`` — read off the data, not a least-squares
fit. Use ``compute_C_constant`` independently of plotting if needed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.figure import Figure


def compute_C_constant(
    df_thresholds: pd.DataFrame,
    n_col: str = "n",
    p_star_col: str = "p_star",
) -> float:
    """Return ``median_n[p_star * n / ln(n)]`` over rows with finite p_star."""
    df = df_thresholds[[n_col, p_star_col]].dropna()
    df = df[df[p_star_col] > 0]
    if df.empty:
        return float("nan")
    n = df[n_col].to_numpy(dtype=float)
    p_star = df[p_star_col].to_numpy(dtype=float)
    return float(np.median(p_star * n / np.log(n)))


def plot_recovery_figure(
    df_agg: pd.DataFrame,
    df_thresholds: pd.DataFrame,
    *,
    title: str,
    C: Optional[float] = None,
    threshold: float = 0.95,
    savepath: Optional[Path] = None,
) -> Figure:
    """Two-panel figure matching the thesis Fig. 1 layout.

    Left: recovery vs p, one curve per n (viridis). Right: empirical p_star vs
    n with optional ``C * ln(n) / n`` overlay. ``df_agg`` columns: ``n, p,
    agreement`` (agreement in [0, 1]). ``df_thresholds`` columns: ``n,
    p_star``.
    """
    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(title)

    n_values = sorted(df_agg["n"].unique())
    cmap = plt.get_cmap("viridis")
    p_star_lookup = dict(zip(df_thresholds["n"], df_thresholds["p_star"]))

    for i, n in enumerate(n_values):
        sub = df_agg[df_agg["n"] == n].sort_values("p")
        color = cmap(i / max(1, len(n_values) - 1))
        p_star = p_star_lookup.get(n, float("nan"))
        label_p = f"{p_star:.4f}" if p_star is not None and not np.isnan(p_star) else "n/a"
        ax_l.plot(
            sub["p"], 100.0 * sub["agreement"],
            "-o", color=color, ms=4, lw=1,
            label=f"n={n}, p̂*={label_p}",
        )
    ax_l.axhline(100.0 * threshold, color="red", ls="--", lw=1, label=f"{int(threshold*100)}% threshold")
    ax_l.set_xscale("log")
    ax_l.set_xlabel("Sampling probability p")
    ax_l.set_ylabel("Agreement (%)")
    ax_l.set_title("Recovery curves by n")
    ax_l.legend(fontsize=8, loc="lower right")
    ax_l.grid(True, which="both", alpha=0.3)

    df_p = df_thresholds.dropna(subset=["p_star"]).sort_values("n")
    if not df_p.empty:
        ax_r.plot(df_p["n"], df_p["p_star"], "-o", label="Empirical p̂*")
        if C is not None and np.isfinite(C):
            n_grid = np.asarray(df_p["n"], dtype=float)
            ax_r.plot(n_grid, C * np.log(n_grid) / n_grid, "--", label=f"C·ln(n)/n, C={C:.2f}")
    ax_r.set_xscale("log")
    ax_r.set_yscale("log")
    ax_r.set_xlabel("n")
    ax_r.set_ylabel("p̂*")
    ax_r.set_title("Scaling: p̂*(n) vs ln(n)/n")
    ax_r.legend()
    ax_r.grid(True, which="both", alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if savepath is not None:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return fig
