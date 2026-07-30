"""The two figures of the operator-comparison benchmark.

Ports the real-data notebook's recovery-curve (cell 22) and scale-plot (cell 24)
so they can be produced headlessly. ``Agg`` is selected before ``pyplot`` is
imported anywhere, so this module is safe to import over SSH with no display.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# (key, legend label, color, marker) — the canonical operator triple.
METHODS = [
    ("G", "Griffing-on-$D$", "#065f46", "s"),
    ("L", "$L$ + k-means", "#4b5563", "o"),
    ("Lsym", "$L_{\\rm sym}$ + k-means", "#b45309", "^"),
]

P_STAR_THRESHOLD = 0.95


def p_star_per_tree(nmi_arr: np.ndarray, p_values: np.ndarray,
                    thresh: float = P_STAR_THRESHOLD) -> np.ndarray:
    """Smallest p whose NMI reaches ``thresh``, per tree; ``inf`` if never."""
    out = []
    for row in np.asarray(nmi_arr, float):
        ok = np.where(row >= thresh)[0]
        out.append(float(p_values[ok[0]]) if ok.size else np.inf)
    return np.array(out)


def plot_recovery_curve(p_values: Sequence[float], curves: Dict[str, np.ndarray],
                        out_path: Path, n_taxa: Optional[int] = None) -> Path:
    """Median NMI vs sub-sampling fraction p, with a ±1 std band per operator."""
    P = np.asarray(p_values, float)
    n_trees = next(iter(curves.values())).shape[0]

    fig, ax = plt.subplots(figsize=(7, 7))
    for key, name, col, _ in METHODS:
        data = curves[key]
        med = np.nanmedian(data, 0)
        sd = np.nanstd(data, 0)
        ax.plot(P, med, color=col, lw=2, marker="o", ms=3, label=name)
        ax.fill_between(P, med - sd, med + sd, color=col, alpha=0.15)
    ax.set_xscale("log")
    ax.set_xlabel("sub-sampling fraction $p$ (log)", fontsize=12)
    ax.set_ylabel("NMI vs full-matrix reference", fontsize=12)
    ax.set_xlim(P.min() * 0.8, 1.2)
    ax.set_ylim(-0.05, 1.05)
    _h, _l = ax.get_legend_handles_labels()
    n_str = f", n={n_taxa}" if n_taxa else ""
    fig.suptitle(
        "Sub-sampling recovery: Griffing-$D$ vs $L$+k-means vs "
        "$L_{\\rm sym}$+k-means\n"
        f"median NMI over {n_trees} trees (band = $\\pm$1 std){n_str}",
        fontsize=12, y=0.99)
    fig.legend(_h, _l, loc="upper center", bbox_to_anchor=(0.5, 0.90),
               ncol=len(_l), fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return Path(out_path)


def plot_scale(p_values: Sequence[float], curves: Dict[str, np.ndarray],
               rT: np.ndarray, out_path: Path) -> Path:
    """p* (smallest p with NMI ≥ 0.95) vs r(T) = tree diameter, log-log."""
    P = np.asarray(p_values, float)
    rT = np.asarray(rT, float)
    n_trees = len(rT)
    pstars = {key: p_star_per_tree(curves[key], P) for key, _, _, _ in METHODS}

    fig, ax = plt.subplots(figsize=(8, 6))
    for key, name, col, mk in METHODS:
        ps = pstars[key]
        m = np.isfinite(ps) & (ps > 0)
        ax.scatter(rT[m], ps[m], color=col, marker=mk, alpha=0.55, s=35, label=name)

    if any((np.isfinite(pstars[k]) & (pstars[k] > 0)).sum() >= 6
           for k, _, _, _ in METHODS):
        edges = np.quantile(rT, [0, 1 / 3, 2 / 3, 1.0])
        centers = 0.5 * (edges[:-1] + edges[1:])
        for key, _, col, mk in METHODS:
            ps = pstars[key]
            m = np.isfinite(ps) & (ps > 0)
            meds = [np.median(ps[m & (rT >= lo) & (rT <= hi)])
                    if (m & (rT >= lo) & (rT <= hi)).any() else np.nan
                    for lo, hi in zip(edges[:-1], edges[1:])]
            ax.plot(centers, meds, color=col, lw=2.5, marker=mk, ms=10,
                    mec="white", mew=1.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$r(T)$  (tree diameter under JC distance)", fontsize=12)
    ax.set_ylabel(r"$p^\star$  per tree (smallest $p$ with NMI$\geq 0.95$)",
                  fontsize=12)
    ax.set_title(f"$p^\\star$ vs $r(T)$ on the cohort (n={n_trees})\n"
                 "(thick markers: tertile medians; thin: per-tree scatter)",
                 fontsize=12)
    ax.legend(fontsize=11, loc="best")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return Path(out_path)


def summarize(p_values: Sequence[float], curves: Dict[str, np.ndarray],
              rT: np.ndarray) -> Dict[str, object]:
    """Numbers worth printing / persisting alongside the two figures."""
    P = np.asarray(p_values, float)
    rT = np.asarray(rT, float)
    out: Dict[str, object] = {
        "n_trees": int(len(rT)),
        "rT_min": float(np.min(rT)), "rT_max": float(np.max(rT)),
        "rT_median": float(np.median(rT)),
        "methods": {},
    }
    for key, name, _, _ in METHODS:
        ps = p_star_per_tree(curves[key], P)
        m = np.isfinite(ps) & (ps > 0)
        out["methods"][key] = {
            "label": name,
            "p_star_finite": int(m.sum()),
            "p_star_median": float(np.nanmedian(ps[m])) if m.any() else None,
            "median_nmi_at_min_p": float(np.nanmedian(curves[key][:, 0])),
        }
    return out
