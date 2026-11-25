#!/usr/bin/env python3
"""Cell occupancy of the iso-``r`` grid, L against L_sym, on the same trees.

The grid needs trees at the ``r`` each cell asks for. Which operator supplies them
is not a matter of taste: ``r = (1+eta)^3 log m/m`` is cubic in the imbalance, so
L's lopsided cuts throw their trees far up the r axis while L_sym's land in the
band the grid is built on. Two panels of the same trees, counted the same way,
is the whole argument.

    python scripts/plot_grid_occupancy.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis.utils.r_grid import LSYM_TARGETS, bin_by_r  # noqa: E402

OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid"
ARM_LABEL = {"S": r"$L = D - S$", "Lsym": r"$L_{sym}$"}


def occupancy(sizes, arm: str, targets=LSYM_TARGETS, tol: float = 0.18,
              seq_len: int = 1000) -> np.ndarray:
    """``(len(targets)+1, len(sizes))`` counts; the last row is "matched no cell"."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from plot_eta_hist import collect

    data = collect(sizes, seq_len)[arm]
    grid = np.zeros((len(targets) + 1, len(sizes)), int)
    for j, m in enumerate(sizes):
        e = data.get(m)
        if e is None or not len(e):
            continue
        lab = bin_by_r(e, np.full(len(e), float(m)), targets=targets, tol=tol)
        for i, t in enumerate(targets):
            grid[i, j] = int((lab == t).sum())
        grid[-1, j] = int(np.isnan(lab).sum())
    return grid


def draw(sizes, targets=LSYM_TARGETS, tol: float = 0.18, seq_len: int = 1000):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    arms = ("S", "Lsym")
    grids = {a: occupancy(sizes, a, targets, tol, seq_len) for a in arms}
    top = max(g.max() for g in grids.values())

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), sharey=True)
    rows = [f"{t:g}" for t in targets] + ["none"]
    for ax, arm in zip(axes, arms):
        g = grids[arm]
        # Zero is masked rather than coloured: an empty cell is the finding, and a
        # dark square at the bottom of a log scale reads as "nearly none" instead.
        ax.imshow(np.where(g > 0, g, np.nan), cmap="viridis",
                  norm=LogNorm(vmin=1, vmax=max(top, 2)), aspect="auto")
        for i in range(g.shape[0]):
            for j in range(g.shape[1]):
                v = g[i, j]
                ax.text(j, i, "-" if v == 0 else str(v), ha="center", va="center",
                        fontsize=9, color="#111827" if v > top / 3 else "#e5e7eb")
        ax.set_xticks(range(len(sizes)), [str(m) for m in sizes])
        ax.set_yticks(range(len(rows)), rows)
        ax.axhline(len(targets) - 0.5, color="#111827", lw=1.5)
        ax.set_xlabel(r"$m$", fontsize=11)
        inside = int(g[:-1].sum())
        ax.set_title(f"{ARM_LABEL[arm]}  —  {inside} trees in cells, "
                     f"{int(g[-1].sum())} in none", fontsize=11)
    axes[0].set_ylabel(rf"iso-$r$ cell   (within {tol:.0%})", fontsize=11)
    fig.tight_layout()
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[500, 1000, 2000, 3000, 4000, 5000, 8000])
    ap.add_argument("--tol", type=float, default=0.18)
    args = ap.parse_args()

    for arm in ("S", "Lsym"):
        g = occupancy(args.sizes, arm, tol=args.tol)
        print(f"{ARM_LABEL[arm]:>14}: {int(g[:-1].sum()):4d} in cells, "
              f"{int(g[-1].sum()):4d} in none, "
              f"{int((g[:-1] == 0).sum())} of {g[:-1].size} cells empty")

    import matplotlib
    matplotlib.use("Agg")
    fig = draw(args.sizes, tol=args.tol)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "grid_occupancy_L_vs_Lsym.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
