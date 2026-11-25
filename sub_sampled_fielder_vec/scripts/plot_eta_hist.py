#!/usr/bin/env python3
"""Eta distribution of Sailfish trees, L = D - S against L_sym.

The iso-``r`` grid wants trees at low eta, and the L arm does not supply them: its
cut is lopsided, so the population sits far above every target and acceptance runs
at ~1%. The symmetric Laplacian cuts the same trees differently, and the question
is whether its eta lands where the grid can use it.

Only VALID splits are counted -- a cut that does not separate a real tree edge has
an eta but not a meaning -- and the two arms are gated separately, so the two
panels are not over the same trees.

    python scripts/plot_eta_hist.py
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CACHE = ROOT.parent / "results" / "real_data" / "_cache"
OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid"
ARMS = ("S", "Lsym")
ARM_LABEL = {"S": r"$L = D - S$", "Lsym": r"$L_{sym}$"}


def collect(sizes, seq_len: int = 1000) -> dict:
    """``{arm: {m: etas}}`` for the i.i.d. Sailfish sets.

    The ``_buckets`` sets are excluded on purpose: they are rejection-sampled to
    eta targets, so their histogram is the target and not the distribution.
    """
    out = {a: {} for a in ARMS}
    for m in sizes:
        # Batches of one size are pooled: "_b2" exists only because generation
        # refuses to top up a dataset that already exists, so it is the same
        # i.i.d. draw continued, not a second population.
        parts = sorted(CACHE.glob(f"sailfish_yule_n{m}_l{seq_len}_g1x4"))
        parts += sorted(CACHE.glob(f"sailfish_yule_n{m}_l{seq_len}_g1x4_b*"))
        rows = []
        for part in parts:
            f = part / "screen.npz"
            if f.exists():
                rows.extend(np.load(f, allow_pickle=True)["rows"])
        if not rows:
            continue
        for arm in ARMS:
            e = np.array([q.get(f"eta_{arm}", np.nan) for q in rows], float)
            v = np.array([bool(q.get(f"valid_{arm}", False)) for q in rows])
            keep = v & np.isfinite(e)
            if keep.any():
                out[arm][m] = e[keep]
    return out


def draw(data: dict, sizes):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, LogLocator

    present = [m for m in sizes if any(m in data[a] for a in ARMS)]
    cmap = plt.get_cmap("viridis")
    color = {m: cmap(i / max(1, len(present) - 1)) for i, m in enumerate(present)}

    everything = np.concatenate([v for a in ARMS for v in data[a].values()])
    bins = np.geomspace(1.0, max(everything.max(), 2.0) * 1.05, 28)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8), sharey=True)
    for ax, arm in zip(axes, ARMS):
        for m in present:
            e = data[arm].get(m)
            if e is None:
                continue
            # Filled and translucent, so overlap reads as overlap: where two sizes
            # put their trees in the same place the colours mix rather than one
            # outline hiding behind another.
            ax.hist(e, bins=bins, histtype="stepfilled", alpha=0.45,
                    color=color[m], edgecolor=color[m], lw=1.6,
                    label=f"m={m} (n={len(e)}, med {np.median(e):.1f})")
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1, 2, 5)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_xlabel(r"$\eta$ of the split (valid splits only)", fontsize=11)
        ax.set_title(ARM_LABEL[arm], fontsize=12)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8.5, framealpha=0.95)
    axes[0].set_ylabel("trees", fontsize=11)
    fig.tight_layout()
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[500, 1000, 2000, 4000, 8000])
    ap.add_argument("--seq-len", type=int, default=1000)
    args = ap.parse_args()

    data = collect(args.sizes, args.seq_len)
    if not any(data[a] for a in ARMS):
        print("no i.i.d. Sailfish screens found")
        return 1

    print(f"{'m':>6}  {'arm':<5} {'n':>4} {'median':>7} {'p90':>7} {'max':>7}  "
          f"{'eta<=8':>7}")
    for m in args.sizes:
        for arm in ARMS:
            e = data[arm].get(m)
            if e is None:
                continue
            print(f"{m:>6}  {arm:<5} {len(e):>4} {np.median(e):>7.1f} "
                  f"{np.percentile(e, 90):>7.1f} {e.max():>7.1f}  "
                  f"{100 * (e <= 8).mean():>6.0f}%")

    import matplotlib
    matplotlib.use("Agg")
    fig = draw(data, args.sizes)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "eta_hist_L_vs_Lsym.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
