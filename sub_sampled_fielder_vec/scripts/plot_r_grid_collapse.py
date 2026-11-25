#!/usr/bin/env python3
"""A second reading of the iso-``r`` sweep: the exponent, not the bands.

``plot_r_grid.py`` asks whether each band is flat, one panel per target. That is the
grid's own question, but it answers it four times and never puts a number on the
disagreement. This script asks the same data one question instead:

    how fast does the threshold actually grow with imbalance?

Panel (a) divides the size term out, so every cell lands on one axis and the slope
*is* the exponent. Panel (b) redraws the map with both families of bands -- the ones
the grid was built from, and the ones the measurement implies -- because the finding
is not that a band is tilted but that the bands run in a different direction.

    python scripts/plot_r_grid_collapse.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis.utils.r_grid import EXPONENT, R_TARGETS, R_MAX, eta_for, r_of  # noqa: E402
from plot_r_grid import IN, THRESHOLD, pstar_table, fit_exponent  # noqa: E402

ETA_COLORS = ["#4b5563", "#1d4ed8", "#ea580c", "#b91c1c", "#065f46",
              "#7c3aed", "#b45309", "#0e7490"]


POOL = IN.parent / "pool_nmi.csv"
SAILFISH = IN.parent / "sailfish_nmi.csv"
LSYM = IN.parent / "lsym_nmi.csv"


def load_trees(path: Path, source: str, min_n: int, L: int | None = None):
    """Per-cell medians of one tree population's p_hat_star.

    One tree is not comparable to one synthetic cell: the synthetic cell is a
    closed-form matrix swept 8 times, while a tree is one draw from a prior whose
    p* varies by ~13x at fixed (eta, m). So trees enter as a median per
    (eta bin, m) cell, carrying the min and max of its trees so the scatter the
    median hides stays visible. ``L`` keeps one alignment length per row.
    """
    import pandas as pd

    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path)
    d = d[~d.censored]
    if L is not None and "L" in d:
        d = d[d.L == L]
    g = (d.groupby(["m", "eta_bin"])
           .agg(eta=("eta", "median"), p_star=("p_star", "median"),
                lo=("p_star", "min"), hi=("p_star", "max"), n=("p_star", "size"))
           .reset_index())
    g = g[g.n >= min_n]
    g["source"] = source
    # The trees behind each cell. Held beside the table rather than on it:
    # DataFrame.attrs does not survive the concat that assembles the populations.
    TREES[source] = d[["m", "eta_bin", "eta", "p_star"]].copy()
    g["reached"] = True
    g["censored"] = False
    g["r_target"] = np.nan
    return g


def load_pool():
    """The Kingman eta-pool, as ``load_trees`` reads it."""
    return load_trees(POOL, "kingman_pool", min_n=5)


def load_table(threshold: float = THRESHOLD, sailfish_L: int = 1000):
    """``(table, exponent, spread, n_measured)`` from the sweep on disk.

    The exponent, spread and count are the synthetic grid's; the tree populations
    are appended for ``draw`` to fit on their own.

    Raises ``FileNotFoundError`` when the sweep has not run, and ``ValueError``
    when too few cells are measured to fit -- a partial sweep is the normal state
    while the large cells are still going, so it has to say which it is.
    """
    import pandas as pd

    if not IN.exists():
        raise FileNotFoundError(f"{IN} -- run scripts/sweep_r_grid.py first")
    tab = pstar_table(pd.read_csv(IN), threshold)
    ok = tab[tab.reached & ~tab.censored]
    if len(ok) < 3:
        raise ValueError(f"only {len(ok)} measured cells -- not enough to fit")
    a_fit, spread, n = fit_exponent(tab)
    trees = [load_pool(),
             load_trees(SAILFISH, "sailfish", min_n=2, L=sailfish_L),
             # The same Sailfish trees read through L_sym, which cuts them far
             # more evenly -- a different operator, so a separate population.
             load_trees(LSYM, "lsym", min_n=2, L=sailfish_L)]
    tab = pd.concat([tab, *[t for t in trees if len(t)]], ignore_index=True)
    return tab, a_fit, spread, n


def collapse_spread(ok, exponent: float) -> float:
    """How badly ``r`` at this exponent disagrees with itself.

    ``r`` is proposed as the coordinate that decides recovery, so if it is the right
    one then ``p*`` is a function of ``r`` alone. Fit that function as a straight
    line in log-log and the residual spread -- ``max/min`` of the per-cell ratios,
    the convention ``median_ratio_C`` already uses -- is what is left over.
    """
    rr = r_of(ok.eta.to_numpy(float), ok.m.to_numpy(float), exponent)
    ps = ok.p_star.to_numpy()
    fit = np.polyfit(np.log(rr), np.log(ps), 1)
    res = np.exp(np.log(ps) - np.polyval(fit, np.log(rr)))
    return float(res.max() / res.min())


#: One row per tree model, in this order. Each row is fitted on its own: the
#: models have different constants and possibly different exponents, so one fit
#: across them would average two laws into a third that neither obeys.
#: Tick positions for panel (a), in eta. Round values only; the axis is
#: logarithmic in 1+eta, so these are not evenly spaced.
ETA_TICKS = (1, 2, 3, 5, 8, 12, 16)

#: Per-tree rows behind each population's cell medians, filled by
#: :func:`load_trees`. The synthetic row never appears here: one of its
#: cells is a closed-form matrix, not a sample.
TREES: dict = {}

ROWS = ("synthetic", "kingman_pool", "lsym")
SOURCE_LABEL = {"synthetic": "synthetic flat CBM",
                "kingman_pool": "Kingman pool (cell median)",
                "sailfish": "Yule: IQ-TREE + Sailfish (cell median)",
                # The same Sailfish trees read through the symmetric Laplacian,
                # which cuts them far more evenly: a different operator, so a
                # separate population rather than a variant of the row above.
                "lsym": r"Yule via $L_{sym}$ (cell median)"}


def source_fits(tab) -> dict:
    """``{source: (cells, exponent)}`` for every source with enough cells to fit."""
    out = {}
    for src in ROWS:
        sub = tab[(tab.source == src) & tab.reached & ~tab.censored]
        if len(sub) >= 3:
            out[src] = (sub, fit_exponent(sub)[0])
    return out


def _group(sub):
    """The same-eta key: the integer eta of a synthetic cell, the bin of a tree."""
    if "eta_bin" in sub and sub.eta_bin.notna().all():
        return sub.eta_bin.to_numpy()
    return sub.eta.round().to_numpy()


#: Size classes for the colour. m is discrete: the trees come at fixed sizes and
#: each class holds exactly one of them, while a synthetic cell joins the class
#: its m falls in. Edges sit between the tree sizes.
M_EDGES = (80, 300, 750, 1500, 2500, 3500, 5000, 7000, 13000)
M_LABELS = ("<300", "300-750", "750-1.5k", "1.5k-2.5k", "2.5k-3.5k",
            "3.5k-5k", "5k-7k", ">7k")


def _scaled(sub):
    m = sub.m.to_numpy(float)
    return m / np.log(m)


def draw(tab, a_fit: float = None):
    """Rows = tree model, columns = (a) exponent, (b) grid's r, (c) row's own r.

    Colour is m throughout, on one log scale for the whole figure: eta is already
    an axis in (a) and inside r in (b)/(c), while m was otherwise invisible. In
    (b)/(c) a thin line joins the cells of one eta, so a stack at one r reads as
    either a size effect (colour changes along it) or an eta effect (lines apart).

    (a) anchors both references at the row's balanced end, where there is nothing
    to disagree about -- a reference floated through the cloud is a slope guide,
    not a prediction. ``a_fit`` is accepted for the old call signature; every row
    is fitted from its own cells.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, NullLocator

    fits = source_fits(tab)
    allc = tab[tab.source.isin(list(fits))]
    k = len(M_EDGES) - 1
    cmap = ListedColormap(plt.get_cmap("viridis")(np.linspace(0, 0.95, k)))
    norm = BoundaryNorm(M_EDGES, k)

    fig, axes = plt.subplots(len(fits), 1, figsize=(7.4, 4.4 * len(fits)),
                             squeeze=False, layout="constrained")
    for row, (src, (sub, a_src)) in zip(axes, fits.items()):
        axA = row[0]
        eta, ps, col = sub.eta.to_numpy(float), sub.p_star.to_numpy(), cmap(norm(sub.m))
        has_bar = "lo" in sub and sub.lo.notna().all()

        # ---- (a) size divided out: the slope is the exponent ----------------
        x, y = 1.0 + eta, ps * _scaled(sub)
        if has_bar:
            axA.vlines(x, sub.lo * _scaled(sub), sub.hi * _scaled(sub),
                       color="#9ca3af", lw=1.2, zorder=1)
        # Every tree, faint, under its cell's median. A cell of the synthetic row
        # is one closed-form matrix rather than a sample of trees, so it has no
        # cloud to draw and shows a single marker -- the legend says which is which.
        trees = TREES.get(src)
        if trees is not None and len(trees):
            tx = 1.0 + trees.eta.to_numpy(float)
            ty = trees.p_star.to_numpy() * (trees.m.to_numpy() / np.log(trees.m.to_numpy()))
            axA.scatter(tx, ty, c=cmap(norm(trees.m)), s=14, alpha=0.30,
                        zorder=2, edgecolors="none",
                        label=f"individual trees (n={len(trees)})")
        axA.scatter(x, y, c=col, s=70, zorder=3, edgecolors="white", linewidths=0.8,
                    label="cell median" if trees is not None and len(trees)
                    else "cell (one matrix)")
        x0 = x.min()
        y0 = float(np.exp(np.mean(np.log(y[x == x0]))))
        xs = np.linspace(x.min() * 0.95, x.max() * 1.05, 60)
        axA.plot(xs, y0 * (xs / x0) ** a_src, "-", color="#111827", lw=2,
                 label=rf"measured: $(1+\eta)^{{{a_src:.2f}}}$")
        axA.plot(xs, y0 * (xs / x0) ** EXPONENT, "--", color="#b91c1c", lw=2,
                 label=rf"theory: $(1+\eta)^{{{EXPONENT:g}}}$")
        axA.set_xscale("log"); axA.set_yscale("log")
        axA.set_ylabel(f"{SOURCE_LABEL[src]}\n" + r"$\hat p^\star \cdot m / \log m$",
                       fontsize=10.5)
        axA.legend(loc="upper left", fontsize=9, framealpha=0.95)


    for ax in axes.ravel():
        ax.grid(True, which="both", alpha=0.3)
    # (a) shares x as well as y, so the rows compare at a glance; the limits span
    # every row, since the pool reaches eta~16, well past the synthetic grid.
    x_all = 1.0 + allc.eta.to_numpy(float)
    y_all = allc.p_star.to_numpy() * _scaled(allc)
    for ax in axes[:, 0]:
        ax.set_xlim(x_all.min() * 0.92, x_all.max() * 1.08)
        ax.set_ylim(y_all.min() * 0.4, y_all.max() * 2.4)
        # Ticked and labelled in eta, though the coordinate is 1+eta: the power
        # law needs the shift, but an axis reading "17" for a tree at eta=16 is
        # a trap for whoever reads the figure.
        ax.xaxis.set_major_locator(FixedLocator([1 + e for e in ETA_TICKS]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v - 1:g}"))
    axes[0, 0].set_title("(a)  size divided out: the slope is the exponent",
                         fontsize=10.5)
    axes[-1, 0].set_xlabel(r"clan imbalance $\eta$", fontsize=10.5)
    cb = fig.colorbar(ScalarMappable(norm, cmap), ax=axes, label="number of taxa $m$",
                      shrink=0.6, spacing="uniform")
    cb.set_ticks([(a * b) ** 0.5 for a, b in zip(M_EDGES, M_EDGES[1:])],
                 labels=list(M_LABELS))
    cb.ax.minorticks_off()
    return fig


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=IN.parent)
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    ap.add_argument("--sailfish-L", type=int, default=1000,
                    help="alignment length of the Yule row (one L per row)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    try:
        tab, a_fit, spread, n = load_table(args.threshold, args.sailfish_L)
    except (FileNotFoundError, ValueError) as exc:
        print(exc)
        return 1

    print(f"{n} measured cells   fitted exponent {a_fit:.2f}   "
          f"constants span {spread:.2f}x   grid built at {EXPONENT:g}")
    for src, (sub, a_src) in source_fits(tab).items():
        print(f"  {src:13s} {len(sub):3d} cells  exponent {a_src:.2f}  "
              f"left over {collapse_spread(sub, EXPONENT):.1f}x at {EXPONENT:g}, "
              f"{collapse_spread(sub, a_src):.1f}x at own")

    fig = draw(tab, a_fit)
    out = args.out / "exponent_collapse.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"wrote {out}")
    return 0


def _m_for_eta1(r: float) -> float:
    """The size at which eta=1 sits on target ``r`` under the grid's exponent."""
    from analysis.utils.r_grid import m_for
    return m_for(r, 1.0, EXPONENT)


if __name__ == "__main__":
    raise SystemExit(main())
