#!/usr/bin/env python3
"""Read off ``p_hat_star`` per iso-``r`` cell and plot it against size.

The grid's claim is that cells sharing an ``r`` share a recovery threshold. So the
figure is one panel per ``r`` target, ``p_hat_star`` against ``m``, and the claim is
read as *flat or not flat*: flat confirms the exponent the grid was built with,
sloped measures the correction.

A band is only evidence if its cells were measured, so a cell whose threshold is
pinned at either end of the ``p`` grid is drawn open and left out of the fit --
there the sweep reports a bound, not a threshold.

    python scripts/plot_r_grid.py
    python scripts/plot_r_grid.py --out /tmp/scratch     # keep the paper dir clean
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis.utils.r_grid import EXPONENT, r_of   # noqa: E402
from src.utils.threshold_utils import find_discrete_threshold  # noqa: E402

IN = ROOT / "analysis" / "notebooks_cache" / "r_grid" / "sweep_nmi.csv"
THRESHOLD = 0.90          # the paper's level; 0.95 lands on the flat top
ETA_COLORS = ["#4b5563", "#1d4ed8", "#ea580c", "#b91c1c", "#065f46",
              "#7c3aed", "#b45309", "#0e7490"]


def pstar_table(df: pd.DataFrame, threshold: float = THRESHOLD) -> pd.DataFrame:
    """One row per cell: its ``p_hat_star`` and whether it was measured.

    ``find_discrete_threshold`` is the estimator the paper's figures use -- the first
    grid ``p`` at or above the level, with no interpolation. Reusing it is what makes
    a number here comparable to one there.
    """
    rows = []
    keys = ["r_target", "eta", "m"] + (["source"] if "source" in df else [])
    agg = (df.groupby(keys + ["p"], as_index=False)["nmi"].mean().sort_values("p"))
    for key, sub in agg.groupby(keys):
        r_t, eta, m = key[0], key[1], key[2]
        p = sub["p"].to_numpy()
        ps = find_discrete_threshold(p, sub["nmi"].to_numpy(), threshold)
        rows.append({
            # Which simulator the cell came from. The flat CBM builds eta exactly;
            # a Sailfish tree only lands near a target, so the two are never
            # pooled silently into one marker.
            "source": str(sub["source"].iloc[0]) if "source" in sub else "synthetic",
            "r_target": r_t, "eta": float(eta), "m": int(m),
            "r_actual": float(r_of(eta, m, EXPONENT)),
            "p_star": np.nan if ps is None else float(ps),
            "reached": ps is not None,
            # Pinned at either end is a bound, not a threshold: at the floor the grid
            # never resolved the transition, at the ceiling recovery needed full data.
            "censored": ps is not None and (ps <= p.min() * 1.001 or ps >= p.max() * 0.999),
            "n_p": len(p),
        })
    return pd.DataFrame(rows).sort_values(["r_target", "m"])


def fit_exponent(tab: pd.DataFrame) -> tuple:
    """Slope ``a`` in ``p* = C (1+eta)^a log m / m``, and the per-point spread.

    Fitted across every measured cell at once, which is the one thing the paper's
    per-eta panels cannot do: each of those absorbs the eta-dependence into its own
    constant, so nothing there compares across eta.
    """
    ok = tab[tab.reached & ~tab.censored]
    if len(ok) < 3:
        return float("nan"), float("nan"), 0
    y = np.log(ok.p_star.to_numpy() * ok.m.to_numpy() / np.log(ok.m.to_numpy()))
    x = np.log1p(ok.eta.to_numpy(float))
    a, b = np.polyfit(x, y, 1)
    ratios = np.exp(y - (a * x + b))
    return float(a), float(ratios.max() / ratios.min()), len(ok)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=IN.parent)
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    args = ap.parse_args()

    if not IN.exists():
        print(f"no sweep output at {IN} -- run scripts/sweep_r_grid.py first")
        return 1

    df = pd.read_csv(IN)
    tab = pstar_table(df, args.threshold)
    args.out.mkdir(parents=True, exist_ok=True)
    tab.to_csv(args.out / "pstar.csv", index=False)

    print(f"NMI >= {args.threshold}, {len(tab)} cells swept\n")
    for r_t, sub in tab.groupby("r_target"):
        print(f"  r={r_t:.2f}")
        for _, q in sub.iterrows():
            mark = "" if q.reached and not q.censored else "  (bound, not threshold)"
            ps = "never" if not q.reached else f"{q.p_star:.5f}"
            print(f"    {q.source:<9} eta={q.eta:5.2f}  m={q.m:5d}  "
                  f"r={q.r_actual:.4f}  p*={ps:>8}{mark}")
        meas = sub[sub.reached & ~sub.censored]
        if len(meas) >= 2:
            spread = meas.p_star.max() / meas.p_star.min()
            verdict = "FLAT" if spread < 1.5 else "SLOPED"
            print(f"    -> p* spans {spread:.2f}x across m  [{verdict}]")
        print()

    a, spread, n = fit_exponent(tab)
    print(f"pooled fit over {n} measured cells: p* ~ (1+eta)^{a:.2f} log m / m")
    print(f"  per-point constants span {spread:.2f}x")
    print(f"  the grid was built at exponent {EXPONENT:g}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    targets = sorted(tab.r_target.unique())
    fig, axes = plt.subplots(1, len(targets), figsize=(3.4 * len(targets), 3.6),
                             sharey=True, squeeze=False)
    for ax, r_t in zip(axes[0], targets):
        sub = tab[tab.r_target == r_t].sort_values("m")
        for _, q in sub.iterrows():
            if not q.reached:
                continue
            c = ETA_COLORS[(int(round(q.eta)) - 1) % len(ETA_COLORS)]
            ax.plot(q.m, q.p_star, "o", ms=8, color=c,
                    mfc="none" if q.censored else c, mew=1.8)
            ax.annotate(rf"$\eta$={q.eta:g}", (q.m, q.p_star), fontsize=8,
                        xytext=(5, 4), textcoords="offset points", color=c)
        meas = sub[sub.reached & ~sub.censored]
        if len(meas) >= 2:
            ax.axhline(meas.p_star.median(), color="#9ca3af", ls="--", lw=1,
                       zorder=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$m$")
        ax.set_title(rf"$r={r_t:g}$", fontsize=11)
        ax.grid(True, which="both", alpha=0.3)
    axes[0][0].set_ylabel(rf"$\hat p^\star$ (NMI $\geq$ {args.threshold:g})")
    fig.tight_layout()
    out = args.out / "pstar_vs_m_per_r.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nwrote {out}\n      {args.out / 'pstar.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
