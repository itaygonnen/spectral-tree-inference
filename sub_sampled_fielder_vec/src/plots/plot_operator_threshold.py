"""Operator / threshold panels shared by the distance-vs-similarity notebooks.

Both ``analysis/supporting/distance_vs_similarity_{real,generated}.ipynb`` drew
these two figures from cells that were identical apart from the caption and the
tree count. Geometry, colours and the bootstrap seed live here so the real and
generated versions stay comparable.
"""
from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

# threshold rule -> colour, in the order the notebooks list them (sigma2, k-means, sign)
THR_PALETTE = ("#1d4ed8", "#065f46", "#b45309")


def plot_eta_histograms(
    panels: Sequence[Tuple[np.ndarray, str]],
    suptitle: str = "Reference-partition imbalance decides recoverability",
):
    """One eta histogram per operator, **each on its own x-range**.

    A shared scale hides the Griffing / L_sym shapes, which is the whole point
    of the comparison. ``panels`` is ``[(etas, operator_label), ...]``.
    """
    fig, axes = plt.subplots(1, len(panels), figsize=(16, 5))
    axes = np.atleast_1d(axes)
    for ax, (data, name) in zip(axes, panels):
        data = np.asarray(data, float)
        ax.hist(data, bins=30, range=(1, max(data.max() * 1.05, 2.0)),
                color="#9ca3af", edgecolor="white")
        ax.set_xlabel(r"$\eta$ = larger clan / smaller clan", fontsize=12)
        ax.set_title(f"{name}\n{len(data)} valid trees, "
                     f"median $\\eta$={np.median(data):.1f}", fontsize=12)
    axes[0].set_ylabel("# valid trees", fontsize=12)
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout()
    plt.show()
    return fig


def _median_ci(x, rng, n_boot: int = 2000, alpha: float = 0.05):
    x = np.asarray(x, float)
    if x.size == 0:
        return np.nan, np.nan, np.nan
    med = float(np.median(x))
    if x.size == 1:
        return med, med, med
    boots = np.median(rng.choice(x, size=(n_boot, x.size), replace=True), axis=1)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return med, float(lo), float(hi)


def plot_operator_threshold_bars(
    thr_etas: Dict[Tuple[str, str], np.ndarray],
    op_names: Sequence[str],
    thr_names: Sequence[str],
    n_total: int,
    suptitle: str,
    *,
    seed: int = 0,
):
    """Median eta with 95% bootstrap CI (left) and valid-partition count (right).

    Each panel is sorted high -> low **independently** — the ranking by
    imbalance is not the ranking by how often the rule fires at all, and reading
    one off the other is the error these two panels exist to prevent.
    """
    rng = np.random.default_rng(seed)
    thr_color = dict(zip(thr_names, THR_PALETTE))

    recs = []
    for op in op_names:
        for thr in thr_names:
            x = thr_etas[(op, thr)]
            med, lo, hi = _median_ci(x, rng)
            recs.append({"lab": f"{op} × {thr}", "thr": thr,
                         "med": med, "lo": lo, "hi": hi, "k": len(x)})

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(18, 6))

    # LEFT — median eta (95% CI)
    sL = sorted([r for r in recs if np.isfinite(r["med"])],
                key=lambda r: r["med"], reverse=True)
    xL = np.arange(len(sL))
    medsL = np.array([r["med"] for r in sL])
    losL = np.array([r["lo"] for r in sL])
    hisL = np.array([r["hi"] for r in sL])
    axL.bar(xL, medsL, yerr=np.vstack([medsL - losL, hisL - medsL]), capsize=4,
            color=[thr_color[r["thr"]] for r in sL], edgecolor="#374151")
    padL = 0.05 * float(hisL.max())
    for xi, r in zip(xL, sL):
        axL.text(xi, r["hi"] + padL, f"{r['med']:.1f}", ha="center", va="bottom",
                 fontsize=14, fontweight="bold")
    axL.set_ylim(0, float(hisL.max()) * 1.28)
    axL.set_xticks(xL)
    axL.set_xticklabels([r["lab"] for r in sL], rotation=30, ha="right")
    axL.set_ylabel(r"median partition imbalance $\eta$")
    axL.set_title("Median $\\eta$ (95% bootstrap CI)", pad=12)
    axL.grid(True, axis="y", alpha=0.3)

    # RIGHT — # valid partitions
    sR = sorted(recs, key=lambda r: r["k"], reverse=True)
    xR = np.arange(len(sR))
    ksR = np.array([r["k"] for r in sR])
    axR.bar(xR, ksR, color=[thr_color[r["thr"]] for r in sR], edgecolor="#374151")
    padR = 0.02 * n_total
    for xi, r in zip(xR, sR):
        axR.text(xi, r["k"] + padR, f"{r['k']}\n({100 * r['k'] / n_total:.0f}%)",
                 ha="center", va="bottom", fontsize=14, fontweight="bold")
    axR.set_ylim(0, n_total * 1.25)
    axR.set_xticks(xR)
    axR.set_xticklabels([r["lab"] for r in sR], rotation=30, ha="right")
    axR.set_ylabel(f"# valid partitions (of {n_total})")
    axR.set_title("Valid partitions", pad=12)
    axR.grid(True, axis="y", alpha=0.3)

    fig.legend(handles=[Patch(color=thr_color[t], label=t) for t in thr_names],
               title="threshold", loc="upper center",
               bbox_to_anchor=(0.5, 0.92), ncol=3, frameon=False)
    fig.suptitle(suptitle, fontsize=14, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    plt.show()
    return fig
