"""Plot NMI recovery curves from the real-data 600-tree benchmark.

Usage
-----
    python scripts/plot_real_data_benchmark.py \\
        --results-dir results/runs/real_data_benchmark \\
        [--out-dir figures/real_data] \\
        [--metric nmi_gt]

Produces:
  1. mean_nmi_curves.pdf  — mean NMI ± std vs p, 200-taxa and 1000-taxa on same axes
  2. eta_binned_nmi.pdf   — NMI curves stratified by imbalance η
  3. phase_transition.pdf — p* comparison: empirical (NMI=0.5 crossing) vs theoretical
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------

def load_agg(results_dir: Path, n_str: str) -> dict | None:
    path = results_dir / f"n{n_str}" / "aggregate_results.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Theory: approximate p* from CBM (λ_min / η analogy)
# ---------------------------------------------------------------------------

def theoretical_p_star(n: int, eta: float) -> float:
    """Rough p* ≈ C * η / sqrt(n) from the main theorem (placeholder scaling)."""
    C = 4.0  # constant from Lemma 0.4 — tune after seeing data
    return min(C * eta / np.sqrt(n), 1.0)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

COLORS = {
    "200":  "#1d4ed8",   # blue
    "1000": "#b45309",   # amber
    "balanced":  "#065f46",
    "moderate":  "#1d4ed8",
    "imbalanced": "#991b1b",
}

LINESTYLES = {
    "200": "-",
    "1000": "--",
}


def plot_mean_curves(results_dir: Path, out_dir: Path, metric: str = "nmi_gt") -> None:
    fig, ax = plt.subplots(figsize=(7, 7))

    for n_str in ["200", "1000"]:
        agg = load_agg(results_dir, n_str)
        if agg is None:
            continue
        p = [row["p"] for row in agg["global_per_p"]]
        mean = [row[f"{metric}_mean"] for row in agg["global_per_p"]]
        std  = [row[f"{metric}_std"]  for row in agg["global_per_p"]]
        p = np.array(p); mean = np.array(mean); std = np.array(std)

        color = COLORS[n_str]
        ls = LINESTYLES[n_str]
        ax.plot(p, mean, color=color, ls=ls, lw=2,
                label=f"n={n_str} (N={agg['n_trees']} trees)")
        ax.fill_between(p, mean - std, mean + std, color=color, alpha=0.15)

    ax.axhline(0.5, color="#6b7280", ls=":", lw=1)
    ax.set_xlabel("Sub-sampling fraction $p$", fontsize=13)
    ylabel = "NMI (vs ground-truth bipartition)" if "nmi" in metric else metric.upper()
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(f"B-matrix clan recovery — real data ({metric})", fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1.05)
    fig.tight_layout()
    out_path = out_dir / "mean_nmi_curves.pdf"
    fig.savefig(str(out_path))
    print(f"Saved: {out_path}")
    plt.close(fig)


def plot_eta_binned(results_dir: Path, out_dir: Path, n_str: str = "200",
                    metric: str = "nmi_gt") -> None:
    agg = load_agg(results_dir, n_str)
    if agg is None:
        print(f"No results for n={n_str}")
        return

    fig, ax = plt.subplots(figsize=(7, 7))
    for bin_data in agg.get("eta_binned", []):
        label = bin_data["label"]
        n_bin = bin_data["n_trees"]
        p = np.array([row["p"] for row in bin_data["per_p"]])
        mean = np.array([row[f"{metric}_mean"] for row in bin_data["per_p"]])
        std  = np.array([row[f"{metric}_std"]  for row in bin_data["per_p"]])
        color = COLORS.get(label, "#374151")
        ax.plot(p, mean, color=color, lw=2,
                label=f"{label} (N={n_bin})")
        ax.fill_between(p, mean - std, mean + std, color=color, alpha=0.15)

    ax.axhline(0.5, color="#6b7280", ls=":", lw=1)
    ax.set_xlabel("Sub-sampling fraction $p$", fontsize=13)
    ax.set_ylabel(f"Mean NMI (vs ground truth)", fontsize=13)
    ax.set_title(f"NMI by imbalance $\\eta$ — n={n_str} taxa", fontsize=14)
    ax.legend(fontsize=11)
    ax.set_xlim(0, 1); ax.set_ylim(-0.05, 1.05)
    fig.tight_layout()
    out_path = out_dir / f"eta_binned_nmi_n{n_str}.pdf"
    fig.savefig(str(out_path))
    print(f"Saved: {out_path}")
    plt.close(fig)


def plot_phase_transition(results_dir: Path, out_dir: Path,
                           metric: str = "nmi_gt") -> None:
    """Scatter: empirical p* (NMI=0.5 crossing) vs theoretical p*."""
    from src.runners.real_data_bpart import _eta_from_partition  # noqa: F401

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, n_str in zip(axes, ["200", "1000"]):
        per_tree_path = results_dir / f"n{n_str}" / "per_tree_results.json"
        if not per_tree_path.exists():
            ax.set_title(f"n={n_str}: no data")
            continue
        with open(per_tree_path) as f:
            data = json.load(f)
        trees = data["trees"]

        emp_p_star, theo_p_star, etas = [], [], []
        for t in trees:
            p_vals = [row["p"] for row in t["per_p"]]
            nmi_means = [row[f"{metric}_mean"] for row in t["per_p"]]
            # Empirical p*: smallest p where NMI >= 0.5
            crossings = [p for p, nmi in zip(p_vals, nmi_means) if nmi >= 0.5]
            if not crossings:
                continue
            emp = float(min(crossings))
            eta = float(t["eta"])
            n = int(t["n_taxa"])
            theo = theoretical_p_star(n, eta)
            emp_p_star.append(emp)
            theo_p_star.append(theo)
            etas.append(eta)

        if not emp_p_star:
            ax.set_title(f"n={n_str}: insufficient crossings")
            continue

        sc = ax.scatter(theo_p_star, emp_p_star, c=np.log(etas),
                        cmap="RdBu_r", alpha=0.6, s=20)
        fig.colorbar(sc, ax=ax, label=r"$\log\,\eta$")
        lim = max(max(theo_p_star), max(emp_p_star)) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.5)
        ax.set_xlabel(r"Theoretical $p^*$", fontsize=12)
        ax.set_ylabel(r"Empirical $p^*$ (NMI$\geq$0.5)", fontsize=12)
        ax.set_title(f"Phase transition — n={n_str} taxa", fontsize=13)

    fig.tight_layout()
    out_path = out_dir / "phase_transition.pdf"
    fig.savefig(str(out_path))
    print(f"Saved: {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default="results/runs/real_data_benchmark")
    parser.add_argument("--out-dir", default="figures/real_data")
    parser.add_argument("--metric", default="nmi_gt",
                        choices=["nmi_gt", "ari_gt", "agr_gt",
                                 "nmi_ref", "ari_ref", "agr_ref"])
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_mean_curves(results_dir, out_dir, args.metric)
    for n_str in ["200", "1000"]:
        plot_eta_binned(results_dir, out_dir, n_str, args.metric)
    plot_phase_transition(results_dir, out_dir, args.metric)


if __name__ == "__main__":
    main()
