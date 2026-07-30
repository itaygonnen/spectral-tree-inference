"""Three-panel NJ-with-subsampling plot.

Panel 1: ||D - D̂||_2 vs p                  (mean ± std)
Panel 2: nRF + branch correctness vs p      (left axis: nRF; right axis: %)
Panel 3: Q-criterion separation per p       (raw Q histograms, adj vs non-adj)

Usage:
    python plot_nj_three_panel.py <results_dir>

``results_dir`` is expected to contain one or more ``n{n}_L{L}/nj_meta.json``
subdirectories (or legacy ``nj_results.json``). One figure is produced per
(n, L).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.nj_io import load_nj_results


def _plot_three_panel(results: Dict, out_path: Path) -> None:
    per_p = results["per_p"]
    p_values = np.array([entry["p"] for entry in per_p])
    order = np.argsort(p_values)
    p_values = p_values[order]
    per_p = [per_p[i] for i in order]

    spec_norm_mean = np.array([entry["spec_norm_mean"] for entry in per_p])
    spec_norm_std = np.array([entry["spec_norm_std"] for entry in per_p])
    rf_mean = np.array([entry["rf_mean"] for entry in per_p])
    rf_std = np.array([entry["rf_std"] for entry in per_p])
    n_taxa = int(results.get("n_taxa", 0))
    max_rf = max(1, 2 * (n_taxa - 3))
    nrf_mean = rf_mean / max_rf
    nrf_std = rf_std / max_rf
    branch_correctness = 100.0 * (1.0 - nrf_mean)

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.35, wspace=0.3)

    # Panel 1: spec norm on D
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.errorbar(p_values, spec_norm_mean, yerr=spec_norm_std, marker="o",
                 capsize=4, color="C0", label=r"$\|D - \hat D\|_2$")
    ax1.set_xlabel("sampling fraction $p$")
    ax1.set_ylabel(r"$\|D - \hat D\|_2$")
    ax1.set_title("Distance-matrix spectral-norm error")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Panel 2: nRF + branch correctness
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.errorbar(p_values, nrf_mean, yerr=nrf_std, marker="s", capsize=4,
                 color="C3", label=r"normalised RF $= \mathrm{RF}/(2(n-3))$")
    ax2.set_xlabel("sampling fraction $p$")
    ax2.set_ylabel("normalised RF")
    ax2.set_ylim(-0.02, 1.02)
    ax2.set_title("Normalised RF vs $p$ (classical NJ)")
    ax2.grid(True, alpha=0.3)

    ax2_right = ax2.twinx()
    ax2_right.plot(p_values, branch_correctness, marker="d", color="C2",
                   linestyle="--", label=r"branch correctness $100(1-\mathrm{nRF})$")
    ax2_right.set_ylabel("branch correctness (%)", color="C2")
    ax2_right.set_ylim(-5, 105)
    ax2_right.tick_params(axis="y", labelcolor="C2")
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_right.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, loc="center left")

    # Panel 3: Q separation — one subplot per p, raw Q on linear axis
    n_p = len(p_values)
    inner_gs = gs[1, :].subgridspec(1, n_p, wspace=0.25)
    for idx, entry in enumerate(per_p):
        ax = fig.add_subplot(inner_gs[0, idx])
        adj = np.asarray(entry["q_adj"], dtype=float)
        nonadj = np.asarray(entry["q_nonadj"], dtype=float)
        all_q = np.concatenate([adj, nonadj]) if (adj.size + nonadj.size) else np.array([0.0])
        if all_q.size == 0:
            continue
        lo, hi = all_q.min(), all_q.max()
        if hi <= lo:
            hi = lo + 1e-9
        bins = np.linspace(lo, hi, 40)
        if adj.size:
            ax.hist(adj, bins=bins, color="C0", alpha=0.6,
                    label="adjacent", density=True)
        if nonadj.size:
            ax.hist(nonadj, bins=bins, color="C1", alpha=0.6,
                    label="non-adj.", density=True)
        ax.set_title(f"$p={entry['p']:.2g}$", fontsize=10)
        ax.set_xlabel(r"$Q$")
        if idx == 0:
            ax.set_ylabel("density")
            ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle(
        f"NJ-subsampling: n={results['n_taxa']}, L={results['seq_len']}, "
        f"{results['tree_model']}, μ={results['mutation_rate']}, "
        f"reps={results['bootstrap_reps']}, sampling={results['sampling_method']}",
        fontsize=11,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_nj_three_panel.py <results_dir>", file=sys.stderr)
        sys.exit(1)
    base_dir = Path(sys.argv[1]).resolve()
    subdirs = sorted([d for d in base_dir.iterdir()
                      if d.is_dir() and (
                          (d / "nj_meta.json").exists()
                          or (d / "nj_results.json").exists())])
    if not subdirs:
        print(f"No NJ results found under {base_dir}", file=sys.stderr)
        sys.exit(2)
    for sd in subdirs:
        results = load_nj_results(sd)
        out_path = sd / "nj_three_panel.png"
        _plot_three_panel(results, out_path)


if __name__ == "__main__":
    main()
