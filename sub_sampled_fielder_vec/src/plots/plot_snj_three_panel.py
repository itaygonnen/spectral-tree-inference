"""Three-panel SNJ-paper analog plot.

Panel 1: ||R - R̂||_2 vs p              (mean ± std)
Panel 2: σ₂ separation (Fig. 3 analog) — KDE/hist per p, adj vs non-adj
Panel 3: RF distance vs p              (mean and success-rate)

Usage:
    python plot_snj_three_panel.py <results_dir>

The results dir is expected to contain one or more ``n{n}_L{L}/snj_results.json``
subdirectories. One figure is produced per (n, L).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.snj_io import load_snj_results


def _load_results(run_dir: Path) -> Dict:
    return load_snj_results(run_dir)


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
    rf_success = np.array([entry["rf_success_rate"] for entry in per_p])
    n_taxa = int(results.get("n_taxa", 0))
    max_rf = max(1, 2 * (n_taxa - 3))
    nrf_mean = rf_mean / max_rf
    nrf_std = rf_std / max_rf
    branch_correctness = 100.0 * (1.0 - nrf_mean)

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.3], hspace=0.35, wspace=0.3)

    # Panel 1: spec norm
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.errorbar(p_values, spec_norm_mean, yerr=spec_norm_std, marker="o",
                 capsize=4, color="C0", label=r"$\|R - \hat R\|_2$")
    ax1.set_xlabel("sampling fraction $p$")
    ax1.set_ylabel(r"$\|R - \hat R\|_2$")
    ax1.set_title("Spectral-norm error (Thm 4.2)")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Panel 3: normalised RF + branch correctness
    ax3 = fig.add_subplot(gs[0, 1])
    ax3.errorbar(p_values, nrf_mean, yerr=nrf_std, marker="s", capsize=4,
                 color="C3", label=r"normalised RF $= \mathrm{RF}/(2(n-3))$")
    ax3.set_xlabel("sampling fraction $p$")
    ax3.set_ylabel("normalised RF")
    ax3.set_ylim(-0.02, 1.02)
    ax3.set_title("Normalised RF vs $p$")
    ax3.grid(True, alpha=0.3)

    ax3_right = ax3.twinx()
    ax3_right.plot(p_values, branch_correctness, marker="d", color="C2",
                   linestyle="--", label=r"branch correctness $100(1-\mathrm{nRF})$")
    ax3_right.set_ylabel("branch correctness (%)", color="C2")
    ax3_right.set_ylim(-5, 105)
    ax3_right.tick_params(axis="y", labelcolor="C2")
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_right.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc="center left")

    # Panel 2: σ₂ separation — one subplot per p
    n_p = len(p_values)
    inner_gs = gs[1, :].subgridspec(1, n_p, wspace=0.25)
    log_floor = 1e-8
    for idx, entry in enumerate(per_p):
        ax = fig.add_subplot(inner_gs[0, idx])
        adj = np.asarray(entry["sigma2_adj"], dtype=float)
        nonadj = np.asarray(entry["sigma2_nonadj"], dtype=float)
        # Plot log10(σ₂) to mirror SNJ Fig. 3.
        adj_log = np.log10(np.maximum(adj, log_floor))
        nonadj_log = np.log10(np.maximum(nonadj, log_floor))
        all_log = np.concatenate([adj_log, nonadj_log]) if (adj_log.size + nonadj_log.size) else np.array([0.0])
        bins = np.linspace(all_log.min() - 0.1, all_log.max() + 0.1, 40)
        if adj_log.size:
            ax.hist(adj_log, bins=bins, color="C0", alpha=0.6,
                    label="adjacent", density=True)
        if nonadj_log.size:
            ax.hist(nonadj_log, bins=bins, color="C1", alpha=0.6,
                    label="non-adj.", density=True)
        ax.set_title(f"$p={entry['p']:.2g}$", fontsize=10)
        ax.set_xlabel(r"$\log_{10}\sigma_2$")
        if idx == 0:
            ax.set_ylabel("density")
            ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"SNJ-subsampling: n={results['n_taxa']}, L={results['seq_len']}, "
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
        print("usage: plot_snj_three_panel.py <results_dir>", file=sys.stderr)
        sys.exit(1)
    base_dir = Path(sys.argv[1]).resolve()
    subdirs = sorted([d for d in base_dir.iterdir()
                      if d.is_dir() and (
                          (d / "snj_meta.json").exists()
                          or (d / "snj_results.json").exists())])
    if not subdirs:
        print(f"No SNJ results found under {base_dir}", file=sys.stderr)
        sys.exit(2)
    for sd in subdirs:
        results = _load_results(sd)
        out_path = sd / "snj_three_panel.png"
        _plot_three_panel(results, out_path)


if __name__ == "__main__":
    main()
