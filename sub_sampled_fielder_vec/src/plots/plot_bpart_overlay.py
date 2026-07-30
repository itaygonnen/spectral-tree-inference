"""Plot the distance-based (B-matrix) clan-partition sweep across n.

Reads all ``n*_L*/bpart_meta.json`` under a sweep directory and overlays, per
n, the partition agreement and leading-eigenvector dot product vs sampling
fraction p. Companion to ``plot_nj_p_star_vs_n.py`` for the B-matrix method
(double-centred distance, STDR Appendix D).

Usage:
    python plot_bpart_overlay.py <sweep_dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt


def load_bpart_results(run_dir: Path) -> Dict:
    with open(Path(run_dir) / "bpart_meta.json") as f:
        return json.load(f)


def _load_sweep(sweep_dir: Path) -> List[Dict]:
    out = [load_bpart_results(sub)
           for sub in sorted(sweep_dir.iterdir())
           if sub.is_dir() and (sub / "bpart_meta.json").exists()]
    if not out:
        raise SystemExit(f"No bpart_meta.json found under {sweep_dir}")
    return sorted(out, key=lambda r: r["n_taxa"])


def _plot_overlay(sweep_dir: Path, out_path: Path) -> None:
    loaded = _load_sweep(sweep_dir)
    n_vals = [r["n_taxa"] for r in loaded]
    cmap = plt.get_cmap("viridis")

    # x = p on a log scale: full data (p=1) on the right, recovery rises to it.
    fig, ((ax_a, ax_d), (ax_ari, ax_nmi)) = plt.subplots(2, 2, figsize=(13, 10))
    for r in loaded:
        color = cmap(n_vals.index(r["n_taxa"]) / max(1, len(n_vals) - 1))
        per_p = sorted(r["per_p"], key=lambda e: e["p"])
        p = np.array([e["p"] for e in per_p])
        label = f"n={r['n_taxa']}"
        if r.get("eta") is not None:
            label += rf", $\eta$={r['eta']:.2f}"
        for ax, mean_key, std_key, marker in (
            (ax_a, "agreement_mean", "agreement_std", "o"),
            (ax_d, "dot_mean", "dot_std", "d"),
            (ax_ari, "ari_mean", "ari_std", "s"),
            (ax_nmi, "nmi_mean", "nmi_std", "^"),
        ):
            mean = np.array([e[mean_key] for e in per_p])
            std = np.array([e.get(std_key, 0.0) for e in per_p])
            ax.errorbar(p, mean, yerr=std, marker=marker, capsize=3,
                        color=color, label=label)

    ax_a.axhline(50.0, color="grey", linestyle=":", linewidth=1,
                 label="50% (random)")
    ax_a.set_ylabel("clan partition agreement (%)")
    ax_a.set_title(r"B-matrix clan agreement vs $p$, per $n$")
    ax_a.set_ylim(45, 102)

    ax_d.set_ylabel(r"$|\langle \hat v_1, v_1\rangle|$ (leading eigvec)")
    ax_d.set_title(r"Leading-eigenvector alignment vs $p$, per $n$")
    ax_d.set_ylim(-0.02, 1.02)

    ax_ari.set_ylabel("Adjusted Rand Index")
    ax_ari.set_title(r"ARI of sign bipartition vs $p$, per $n$")
    ax_ari.set_ylim(-0.05, 1.05)

    ax_nmi.set_ylabel("Normalized Mutual Information")
    ax_nmi.set_title(r"NMI of sign bipartition vs $p$, per $n$")
    ax_nmi.set_ylim(-0.05, 1.05)

    for ax in (ax_a, ax_d, ax_ari, ax_nmi):
        ax.set_xscale("log")
        ax.set_xlabel(r"sampling fraction $p$")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8, ncol=2)

    tree_models = sorted({str(r.get("tree_model", "?")) for r in loaded})
    fig.suptitle(
        "Distance double-centering ($B=(I-11^T/m)D(I-11^T/m)$) clan partition "
        f"under subsampling  —  tree model: {', '.join(tree_models)}",
        fontsize=11,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_bpart_overlay.py <sweep_dir>", file=sys.stderr)
        sys.exit(1)
    sweep_dir = Path(sys.argv[1]).resolve()
    _plot_overlay(sweep_dir, sweep_dir / "bpart_overlay.png")


if __name__ == "__main__":
    main()
