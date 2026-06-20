"""Plot the η-binned Kingman B-method sweep.

Reads ``bpart_eta_pool_n*.json`` (written by ``runners.bpart_eta_pool``) under a
directory and overlays, per metric, one curve per (n, η-bin): colour encodes n,
line style encodes the η bin (solid = smallest target, dashed = next, ...). The
legend reports the bin's realized η as mean ± std, so the figure is directly
comparable to the η-controlled theory sweeps.

Usage:
    python plot_bpart_eta_pool.py <dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt

_STYLES = ["-", "--", ":", "-."]


def _load(dir_path: Path) -> List[Dict]:
    out = [json.loads(p.read_text()) for p in sorted(dir_path.glob("bpart_eta_pool_n*.json"))]
    if not out:
        raise SystemExit(f"No bpart_eta_pool_n*.json under {dir_path}")
    return sorted(out, key=lambda r: r["n"])


def _plot(dir_path: Path, out_path: Path) -> None:
    loaded = _load(dir_path)
    n_vals = [r["n"] for r in loaded]
    targets = sorted({int(t) for r in loaded for t in r["bins"]})
    style_of = {t: _STYLES[i % len(_STYLES)] for i, t in enumerate(targets)}
    cmap = plt.get_cmap("viridis")

    fig, ((ax_a, ax_d), (ax_ari, ax_nmi)) = plt.subplots(2, 2, figsize=(13, 10))
    panels = [(ax_a, "agreement_mean", "agreement_std"),
              (ax_d, "dot_mean", "dot_std"),
              (ax_ari, "ari_mean", "ari_std"),
              (ax_nmi, "nmi_mean", "nmi_std")]

    for r in loaded:
        color = cmap(n_vals.index(r["n"]) / max(1, len(n_vals) - 1))
        for t in targets:
            b = r["bins"].get(str(t))
            if b is None:
                continue
            per_p = sorted(b["per_p"], key=lambda e: e["p"])
            p = np.array([e["p"] for e in per_p])
            label = rf"n={r['n']}, $\eta$={b['eta_mean']:.1f}$\pm${b['eta_std']:.1f}"
            for ax, mk, sk in panels:
                ax.errorbar(p, [e[mk] for e in per_p], yerr=[e[sk] for e in per_p],
                            marker="o", ms=4, capsize=2, color=color,
                            linestyle=style_of[t], label=label)

    ax_a.axhline(50.0, color="grey", linestyle=":", linewidth=1)
    ax_a.set_ylabel("clan partition agreement (%)"); ax_a.set_ylim(45, 102)
    ax_a.set_title(r"Agreement vs $p$ (Kingman, by $\eta$ bin)")
    ax_d.set_ylabel(r"$|\langle \hat v_1, v_1\rangle|$"); ax_d.set_ylim(-0.02, 1.02)
    ax_d.set_title(r"Leading-eigenvector alignment vs $p$")
    ax_ari.set_ylabel("Adjusted Rand Index"); ax_ari.set_ylim(-0.05, 1.05)
    ax_ari.set_title(r"ARI vs $p$")
    ax_nmi.set_ylabel("Normalized Mutual Information"); ax_nmi.set_ylim(-0.05, 1.05)
    ax_nmi.set_title(r"NMI vs $p$")
    for ax, _mk, _sk in panels:
        ax.set_xscale("log"); ax.set_xlabel(r"sampling fraction $p$")
        ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7, ncol=2)

    fig.suptitle(
        r"Kingman B-matrix clan partition under subsampling, by $\eta$ bin "
        r"(solid/dashed = $\eta$ target); line style $\eta\!\approx\!1$ vs $\eta\!\approx\!5$",
        fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_bpart_eta_pool.py <dir>", file=sys.stderr)
        sys.exit(1)
    d = Path(sys.argv[1]).resolve()
    _plot(d, d / "bpart_eta_pool_overlay.png")


if __name__ == "__main__":
    main()
