"""Per-$n$ / per-$\\eta$ subplot grids for the bpart η-pool / synthetic sweeps.

Reads ``bpart_eta_pool_n*.json`` (the schema written by
``runners.bpart_eta_pool`` and ``runners.bpart_synthetic``) and renders, for a
chosen metric, either layout used in
``02_real_data_sweeps/eta_pool_sweep_per_n_median.ipynb``:

  * ``plot_per_n``  — one subplot per ``n``, one curve per η target  (Fig 5 style)
  * ``plot_per_eta`` — one subplot per η target, one curve per ``n`` (Fig 6 style)

This replaces the single all-in-one overlay (``plot_bpart_eta_pool``) when the
(n × η) grid is too dense to read in one axes.

Usage:
    python plot_bpart_eta_grid.py <dir> [metric]   # metric ∈ agreement|dot|ari|nmi
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt

# metric -> (mean_field, std_field, ylabel, ylim, reference_line)
_METRICS = {
    "agreement": ("agreement_mean", "agreement_std", "clan partition agreement (%)", (45, 102), 50.0),
    "dot":       ("dot_mean", "dot_std", r"$|\langle \hat v_1, v_1\rangle|$", (-0.02, 1.02), None),
    "ari":       ("ari_mean", "ari_std", "Adjusted Rand Index", (-0.05, 1.05), 0.0),
    "nmi":       ("nmi_mean", "nmi_std", "Normalized Mutual Information", (-0.05, 1.05), 0.0),
}


def _load(dir_path: Path) -> List[Dict]:
    out = [json.loads(p.read_text()) for p in sorted(Path(dir_path).glob("bpart_eta_pool_n*.json"))]
    if not out:
        raise SystemExit(f"No bpart_eta_pool_n*.json under {dir_path}")
    return sorted(out, key=lambda r: r["n"])


def _targets(loaded: List[Dict]) -> List[int]:
    return sorted({int(t) for r in loaded for t in r["bins"]})


def _grid(n_panels: int):
    ncols = 2
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows),
                             sharey=True, squeeze=False)
    return fig, axes


def _per_p(b: Dict):
    pp = sorted(b["per_p"], key=lambda e: e["p"])
    return [e["p"] for e in pp], pp


def _finish(ax, ylabel, ylim, ref) -> None:
    if ref is not None:
        ax.axhline(ref, color="grey", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel(r"sampling fraction $p$")
    ax.set_ylim(*ylim)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)


def _save(fig, out_path: Path) -> Path:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_per_n(dir_path, out_path, metric: str = "agreement", title_extra: str = "") -> Path:
    """One subplot per ``n``; curves coloured by η target."""
    loaded = _load(Path(dir_path))
    mk, _sk, ylabel, ylim, ref = _METRICS[metric]
    targets = _targets(loaded)
    cmap = plt.get_cmap("viridis")
    color_of = {t: cmap(i / max(1, len(targets) - 1)) for i, t in enumerate(targets)}

    fig, axes = _grid(len(loaded))
    flat = axes.flatten()
    for ax, r in zip(flat, loaded):
        for t in targets:
            b = r["bins"].get(str(t))
            if b is None:
                continue
            p, pp = _per_p(b)
            ax.plot(p, [e[mk] for e in pp], "-o", ms=4, lw=1.5, color=color_of[t],
                    label=rf"$\eta$={b['eta_mean']:.1f}")
        ax.set_title(rf"$n$ = {r['n']}")
        _finish(ax, ylabel, ylim, ref)
    for ax in flat[len(loaded):]:
        ax.set_visible(False)
    for row in axes:
        row[0].set_ylabel(ylabel)
    fig.suptitle(rf"{metric} vs $p$ — one panel per $n$, curves by $\eta${title_extra}")
    return _save(fig, out_path)


def plot_per_eta(dir_path, out_path, metric: str = "agreement", title_extra: str = "") -> Path:
    """One subplot per η target; curves coloured by ``n``."""
    loaded = _load(Path(dir_path))
    mk, _sk, ylabel, ylim, ref = _METRICS[metric]
    targets = _targets(loaded)
    n_vals = [r["n"] for r in loaded]
    cmap = plt.get_cmap("viridis")
    color_of = {n: cmap(i / max(1, len(n_vals) - 1)) for i, n in enumerate(n_vals)}

    fig, axes = _grid(len(targets))
    flat = axes.flatten()
    for ax, t in zip(flat, targets):
        for r in loaded:
            b = r["bins"].get(str(t))
            if b is None:
                continue
            p, pp = _per_p(b)
            ax.plot(p, [e[mk] for e in pp], "-o", ms=4, lw=1.5, color=color_of[r["n"]],
                    label=rf"$n$={r['n']} ($\eta$={b['eta_mean']:.1f})")
        ax.set_title(rf"$\eta \approx {t}$")
        _finish(ax, ylabel, ylim, ref)
    for ax in flat[len(targets):]:
        ax.set_visible(False)
    for row in axes:
        row[0].set_ylabel(ylabel)
    fig.suptitle(rf"{metric} vs $p$ — one panel per $\eta$, curves by $n${title_extra}")
    return _save(fig, out_path)


def main():
    if len(sys.argv) < 2:
        print("usage: plot_bpart_eta_grid.py <dir> [metric]", file=sys.stderr)
        sys.exit(1)
    d = Path(sys.argv[1]).resolve()
    metric = sys.argv[2] if len(sys.argv) > 2 else "agreement"
    plot_per_n(d, d / f"bpart_per_n_{metric}.png", metric=metric)
    plot_per_eta(d, d / f"bpart_per_eta_{metric}.png", metric=metric)


if __name__ == "__main__":
    main()
