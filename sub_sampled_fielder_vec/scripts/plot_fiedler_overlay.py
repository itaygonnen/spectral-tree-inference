"""First-layer Fiedler partition-agreement overlay.

Matches the visual style of
``02_real_data_sweeps/kingman_threshold_vs_theory.ipynb``:

  * Left:  ``partition_agreement_M`` (%) vs $p$, **log-x**, one curve per
           $n$ (viridis), red dashed 95% threshold. Each $n$ is labelled
           with its empirical $\\hat p^*$.
  * Right: $\\hat p^*(n)$ vs $n$ (**log-log**). Empirical only here — a
           theory overlay (non-balanced CBM bound) needs cached
           ``similarity_matrix.npz`` / ``fiedler_ref.npz`` from a
           production sweep; not produced by this lightweight runner.

Usage:
    python plot_fiedler_overlay.py <sweep_dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.fiedler_io import load_fiedler_results


RECOVERY_THRESHOLD = 0.95
RECOVERY_THRESHOLDS = [0.95, 0.9, 0.8]


def find_plateau_threshold(p_vals: np.ndarray, agreement_frac: np.ndarray,
                           threshold: float) -> float | None:
    """Smallest p such that agreement >= threshold for all p' >= p.

    Guards against IPW low-p spikes where a degenerate Fiedler can
    coincidentally match the reference.
    """
    order = np.argsort(p_vals)[::-1]
    p_sorted = p_vals[order]
    a_sorted = agreement_frac[order]
    last_good = None
    for p, a in zip(p_sorted, a_sorted):
        if a >= threshold:
            last_good = p
        else:
            break
    return last_good


def _plot_overlay(sweep_dir: Path, out_path: Path) -> None:
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "fiedler_meta.json").exists()]
    subs = sorted(subs, key=lambda d: load_fiedler_results(d)["n_taxa"])
    loaded = [load_fiedler_results(s) for s in subs]
    n_vals = [r["n_taxa"] for r in loaded]
    tree_models = sorted({str(r.get("tree_model", "?")) for r in loaded})
    tree_label = ", ".join(tree_models)
    seq_len = next((r.get("seq_len") for r in loaded if r.get("seq_len")), "?")

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f"First-layer Fiedler — tree={tree_label}, L={seq_len}, "
        f"sampling=uniform IPW, metric=partition_agreement_M"
    )

    cmap = plt.get_cmap("viridis")
    n_to_color = {n: cmap(i / max(1, len(n_vals) - 1)) for i, n in enumerate(n_vals)}

    # ---------- Left: recovery curves ----------
    p_star_per_n: Dict[int, float | None] = {}
    for r in loaded:
        n = r["n_taxa"]
        per_p = sorted(r["per_p"], key=lambda e: e["p"])
        p = np.array([e["p"] for e in per_p])
        agreement_frac = np.array([e["recovery_mean"] for e in per_p])
        # 'recovery_mean' in the runner JSON is already max(s, 1-s) in [0.5, 1.0]
        agreement_pct = 100.0 * agreement_frac

        p_star = find_plateau_threshold(p, agreement_frac, RECOVERY_THRESHOLD)
        p_star_per_n[n] = p_star
        label_p = (f"{p_star:.4f}"
                   if p_star is not None and np.isfinite(p_star) else "n/a")
        ax_l.plot(
            p, agreement_pct, "-o", color=n_to_color[n], ms=4, lw=1,
            label=fr"n={n}, $\hat p^*$={label_p}",
        )

    ax_l.axhline(100.0 * RECOVERY_THRESHOLD, color="red", ls="--", lw=1,
                 label=f"{int(RECOVERY_THRESHOLD * 100)}% threshold")
    ax_l.set_xscale("log")
    ax_l.set_xlabel(r"Sampling probability $p$")
    ax_l.set_ylabel("partition_agreement_M (%)")
    ax_l.set_title(fr"Recovery curves by $n$ — {tree_label}")
    ax_l.set_ylim(48, 102)
    ax_l.legend(fontsize=8, loc="lower right")
    ax_l.grid(True, which="both", alpha=0.3)

    # ---------- Right: p*(n) scaling ----------
    n_arr = np.array([n for n in n_vals if p_star_per_n[n] is not None
                       and np.isfinite(p_star_per_n[n])], dtype=float)
    pstar_arr = np.array([p_star_per_n[int(n)] for n in n_arr], dtype=float)
    if n_arr.size > 0:
        ax_r.plot(n_arr, pstar_arr, "-o", label=r"Empirical $\hat p^*$")
    ax_r.set_xscale("log")
    ax_r.set_yscale("log")
    ax_r.set_xlabel(r"$n$ (number of taxa)")
    ax_r.set_ylabel(r"$\hat p^*$")
    ax_r.set_title(fr"Scaling: $\hat p^*(n)$ — {tree_label}")
    ax_r.legend(fontsize=9, loc="best")
    ax_r.grid(True, which="both", alpha=0.3)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def _compute_p_stars(sweep_dir: Path) -> List[Dict]:
    """Per (n, threshold): smallest p with partition_agreement >= threshold AND
    every higher p also at that level (plateau definition)."""
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "fiedler_meta.json").exists()]
    subs = sorted(subs, key=lambda d: load_fiedler_results(d)["n_taxa"])
    rows = []
    for s in subs:
        r = load_fiedler_results(s)
        n = int(r["n_taxa"])
        per_p = sorted(r["per_p"], key=lambda e: e["p"])
        p_arr = np.array([e["p"] for e in per_p])
        rec_arr = np.array([e["recovery_mean"] for e in per_p])
        thr_p_star = {
            t: find_plateau_threshold(p_arr, rec_arr, t)
            for t in RECOVERY_THRESHOLDS
        }
        rows.append({"n": n, "p_star_by_threshold": thr_p_star})
    return rows


def main():
    if len(sys.argv) < 2:
        print("usage: plot_fiedler_overlay.py <sweep_dir>", file=sys.stderr)
        sys.exit(1)
    sweep_dir = Path(sys.argv[1]).resolve()
    rows = _compute_p_stars(sweep_dir)
    header = f"{'n':>5}  " + "  ".join(f">={int(t*100)}%" for t in RECOVERY_THRESHOLDS)
    print(header)
    for r in rows:
        cells = [f"{r['n']:>5}"]
        for t in RECOVERY_THRESHOLDS:
            p = r["p_star_by_threshold"].get(t)
            cells.append(f"{p:>5.3f}" if (p is not None and np.isfinite(p)) else "  n/a")
        print("  ".join(cells))
    _plot_overlay(sweep_dir, sweep_dir / "fiedler_overlay.png")

    summary = [{
        "n": r["n"],
        "p_star_by_threshold": {f"{int(t*100):d}": r["p_star_by_threshold"][t]
                                 for t in RECOVERY_THRESHOLDS},
    } for r in rows]
    with open(sweep_dir / "fiedler_p_star_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"wrote {sweep_dir / 'fiedler_p_star_summary.json'}")


if __name__ == "__main__":
    main()
