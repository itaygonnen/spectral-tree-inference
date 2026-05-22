"""Griffing-sweep overlay across n.

Two panels:
  * Left:  raw ``||D − D̂||_2`` vs p (log-y) — included for direct
           comparison with figures 8/9. Same dimensional artifact (grows
           with n).
  * Right: **recovery %** (sign-agnostic agreement of leading eigvec of
           ``J D̂ J`` vs leading eigvec of ``J D J``) vs p, linear y in
           [50, 100]. This is the focused panel — does subsampling
           preserve the dominant partition?

Usage:
    python plot_griffing_overlay.py <sweep_dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.griffing_io import load_griffing_results


RECOVERY_THRESHOLDS = [0.95, 0.9, 0.8]


def find_plateau_threshold(p_vals: np.ndarray, recovery: np.ndarray,
                           threshold: float) -> float | None:
    """Smallest p such that recovery(p') >= threshold for all p' >= p.

    Same convention as the NJ overlay's plateau finder — guards against
    IPW spikes at very small p that can coincidentally give high recovery
    on degenerate eigenvectors.
    """
    order = np.argsort(p_vals)[::-1]
    p_sorted = p_vals[order]
    r_sorted = recovery[order]
    last_good = None
    for p, r in zip(p_sorted, r_sorted):
        if r >= threshold:
            last_good = p
        else:
            break
    return last_good


def _plot_overlay(sweep_dir: Path, out_path: Path) -> None:
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "griffing_meta.json").exists()]
    subs = sorted(subs, key=lambda d: load_griffing_results(d)["n_taxa"])
    loaded = [load_griffing_results(s) for s in subs]
    n_vals = [r["n_taxa"] for r in loaded]

    fig, (ax_s, ax_r) = plt.subplots(1, 2, figsize=(13, 5))
    cmap = plt.get_cmap("viridis")

    for r in loaded:
        n = r["n_taxa"]
        color = cmap(n_vals.index(n) / max(1, len(n_vals) - 1))
        per_p = sorted(r["per_p"], key=lambda e: e["p"])
        p = np.array([e["p"] for e in per_p])
        sn_mean = np.array([e["spec_norm_mean"] for e in per_p])
        sn_std = np.array([e["spec_norm_std"] for e in per_p])
        rec_mean = 100.0 * np.array([e["recovery_mean"] for e in per_p])
        rec_std = 100.0 * np.array([e["recovery_std"] for e in per_p])

        ax_s.errorbar(p, sn_mean, yerr=sn_std, marker="o", capsize=3,
                      color=color, label=f"n={n}")
        ax_r.errorbar(p, rec_mean, yerr=rec_std, marker="s", capsize=3,
                      color=color, label=f"n={n}")

    ax_s.set_xlabel("sampling fraction $p$")
    ax_s.set_ylabel(r"$\|D - \hat D\|_2$")
    ax_s.set_title("Distance spectral-norm error (raw — grows with $n$)")
    ax_s.set_yscale("log")
    ax_s.grid(True, which="both", alpha=0.3)
    ax_s.legend(fontsize=8)

    ax_r.set_xlabel("sampling fraction $p$")
    ax_r.set_ylabel("recovery % (sign-agnostic)")
    ax_r.set_title(r"Leading eigvec of $J\hat D J$ vs $JDJ$ — partition recovery")
    ax_r.axhline(95.0, color="grey", ls=":", lw=1, label="95% level")
    ax_r.set_ylim(48, 102)
    ax_r.grid(True, alpha=0.3)
    ax_r.legend(fontsize=8)

    tree_models = sorted({str(r.get("tree_model", "?")) for r in loaded})
    tree_label = ", ".join(tree_models)
    imp_set = sorted({str(r.get("imputation", "zero")) for r in loaded})
    imp_label = ", ".join(imp_set)
    fig.suptitle(
        f"Griffing partition recovery: overlay across $n$  —  "
        f"tree model: {tree_label}  |  imputation: {imp_label}",
        fontsize=11,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def _compute_p_stars(sweep_dir: Path) -> List[Dict]:
    """For each (n, threshold) report the smallest p with mean recovery >= threshold
    AND every higher p also at >= threshold (plateau definition)."""
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "griffing_meta.json").exists()]
    subs = sorted(subs, key=lambda d: load_griffing_results(d)["n_taxa"])
    rows = []
    for s in subs:
        r = load_griffing_results(s)
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
        print("usage: plot_griffing_overlay.py <sweep_dir>", file=sys.stderr)
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
    _plot_overlay(sweep_dir, sweep_dir / "griffing_overlay.png")

    summary = [{
        "n": r["n"],
        "p_star_by_threshold": {f"{int(t*100):d}": r["p_star_by_threshold"][t]
                                 for t in RECOVERY_THRESHOLDS},
    } for r in rows]
    with open(sweep_dir / "griffing_p_star_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"wrote {sweep_dir / 'griffing_p_star_summary.json'}")


if __name__ == "__main__":
    main()
