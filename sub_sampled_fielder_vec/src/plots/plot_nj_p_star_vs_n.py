"""Plot classical-NJ critical p* as a function of n.

Reads all ``n*_L*/nj_meta.json`` files under a sweep directory and computes a
family of branch-correctness thresholds. Mirrors ``plot_snj_p_star_vs_n.py``
so SNJ and NJ scaling can be compared head-to-head.

Usage:
    python plot_nj_p_star_vs_n.py <sweep_dir>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.threshold_utils import (
    fit_power_law,
    evaluate_power_law,
)
from src.utils.nj_io import load_nj_results


def find_plateau_threshold(p_vals: np.ndarray, correctness: np.ndarray,
                           threshold: float) -> float | None:
    """Smallest p such that correctness(p') ≥ threshold for **all** p' ≥ p.

    Avoids IPW-artifact spikes at extremely small p, where a near-zero
    sub-sampled matrix with a few entries scaled by 1/p can coincidentally
    produce a non-trivial tree.
    """
    order = np.argsort(p_vals)[::-1]
    p_sorted = p_vals[order]
    c_sorted = correctness[order]
    last_good = None
    for p, c in zip(p_sorted, c_sorted):
        if c >= threshold:
            last_good = p
        else:
            break
    return last_good


P_STAR_THRESHOLD = 90.0
P_STAR_THRESHOLDS = [50.0, 70.0, 80.0, 90.0, 95.0]


def _branch_correctness_pct(rf_mean: np.ndarray, n_taxa: int) -> np.ndarray:
    max_rf = max(1, 2 * (n_taxa - 3))
    return 100.0 * (1.0 - rf_mean / max_rf)


def _load_sweep(sweep_dir: Path) -> Dict[int, List[Tuple[float, float]]]:
    """Return mapping n_taxa -> sorted list of (p, branch_correctness_percent)."""
    out: Dict[int, List[Tuple[float, float]]] = {}
    for sub in sorted(sweep_dir.iterdir()):
        if not sub.is_dir() or not sub.name.startswith("n"):
            continue
        if not ((sub / "nj_meta.json").exists()
                or (sub / "nj_results.json").exists()):
            continue
        r = load_nj_results(sub)
        n = int(r["n_taxa"])
        rf_mean = np.array([float(e["rf_mean"]) for e in r["per_p"]])
        p_arr = np.array([float(e["p"]) for e in r["per_p"]])
        correctness = _branch_correctness_pct(rf_mean, n)
        pairs = sorted(zip(p_arr.tolist(), correctness.tolist()),
                       key=lambda x: x[0])
        out[n] = pairs
    if not out:
        raise SystemExit(f"No NJ results found under {sweep_dir}")
    return out


def _compute_p_stars(data: Dict[int, List[Tuple[float, float]]]):
    rows = []
    for n in sorted(data.keys()):
        pairs = data[n]
        p_vals = np.array([p for p, _ in pairs])
        correctness = np.array([s for _, s in pairs])
        thr_p_star = {}
        for t in P_STAR_THRESHOLDS:
            p = find_plateau_threshold(p_vals, correctness, threshold=t)
            thr_p_star[t] = p if p is not None else np.nan
        rows.append({
            "n": n,
            "p_star": thr_p_star[P_STAR_THRESHOLD],
            "p_star_by_threshold": thr_p_star,
            "p_values": p_vals.tolist(),
            "correctness": correctness.tolist(),
        })
    return rows


def _plot_q_overlay(sweep_dir: Path, out_path: Path,
                    max_samples: int = 5000) -> None:
    """One panel per p; KDE of raw Q for adj (solid) and non-adj (dashed),
    one colour per n."""
    from scipy.stats import gaussian_kde
    subs = sorted([d for d in sweep_dir.iterdir()
                   if d.is_dir() and ((d / "nj_meta.json").exists()
                                      or (d / "nj_results.json").exists())])
    per_n = {load_nj_results(s)["n_taxa"]: load_nj_results(s) for s in subs}
    n_vals = sorted(per_n.keys())
    p_values = sorted({e["p"] for e in per_n[n_vals[-1]]["per_p"]})

    cmap = plt.get_cmap("viridis")
    n_p = len(p_values)
    ncols = int(np.ceil(np.sqrt(n_p)))
    nrows = int(np.ceil(n_p / ncols))
    fig, axes_grid = plt.subplots(nrows, ncols,
                                  figsize=(3.2 * ncols, 3.0 * nrows),
                                  squeeze=False)
    axes = axes_grid.flatten()
    for extra_ax in axes[n_p:]:
        extra_ax.set_visible(False)
    rng = np.random.default_rng(0)

    for col, p in enumerate(p_values):
        ax = axes[col]
        all_vals = []
        for n in n_vals:
            entry = next((e for e in per_n[n]["per_p"] if e["p"] == p), None)
            if entry is None:
                continue
            for arr in (entry["q_adj"], entry["q_nonadj"]):
                arr = np.asarray(arr, dtype=float)
                if arr.size > 0:
                    all_vals.append(arr)
        if not all_vals:
            continue
        cat = np.concatenate(all_vals)
        lo, hi = cat.min(), cat.max()
        if hi <= lo:
            hi = lo + 1e-9
        pad = 0.05 * (hi - lo)
        x_grid = np.linspace(lo - pad, hi + pad, 200)

        for n in n_vals:
            entry = next((e for e in per_n[n]["per_p"] if e["p"] == p), None)
            if entry is None:
                continue
            color = cmap(n_vals.index(n) / max(1, len(n_vals) - 1))
            for arr, style, tag in (
                (entry["q_adj"], "-", "adj"),
                (entry["q_nonadj"], "--", "non-adj"),
            ):
                arr = np.asarray(arr, dtype=float)
                if arr.size < 2:
                    continue
                if arr.size > max_samples:
                    arr = rng.choice(arr, size=max_samples, replace=False)
                # gaussian_kde fails on zero-variance arrays
                if np.allclose(arr, arr[0]):
                    continue
                try:
                    kde = gaussian_kde(arr)
                    ax.plot(x_grid, kde(x_grid), color=color, linestyle=style,
                            linewidth=1.2,
                            label=f"n={n} ({tag})" if col == 0 else None)
                except Exception:
                    continue
        ax.set_title(f"$p={p:.4g}$", fontsize=10)
        ax.set_xlabel(r"$Q$")
        if col == 0:
            ax.set_ylabel("density")
            ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle(r"NJ $Q$-criterion: adj (solid) vs non-adj (dashed), per $n$",
                 fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def _plot_overlay(sweep_dir: Path, out_path: Path) -> None:
    """Up to four panels — raw spec-norm + nRF on the top row; max-entry
    error and relative spec-norm on the bottom row (if available).

    The raw spec-norm panel keeps the dimensional view (||·||_2 of an n×n
    matrix grows with n, so this panel *should* spread with n at fixed p —
    it's not a bug). The two bottom-row panels show the per-entry truth:
    ``||D - D̂||_∞`` is governed by the worst sampled/imputed entry and
    should be roughly flat across n; ``||D - D̂||_2 / ||D||_2`` divides out
    the dimensional growth and should collapse the curves.

    Bottom row is rendered iff every loaded result has the extra metrics
    merged in (i.e. ``nj_meta_extra.json`` exists for each subdir — produced
    by ``scripts/nj_recompute_normalized_metrics.py``). Otherwise we fall
    back to the legacy 1×2 layout so this script still works on
    pre-diagnostic sweeps.
    """
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and ((d / "nj_meta.json").exists()
                               or (d / "nj_results.json").exists())]
    # Sort by n_taxa so the legend is monotone in n.
    subs = sorted(subs, key=lambda d: load_nj_results(d)["n_taxa"])
    loaded = [load_nj_results(s) for s in subs]
    n_vals = [r["n_taxa"] for r in loaded]
    has_extras = all(
        any("linf_mean" in e for e in r["per_p"]) for r in loaded
    )

    if has_extras:
        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        (ax_s, ax_r), (ax_linf, ax_rel) = axes
    else:
        fig, (ax_s, ax_r) = plt.subplots(1, 2, figsize=(12, 5))
        ax_linf = ax_rel = None

    cmap = plt.get_cmap("viridis")

    for r in loaded:
        n = r["n_taxa"]
        color = cmap(n_vals.index(n) / max(1, len(n_vals) - 1))
        per_p = sorted(r["per_p"], key=lambda e: e["p"])
        p = np.array([e["p"] for e in per_p])
        sn_mean = np.array([e["spec_norm_mean"] for e in per_p])
        sn_std = np.array([e["spec_norm_std"] for e in per_p])
        rf_mean = np.array([e["rf_mean"] for e in per_p])
        rf_std = np.array([e["rf_std"] for e in per_p])
        max_rf = max(1, 2 * (n - 3))
        nrf_mean = rf_mean / max_rf
        nrf_std = rf_std / max_rf

        ax_s.errorbar(p, sn_mean, yerr=sn_std, marker="o", capsize=3,
                      color=color, label=f"n={n}")
        ax_r.errorbar(p, nrf_mean, yerr=nrf_std, marker="s", capsize=3,
                      color=color, label=f"n={n}")

        if has_extras:
            linf_mean = np.array([e.get("linf_mean", np.nan) for e in per_p])
            linf_std = np.array([e.get("linf_std", 0.0) for e in per_p])
            rel_mean = np.array([e.get("rel_specnorm_mean", np.nan) for e in per_p])
            rel_std = np.array([e.get("rel_specnorm_std", 0.0) for e in per_p])
            ax_linf.errorbar(p, linf_mean, yerr=linf_std, marker="^",
                             capsize=3, color=color, label=f"n={n}")
            ax_rel.errorbar(p, rel_mean, yerr=rel_std, marker="d", capsize=3,
                            color=color, label=f"n={n}")

    ax_s.set_xlabel("sampling fraction $p$")
    ax_s.set_ylabel(r"$\|D - \hat D\|_2$")
    ax_s.set_title("Distance spectral-norm error (raw — grows with $n$)")
    ax_s.set_yscale("log")
    ax_s.grid(True, which="both", alpha=0.3)
    ax_s.legend(fontsize=8)

    ax_r.set_xlabel("sampling fraction $p$")
    ax_r.set_ylabel(r"normalised RF $= \mathrm{RF}/(2(n-3))$")
    ax_r.set_title("Normalised RF vs $p$ (classical NJ)")
    ax_r.set_ylim(-0.02, 1.02)
    ax_r.grid(True, alpha=0.3)
    ax_r.legend(fontsize=8)

    if has_extras:
        ax_linf.set_xlabel("sampling fraction $p$")
        ax_linf.set_ylabel(r"$\|D - \hat D\|_\infty$")
        ax_linf.set_title(r"Max-entry error (per-entry — roughly flat in $n$)")
        ax_linf.set_yscale("log")
        ax_linf.grid(True, which="both", alpha=0.3)
        ax_linf.legend(fontsize=8)

        ax_rel.set_xlabel("sampling fraction $p$")
        ax_rel.set_ylabel(r"$\|D - \hat D\|_2 \,/\, \|D\|_2$")
        ax_rel.set_title(r"Relative spectral-norm error (curves should collapse)")
        ax_rel.set_yscale("log")
        ax_rel.grid(True, which="both", alpha=0.3)
        ax_rel.legend(fontsize=8)

    tree_models = sorted({str(r.get("tree_model", "?")) for r in loaded})
    tree_label = ", ".join(tree_models)
    fig.suptitle(
        f"Classical NJ subsampling: overlay across $n$  —  tree model: {tree_label}",
        fontsize=11,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def _plot(rows: List[Dict], out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    ax0 = axes[0]
    cmap = plt.get_cmap("viridis")
    n_vals_sorted = sorted({r["n"] for r in rows})
    eps = 1e-12
    for r in rows:
        color = cmap(n_vals_sorted.index(r["n"]) / max(1, len(n_vals_sorted) - 1))
        p_arr = np.asarray(r["p_values"])
        one_minus_p = np.maximum(1.0 - p_arr, eps)
        ax0.plot(one_minus_p, r["correctness"], marker="o",
                 color=color, label=f"n={r['n']}")
    ax0.axhline(P_STAR_THRESHOLD, color="grey", linestyle=":", linewidth=1,
                label=f"{P_STAR_THRESHOLD:.0f}% level")
    ax0.set_xscale("log")
    ax0.set_xlabel(r"$1 - p$ (fraction dropped)")
    ax0.set_ylabel(r"branch correctness $100\,(1 - \mathrm{nRF})$ (%)")
    ax0.set_title("Branch correctness vs $1-p$, per $n$ (NJ)")
    ax0.set_ylim(-5, 105)
    ax0.grid(True, alpha=0.3)
    ax0.legend(fontsize=8, ncol=2)

    ax1 = axes[1]
    n_arr = np.array([r["n"] for r in rows], dtype=float)
    thr_cmap = plt.get_cmap("plasma")
    for t_idx, t in enumerate(P_STAR_THRESHOLDS):
        color = thr_cmap(t_idx / max(1, len(P_STAR_THRESHOLDS) - 1))
        p_star = np.array([r["p_star_by_threshold"][t] for r in rows],
                          dtype=float)
        one_minus_pstar = 1.0 - p_star
        valid = np.isfinite(one_minus_pstar) & (one_minus_pstar > 0)
        if valid.sum() < 1:
            continue
        ax1.scatter(n_arr[valid], one_minus_pstar[valid], marker="o", s=70,
                    color=color)
        if valid.sum() >= 2:
            alpha, A, _ = fit_power_law(n_arr[valid], one_minus_pstar[valid])
            n_fit = np.logspace(np.log10(n_arr.min()), np.log10(n_arr.max()), 50)
            ax1.plot(n_fit, evaluate_power_law(n_fit, alpha, A),
                     linestyle="--", color=color,
                     label=f"$\\geq{t:.0f}$%: $1-p^* \\propto n^{{{alpha:.2f}}}$")
        else:
            ax1.plot([], [], color=color, label=f"$\\geq{t:.0f}$% (only 1 pt)")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("$n$ (taxa)")
    ax1.set_ylabel(r"$1 - p^*$")
    ax1.set_title(r"NJ critical $1-p^*$ vs $n$ (per correctness threshold)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8, loc="best")

    fig.suptitle("Classical NJ subsampling phase transition", fontsize=12)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_nj_p_star_vs_n.py <sweep_dir>", file=sys.stderr)
        sys.exit(1)
    sweep_dir = Path(sys.argv[1]).resolve()
    data = _load_sweep(sweep_dir)
    rows = _compute_p_stars(data)

    header = f"{'n':>5}  " + "  ".join(f"≥{t:>3.0f}%" for t in P_STAR_THRESHOLDS)
    print(header)
    for r in rows:
        cells = [f"{r['n']:>5}"]
        for t in P_STAR_THRESHOLDS:
            p = r["p_star_by_threshold"].get(t)
            cells.append(f"{p:>5.3f}" if np.isfinite(p) else "  n/a")
        print("  ".join(cells))

    out_path = sweep_dir / "nj_p_star_vs_n.png"
    _plot(rows, out_path)
    _plot_overlay(sweep_dir, sweep_dir / "nj_overlay_specnorm_rf.png")
    _plot_q_overlay(sweep_dir, sweep_dir / "nj_overlay_q.png")

    summary = [{
        "n": r["n"],
        "p_star": r["p_star"],
        "p_star_by_threshold": {f"{t:g}": r["p_star_by_threshold"][t]
                                 for t in P_STAR_THRESHOLDS},
    } for r in rows]
    with open(sweep_dir / "nj_p_star_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"wrote {sweep_dir / 'nj_p_star_summary.json'}")


if __name__ == "__main__":
    main()
