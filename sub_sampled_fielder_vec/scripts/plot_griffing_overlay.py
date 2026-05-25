"""Griffing first-partition recovery — paper-style overlay (3 panels).

Matches `02_real_data_sweeps/kingman_threshold_vs_theory.ipynb` style.

  * Left:   ``partition_agreement_M`` (%) vs $p$, **log-x**, one curve
            per $n$ (viridis), red dashed 95% threshold. Per-$n$ legend
            entry includes empirical $\\hat p^*$ for agreement.
  * Middle: **ARI** vs $p$, log-x, red dashed 95% threshold. Same
            viridis colours; legend includes empirical $\\hat p^*$ for
            ARI (the more demanding metric).
  * Right:  $\\hat p^*(n)$ vs $n$ on **log-log** with two curves
            (agreement vs ARI thresholds) so the n-scaling can be
            compared head-to-head.

Output: PNG + p̂*-summary JSON land under
``results/notebooks/05_nj_distance/griffing_distance_partition/`` per
the unified ``cache_io.notebook_dir`` convention.

Usage:
    python plot_griffing_overlay.py <sweep_dir> [<sweep_dir> ...]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.utils.griffing_io import load_griffing_results
from src.cache_io import notebook_dir


NOTEBOOK_REL = "05_nj_distance/griffing_distance_partition"

# Plateau thresholds. Agreement is in [0.5, 1.0] (sign-agnostic, so 0.5 == random);
# ARI is in [-1, 1] with 0 == random. We use 0.95 for both, which is a more
# stringent criterion for ARI than for agreement.
AGREEMENT_THRESHOLD = 0.95
ARI_THRESHOLD = 0.95
RECOVERY_THRESHOLDS = [0.95, 0.9, 0.8]   # used by _compute_p_stars table


def find_plateau_threshold(p_vals: np.ndarray, metric: np.ndarray,
                           threshold: float) -> float | None:
    """Smallest p such that metric(p') >= threshold for all p' >= p.

    Guards against IPW low-p spikes where a degenerate eigvec can
    coincidentally match the reference.
    """
    order = np.argsort(p_vals)[::-1]
    p_sorted = p_vals[order]
    m_sorted = metric[order]
    last_good = None
    for p, m in zip(p_sorted, m_sorted):
        if m >= threshold:
            last_good = p
        else:
            break
    return last_good


def _series_for(r: Dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-p arrays (p, agreement_frac, agreement_std, ari_mean, ari_std).

    Backward-compat: if ``ari_mean`` isn't present in a meta JSON (older
    runs), the ARI arrays come back as NaN so the plot degrades gracefully.
    """
    per_p = sorted(r["per_p"], key=lambda e: e["p"])
    p = np.array([e["p"] for e in per_p])
    agreement = np.array([e["recovery_mean"] for e in per_p])
    agreement_std = np.array([e.get("recovery_std", 0.0) for e in per_p])
    ari_mean = np.array([e.get("ari_mean", np.nan) for e in per_p])
    ari_std = np.array([e.get("ari_std", 0.0) for e in per_p])
    return p, agreement, agreement_std, ari_mean, ari_std


def _plot_overlay(sweep_dir: Path, out_path: Path) -> None:
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "griffing_meta.json").exists()]
    subs = sorted(subs, key=lambda d: load_griffing_results(d)["n_taxa"])
    loaded = [load_griffing_results(s) for s in subs]
    n_vals = [r["n_taxa"] for r in loaded]
    tree_models = sorted({str(r.get("tree_model", "?")) for r in loaded})
    tree_label = ", ".join(tree_models)
    seq_len = next((r.get("seq_len") for r in loaded if r.get("seq_len")), "?")
    imp_set = sorted({str(r.get("imputation", "zero")) for r in loaded})
    imp_label = ", ".join(imp_set)
    has_ari = all(not np.isnan(_series_for(r)[3]).all() for r in loaded)

    fig, (ax_l, ax_m, ax_r) = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(
        f"Griffing first partition — tree={tree_label}, L={seq_len}, "
        f"sampling=uniform IPW (imputation={imp_label})"
    )

    cmap = plt.get_cmap("viridis")
    n_to_color = {n: cmap(i / max(1, len(n_vals) - 1)) for i, n in enumerate(n_vals)}

    pstar_agreement: Dict[int, float | None] = {}
    pstar_ari: Dict[int, float | None] = {}
    for r in loaded:
        n = r["n_taxa"]
        p, agreement, agreement_std, ari_mean, ari_std = _series_for(r)
        color = n_to_color[n]

        # ---- Left: partition agreement (%) ----
        ps_a = find_plateau_threshold(p, agreement, AGREEMENT_THRESHOLD)
        pstar_agreement[n] = ps_a
        label_a = (f"{ps_a:.4f}" if ps_a is not None and np.isfinite(ps_a) else "n/a")
        ax_l.errorbar(p, 100.0 * agreement, yerr=100.0 * agreement_std,
                      marker="o", capsize=3, color=color, ms=4, lw=1,
                      label=fr"n={n}, $\hat p^*$={label_a}")

        # ---- Middle: ARI ----
        if has_ari:
            ps_b = find_plateau_threshold(p, ari_mean, ARI_THRESHOLD)
            pstar_ari[n] = ps_b
            label_b = (f"{ps_b:.4f}" if ps_b is not None and np.isfinite(ps_b) else "n/a")
            ax_m.errorbar(p, ari_mean, yerr=ari_std,
                          marker="s", capsize=3, color=color, ms=4, lw=1,
                          label=fr"n={n}, $\hat p^*$={label_b}")

    # Left axis cosmetics
    ax_l.axhline(100.0 * AGREEMENT_THRESHOLD, color="red", ls="--", lw=1,
                 label=f"{int(AGREEMENT_THRESHOLD * 100)}% threshold")
    ax_l.set_xscale("log")
    ax_l.set_xlabel(r"Sampling probability $p$")
    ax_l.set_ylabel("partition_agreement_M (%)")
    ax_l.set_title(fr"Partition agreement — {tree_label}")
    ax_l.set_ylim(48, 102)
    ax_l.legend(fontsize=8, loc="lower right")
    ax_l.grid(True, which="both", alpha=0.3)

    # Middle axis cosmetics
    if has_ari:
        ax_m.axhline(ARI_THRESHOLD, color="red", ls="--", lw=1,
                     label=f"{int(ARI_THRESHOLD * 100)}% threshold")
        ax_m.set_xscale("log")
        ax_m.set_xlabel(r"Sampling probability $p$")
        ax_m.set_ylabel("ARI")
        ax_m.set_title(fr"Adjusted Rand Index — {tree_label}")
        ax_m.set_ylim(-0.04, 1.04)
        ax_m.legend(fontsize=8, loc="lower right")
        ax_m.grid(True, which="both", alpha=0.3)
    else:
        ax_m.text(0.5, 0.5, "ari_mean not present in meta\n(re-run sweep)",
                  ha="center", va="center", transform=ax_m.transAxes)
        ax_m.set_axis_off()

    # ---------- Right: p̂*(n) scaling, two curves ----------
    def _xy(pstar_dict):
        ns = [n for n in n_vals if pstar_dict.get(n) is not None
              and np.isfinite(pstar_dict[n])]
        return np.array(ns, dtype=float), np.array([pstar_dict[int(n)] for n in ns],
                                                    dtype=float)

    na, pa = _xy(pstar_agreement)
    if na.size > 0:
        ax_r.plot(na, pa, "-o",
                  label=fr"agreement $\geq$ {int(AGREEMENT_THRESHOLD*100)}%")
    if has_ari:
        nb, pb = _xy(pstar_ari)
        if nb.size > 0:
            ax_r.plot(nb, pb, "-s",
                      label=fr"ARI $\geq$ {int(ARI_THRESHOLD*100)}%")
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
    """Per (n, threshold): smallest p with mean partition_agreement >= threshold
    AND every higher p also at that level (plateau definition).

    Backward-compatible — only uses ``recovery_mean`` (always present).
    """
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
        ari_arr = np.array([e.get("ari_mean", np.nan) for e in per_p])
        thr_p_star = {
            t: find_plateau_threshold(p_arr, rec_arr, t)
            for t in RECOVERY_THRESHOLDS
        }
        ari_p_star = {
            t: find_plateau_threshold(p_arr, ari_arr, t)
            for t in RECOVERY_THRESHOLDS
        } if not np.isnan(ari_arr).all() else None
        rows.append({
            "n": n,
            "p_star_by_threshold": thr_p_star,
            "ari_p_star_by_threshold": ari_p_star,
        })
    return rows


def _tree_model_for(sweep_dir: Path) -> str:
    subs = [d for d in sweep_dir.iterdir()
            if d.is_dir() and (d / "griffing_meta.json").exists()]
    if not subs:
        return "unknown"
    meta = load_griffing_results(subs[0])
    return str(meta.get("tree_model", "unknown"))


def _process_one(sweep_dir: Path) -> None:
    rows = _compute_p_stars(sweep_dir)
    print(f"\n=== {sweep_dir.name} ===")
    print(f"{'n':>5}  " + "  ".join(f"agr>={int(t*100)}%" for t in RECOVERY_THRESHOLDS)
          + "  | " + "  ".join(f"ari>={int(t*100)}%" for t in RECOVERY_THRESHOLDS))
    for r in rows:
        cells = [f"{r['n']:>5}"]
        for t in RECOVERY_THRESHOLDS:
            p = r["p_star_by_threshold"].get(t)
            cells.append(f"{p:>7.3f}" if (p is not None and np.isfinite(p)) else "    n/a")
        cells.append(" |")
        if r["ari_p_star_by_threshold"] is not None:
            for t in RECOVERY_THRESHOLDS:
                p = r["ari_p_star_by_threshold"].get(t)
                cells.append(f"{p:>7.3f}" if (p is not None and np.isfinite(p)) else "    n/a")
        else:
            cells.extend(["    n/a"] * len(RECOVERY_THRESHOLDS))
        print("  ".join(cells))

    tree_model = _tree_model_for(sweep_dir)
    out_dir = notebook_dir(NOTEBOOK_REL)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"griffing_overlay_{tree_model}.png"
    _plot_overlay(sweep_dir, png_path)

    summary_path = out_dir / f"griffing_p_star_summary_{tree_model}.json"
    summary = []
    for r in rows:
        entry = {
            "n": r["n"],
            "p_star_agreement": {f"{int(t*100):d}": r["p_star_by_threshold"][t]
                                  for t in RECOVERY_THRESHOLDS},
        }
        if r["ari_p_star_by_threshold"] is not None:
            entry["p_star_ari"] = {f"{int(t*100):d}": r["ari_p_star_by_threshold"][t]
                                    for t in RECOVERY_THRESHOLDS}
        summary.append(entry)
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"wrote {summary_path}")


def main():
    if len(sys.argv) < 2:
        print("usage: plot_griffing_overlay.py <sweep_dir> [<sweep_dir> ...]",
              file=sys.stderr)
        sys.exit(1)
    for raw in sys.argv[1:]:
        sweep_dir = Path(raw).resolve()
        if not sweep_dir.is_dir():
            print(f"skipping {sweep_dir} (not a directory)", file=sys.stderr)
            continue
        _process_one(sweep_dir)


if __name__ == "__main__":
    main()
