"""Human-readable exports and run logs beside the real-cohort ``.npz`` results.

``results/real_data/<cohort>/`` is the one directory a collaborator has to find, so it
carries both machine-readable and readable forms of the same numbers:

    screen.npz        one row per tree: eta and validity per operator
    screen.csv        the same, openable in anything
    sweep/<tree>.npz  per tree: every metric, both arms, one value per p
    sweep_summary.csv median and std across trees, per p, per metric, per arm
    summary.json      headline numbers: screening verdicts, median curve per p
    config.json       cohort, gate, p-grid, reps, git commit, timestamp
    recovery_curve.png  median NMI vs p, both arms, both references
    screen.log        what the run printed
    sweep.log

Nothing here is required to *run* the experiments -- it exists so results can be read,
mailed and committed without a Python session.
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Sequence

import numpy as np

from .real_cohorts import (_slug, cohort_results_dir, screen_cache_path,
                           sweep_cache_dir)
from .real_recovery_sweep import METRICS


class Tee:
    """Duplicate stdout into a log file, so an interactive run leaves a record too."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(path, "a", buffering=1)
        self._out = sys.stdout

    def write(self, s):                      # noqa: D102
        self._out.write(s)
        self._fh.write(s)
        return len(s)

    def flush(self):                         # noqa: D102
        self._out.flush()
        self._fh.flush()

    def close(self):                         # noqa: D102
        try:
            self._fh.close()
        finally:
            if sys.stdout is self:
                sys.stdout = self._out


def export_screen_csv(cohort_name: str) -> Path | None:
    """Write ``screen.csv`` next to ``screen.npz``. Returns the path, or None if absent."""
    npz = screen_cache_path(cohort_name)
    if not npz.exists():
        return None
    rows = [r for r in np.load(npz, allow_pickle=True)["rows"] if "error" not in r]
    out = npz.with_suffix(".csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["tree", "m", "eta_L", "valid_L", "eta_B", "valid_B"])
        for r in sorted(rows, key=lambda r: r["tree"]):
            w.writerow([r["tree"], r["m"], f"{r['eta_S']:.4f}", int(r["valid_S"]),
                        f"{r['eta_B']:.4f}", int(r["valid_B"])])
    return out


def export_sweep_summary_csv(cohort_name: str, tree_ids: Sequence[str] | None = None,
                             prefix: str = "") -> Path | None:
    """Write ``sweep_summary.csv``: per p, the median and std across trees, every metric."""
    sweep_dir = sweep_cache_dir(cohort_name, prefix)
    files = sorted(sweep_dir.glob("*.npz")) if sweep_dir.is_dir() else []
    if tree_ids is not None:
        keep = set(tree_ids)
        files = [f for f in files if f.stem in keep]
    if not files:
        return None

    # Trees swept under different grids cannot share one table. Group by grid and
    # summarise the one most trees have, naming the others in the header comment.
    by_grid: dict = {}
    for f in files:
        z = np.load(f, allow_pickle=True)
        key = tuple(np.round(np.asarray(z["p_values"], float), 6))
        by_grid.setdefault(key, []).append(f)
    grid = max(by_grid, key=lambda k: len(by_grid[k]))
    chosen, others = by_grid[grid], sum(len(v) for k, v in by_grid.items() if k != grid)

    per_metric: dict = {f"{met}_{arm}": [] for met in METRICS for arm in ("L", "B")}
    per_metric["sign_L"] = []
    for f in chosen:
        z = np.load(f, allow_pickle=True)
        for key in per_metric:
            if key in z.files:
                per_metric[key].append(np.asarray(z[key], float))

    p_values = np.asarray(grid, float)
    cols = [k for k, v in per_metric.items() if v]
    out = cohort_results_dir(cohort_name) / (
        f"sweep_summary_{prefix.replace(' ', '_').lower()}.csv" if prefix
        else "sweep_summary.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        if others:
            fh.write(f"# {others} tree(s) swept on a different p-grid are not included\n")
        w.writerow(["p", "n_trees"]
                   + [f"{c}_median" for c in cols] + [f"{c}_std" for c in cols])
        stacks = {c: np.vstack(per_metric[c]) for c in cols}
        for i, pv in enumerate(p_values):
            med = [f"{np.nanmedian(stacks[c][:, i]):.6f}" for c in cols]
            sd = [f"{np.nanstd(stacks[c][:, i]):.6f}" for c in cols]
            w.writerow([f"{pv:.6g}", stacks[cols[0]].shape[0]] + med + sd)
    return out


def export_all(cohort_name: str, prefix: str = "") -> List[Path]:
    """Every readable artefact for one cohort; returns the files actually written."""
    return [p for p in (export_screen_csv(cohort_name),
                        export_sweep_summary_csv(cohort_name, prefix=prefix),
                        write_summary_json(cohort_name, prefix),
                        plot_recovery_curve(cohort_name, prefix))
            if p is not None]


def write_config_json(cohort_name: str, cfg: dict, prefix: str = "") -> Path:
    """Record what produced a run, the way ``results/runs/*/config.json`` does."""
    import json
    import subprocess

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        commit = ""
    out = cohort_results_dir(cohort_name) / (
        f"config_{_slug(prefix)}.json" if prefix else "config.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(cfg, cohort=cohort_name, prefix=prefix, commit=commit,
                   written=datetime.now().isoformat(timespec="seconds"))
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n")
    return out


def write_summary_json(cohort_name: str, prefix: str = "") -> Path | None:
    """Headline numbers: screening verdicts, and the median curve per p for both arms."""
    import json

    npz = screen_cache_path(cohort_name)
    summary: dict = {"cohort": cohort_name}
    if npz.exists():
        rows = [r for r in np.load(npz, allow_pickle=True)["rows"] if "error" not in r]
        summary["screening"] = dict(
            trees=len(rows),
            m=int(rows[0]["m"]) if rows else 0,
            valid_L=sum(bool(r["valid_S"]) for r in rows),
            valid_B=sum(bool(r["valid_B"]) for r in rows),
            valid_both=sum(bool(r["valid_S"] and r["valid_B"]) for r in rows),
            median_eta_L=float(np.median([r["eta_S"] for r in rows])) if rows else None,
            median_eta_B=float(np.median([r["eta_B"] for r in rows])) if rows else None,
        )

    sweep_dir = sweep_cache_dir(cohort_name, prefix)
    files = sorted(sweep_dir.glob("*.npz")) if sweep_dir.is_dir() else []
    if files:
        stacks, p_values = {}, None
        for f in files:
            z = np.load(f, allow_pickle=True)
            pv = np.asarray(z["p_values"], float)
            if p_values is None:
                p_values = pv
            elif pv.shape != p_values.shape or not np.allclose(pv, p_values):
                continue
            for k in z.files:
                if k not in ("meta", "p_values"):
                    stacks.setdefault(k, []).append(np.asarray(z[k], float))
        summary["sweep"] = dict(
            trees=len(files), p_values=[float(x) for x in p_values],
            median={k: [float(x) for x in np.nanmedian(np.vstack(v), axis=0)]
                    for k, v in stacks.items()})

    out = cohort_results_dir(cohort_name) / (
        f"summary_{_slug(prefix)}.json" if prefix else "summary.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    return out


def plot_recovery_curve(cohort_name: str, prefix: str = "") -> Path | None:
    """``recovery_curve.png``: median NMI vs p, both arms, both references."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sweep_dir = sweep_cache_dir(cohort_name, prefix)
    files = sorted(sweep_dir.glob("*.npz")) if sweep_dir.is_dir() else []
    if not files:
        return None
    stacks, p_values = {}, None
    for f in files:
        z = np.load(f, allow_pickle=True)
        pv = np.asarray(z["p_values"], float)
        if p_values is None:
            p_values = pv
        elif pv.shape != p_values.shape or not np.allclose(pv, p_values):
            continue
        for k in ("nmi_L", "nmi_B", "nmi_gt_L", "nmi_gt_B"):
            if k in z.files:
                stacks.setdefault(k, []).append(np.asarray(z[k], float))

    fig, ax = plt.subplots(figsize=(7, 7))
    style = {"nmi_L": ("#4b5563", "-", r"$L(S)$ vs full matrix"),
             "nmi_B": ("#065f46", "-", r"$B=H\mathcal{D}H$ vs full matrix"),
             "nmi_gt_L": ("#4b5563", "--", r"$L(S)$ vs true tree"),
             "nmi_gt_B": ("#065f46", "--", r"$B=H\mathcal{D}H$ vs true tree")}
    n_trees = 0
    for key, (color, ls, label) in style.items():
        if key not in stacks:
            continue
        arr = np.vstack(stacks[key])
        n_trees = arr.shape[0]
        med, sd = np.nanmedian(arr, 0), np.nanstd(arr, 0)
        ax.plot(p_values, med, color=color, ls=ls, lw=2, marker="o", ms=3, label=label)
        if ls == "-":
            ax.fill_between(p_values, med - sd, med + sd, color=color, alpha=0.12)
    ax.set_xscale("log")
    ax.set_xlabel(r"sub-sampling fraction $p$ (log)", fontsize=12)
    ax.set_ylabel("NMI", fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=10, loc="upper left")
    ax.set_title(f"{cohort_name}: recovery over {n_trees} tree(s)"
                 + (f" [{prefix}]" if prefix else ""), fontsize=12)
    fig.tight_layout()
    out = cohort_results_dir(cohort_name) / (
        f"recovery_curve_{_slug(prefix)}.png" if prefix else "recovery_curve.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
