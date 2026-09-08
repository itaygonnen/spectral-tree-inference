"""What one run leaves behind: a single directory holding every cohort it covered.

    results/real_data/runs/<timestamp>-<name>/
        config.json         what was asked for: cohorts, gate, p-grid, reps, commit, time
        summary.json        headline numbers per cohort
        run.log             everything the run printed
        screening.csv       every tree of every cohort: eta and validity per operator
        per_tree.csv        every tree x every p: all metrics, both arms   (long format)
        curves.csv          per cohort x p: median, std, quartiles and a bootstrap
                            interval for the median, across trees          (plot-ready)
        recovery_curve.png  median NMI vs p, every cohort, both arms, both references

``curves.csv`` is the file to plot from: one row per (cohort, p), one column per
metric x arm x {median, std}. ``per_tree.csv`` is the same numbers before aggregation,
for a different cut (per-tree spread, tree-by-tree outliers).

The per-tree ``.npz`` files stay in ``_cache/<cohort>/`` so a killed run resumes; a run
directory is an export of them and never changes after the run.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from src.plots.plot_operator_threshold import _median_ci

from .real_cohorts import screen_cache_path, sweep_cache_dir
from .real_recovery_sweep import GT_METRICS, METRICS

ARMS = ("L", "B")
ALL_METRICS = tuple(f"{m}_{a}" for m in METRICS + GT_METRICS for a in ARMS) + ("sign_L",)

# Spread across TREES, per p. std alone is a poor summary here: NMI is bounded and often
# bimodal (a tree either recovers its split or does not), so the mean +- std band leaves
# the data and hides the split. Quartiles say where the middle half of the trees sit;
# lo95/hi95 are a bootstrap interval for the median itself, i.e. how well this cohort
# pins the curve down, and they narrow as trees are added while the quartiles do not.
STATS = ("median", "std", "q25", "q75", "lo95", "hi95", "n")


def _stats(vals: np.ndarray, rng) -> Dict[str, float]:
    vals = np.asarray(vals, float)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return {k: float("nan") for k in STATS}
    med, lo, hi = _median_ci(vals, rng)
    return dict(median=med, std=float(np.std(vals)),
                q25=float(np.percentile(vals, 25)), q75=float(np.percentile(vals, 75)),
                lo95=lo, hi95=hi, n=float(vals.size))


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


def _commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return ""


def write_config(run_dir: Path, cfg: dict, cohorts: Sequence[str], stage: str) -> Path:
    """What was asked for, written before the work starts."""
    out = run_dir / "config.json"
    out.write_text(json.dumps(dict(cfg, cohorts=list(cohorts), stage=stage,
                                   commit=_commit(),
                                   started=datetime.now().isoformat(timespec="seconds")),
                              indent=2, default=str) + "\n")
    return out


def _screen_rows(cohort_name: str) -> List[dict]:
    npz = screen_cache_path(cohort_name)
    if not npz.exists():
        return []
    return [r for r in np.load(npz, allow_pickle=True)["rows"] if "error" not in r]


def _sweep_arrays(cohort_name: str, tree_ids: Sequence[str] | None = None):
    """``(p_values, {tree: {key: array}})`` for the p-grid most trees share."""
    d = sweep_cache_dir(cohort_name)
    files = sorted(d.glob("*.npz")) if d.is_dir() else []
    if tree_ids is not None:
        keep = set(tree_ids)
        files = [f for f in files if f.stem in keep]
    if not files:
        return None, {}
    by_grid: Dict[tuple, list] = {}
    for f in files:
        z = np.load(f, allow_pickle=True)
        by_grid.setdefault(tuple(np.round(np.asarray(z["p_values"], float), 6)),
                           []).append(f)
    grid = max(by_grid, key=lambda k: len(by_grid[k]))
    out = {}
    for f in by_grid[grid]:
        z = np.load(f, allow_pickle=True)
        out[f.stem] = {k: np.asarray(z[k], float) for k in z.files
                       if k not in ("meta", "p_values")}
    return np.asarray(grid, float), out


def write_screening_csv(run_dir: Path, cohorts: Sequence[str]) -> Path | None:
    """Every tree of every cohort, one row: eta and validity per operator."""
    rows = [(c, r) for c in cohorts for r in _screen_rows(c)]
    if not rows:
        return None
    out = run_dir / "screening.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "m", "tree", "eta_L", "valid_L", "eta_B", "valid_B"])
        for c, r in sorted(rows, key=lambda x: (x[0], x[1]["tree"])):
            w.writerow([c, r["m"], r["tree"], f"{r['eta_S']:.4f}", int(r["valid_S"]),
                        f"{r['eta_B']:.4f}", int(r["valid_B"])])
    return out


def write_curve_csvs(run_dir: Path, cohort_ids: Dict[str, Sequence[str]]) -> List[Path]:
    """``per_tree.csv`` (raw) and ``curves.csv`` (median/std), all cohorts in one file."""
    written: List[Path] = []
    per_tree_rows: List[list] = []
    curve_rows: List[list] = []
    cols: List[str] = []

    for cohort, ids in cohort_ids.items():
        p_values, trees = _sweep_arrays(cohort, ids)
        if p_values is None:
            continue
        keys = [k for k in ALL_METRICS if any(k in v for v in trees.values())]
        cols = cols or keys
        for tree, arrays in sorted(trees.items()):
            for i, p in enumerate(p_values):
                per_tree_rows.append(
                    [cohort, tree, f"{p:.6g}"]
                    + [f"{arrays[k][i]:.6f}" if k in arrays else "" for k in cols])
        rng = np.random.default_rng(0)      # fixed, so re-exporting a run is stable
        for i, p in enumerate(p_values):
            row: list = [cohort, f"{p:.6g}", len(trees)]
            for k in cols:
                vals = np.array([v[k][i] for v in trees.values() if k in v], float)
                st = _stats(vals, rng)
                row += [f"{st[stat]:.6f}" for stat in STATS]
            curve_rows.append(row)

    if per_tree_rows:
        out = run_dir / "per_tree.csv"
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cohort", "tree", "p"] + cols)
            w.writerows(per_tree_rows)
        written.append(out)
    if curve_rows:
        out = run_dir / "curves.csv"
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["cohort", "p", "n_trees"]
                       + [f"{k}_{stat}" for k in cols for stat in STATS])
            w.writerows(curve_rows)
        written.append(out)
    return written


def write_summary(run_dir: Path, cohort_ids: Dict[str, Sequence[str]],
                  status: str = "completed", note: str = "") -> Path:
    """Headline numbers per cohort, plus how the run ended.

    ``status`` matters: an interrupted run still exports every tree that finished, so its
    CSVs look exactly like a complete run's, only shorter. Without this the difference is
    invisible.
    """
    summary: dict = {"status": status, "cohorts": {}}
    if note:
        summary["note"] = note
    summary["finished"] = datetime.now().isoformat(timespec="seconds")
    for cohort, ids in cohort_ids.items():
        rows = _screen_rows(cohort)
        p_values, trees = _sweep_arrays(cohort, ids)
        entry: dict = {"trees_selected": len(list(ids))}
        if rows:
            entry["screening"] = dict(
                trees=len(rows), m=int(rows[0]["m"]),
                valid_L=sum(bool(r["valid_S"]) for r in rows),
                valid_B=sum(bool(r["valid_B"]) for r in rows),
                valid_both=sum(bool(r["valid_S"] and r["valid_B"]) for r in rows),
                median_eta_L=float(np.median([r["eta_S"] for r in rows])),
                median_eta_B=float(np.median([r["eta_B"] for r in rows])))
        if p_values is not None and trees:
            entry["sweep"] = {"trees_requested": len(list(ids)),
                              "trees_done": len(trees),
                              "p_values": [float(x) for x in p_values]}
            for key in ("nmi_L", "nmi_B", "nmi_gt_L", "nmi_gt_B"):
                arr = [v[key] for v in trees.values() if key in v]
                if arr:
                    entry["sweep"][f"median_{key}"] = [
                        float(x) for x in np.nanmedian(np.vstack(arr), 0)]
        summary["cohorts"][cohort] = entry
    out = run_dir / "summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    return out


def plot_recovery(run_dir: Path, cohort_ids: Dict[str, Sequence[str]]) -> Path | None:
    """One figure for the whole run: median NMI vs p, every cohort, both references."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = []
    for cohort, ids in cohort_ids.items():
        p_values, trees = _sweep_arrays(cohort, ids)
        if p_values is not None and trees:
            series.append((cohort, p_values, trees))
    if not series:
        return None

    palette = ["#065f46", "#1d4ed8", "#b45309", "#991b1b"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax, ref, ttl in ((axes[0], "", "vs full-matrix split"),
                         (axes[1], "_gt", "vs true tree's split")):
        for i, (cohort, p_values, trees) in enumerate(series):
            color = palette[i % len(palette)]
            for arm, ls in (("L", "-"), ("B", "--")):
                key = f"nmi{ref}_{arm}"
                arr = [v[key] for v in trees.values() if key in v]
                if not arr:
                    continue
                arr = np.vstack(arr)
                med = np.nanmedian(arr, 0)
                q25, q75 = (np.nanpercentile(arr, 25, axis=0),
                            np.nanpercentile(arr, 75, axis=0))
                ax.plot(p_values, med, color=color, ls=ls, lw=2, marker="o", ms=3,
                        label=f"{cohort} - {'L(S)' if arm == 'L' else 'B=HDH'} "
                              f"({arr.shape[0]} trees)")
                # the middle half of the trees; a +-std band would leave [0, 1]
                ax.fill_between(p_values, q25, q75, color=color,
                                alpha=0.16 if arm == "L" else 0.08)
        ax.set_xscale("log")
        ax.set_xlabel(r"sub-sampling fraction $p$ (log)", fontsize=12)
        ax.set_title(ttl, fontsize=12)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("median NMI across trees (band: interquartile range)",
                       fontsize=11)
    axes[0].legend(fontsize=9, loc="upper left")
    fig.tight_layout()
    out = run_dir / "recovery_curve.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def export_run(run_dir: Path, cohort_ids: Dict[str, Sequence[str]],
               status: str = "completed", note: str = "") -> List[Path]:
    """Everything readable for a finished run. Returns the files written."""
    written = [p for p in [write_screening_csv(run_dir, list(cohort_ids))] if p]
    written += write_curve_csvs(run_dir, cohort_ids)
    written.append(write_summary(run_dir, cohort_ids, status, note))
    plot = plot_recovery(run_dir, cohort_ids)
    if plot:
        written.append(plot)
    return written
