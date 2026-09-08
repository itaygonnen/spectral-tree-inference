"""Human-readable exports and run logs beside the real-cohort ``.npz`` results.

``results/real_data/<cohort>/`` is the one directory a collaborator has to find, so it
carries both machine-readable and readable forms of the same numbers:

    screen.npz        one row per tree: eta and validity per operator
    screen.csv        the same, openable in anything
    sweep/<tree>.npz  per tree: every metric, both arms, one value per p
    sweep_summary.csv median and std across trees, per p, per metric, per arm
    screen.log        what the run printed
    sweep.log

Nothing here is required to *run* the experiments -- it exists so results can be read,
mailed and committed without a Python session.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import List, Sequence

import numpy as np

from .real_cohorts import cohort_results_dir, screen_cache_path, sweep_cache_dir
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
    """Both CSVs for one cohort; returns the files actually written."""
    return [p for p in (export_screen_csv(cohort_name),
                        export_sweep_summary_csv(cohort_name, prefix=prefix))
            if p is not None]
