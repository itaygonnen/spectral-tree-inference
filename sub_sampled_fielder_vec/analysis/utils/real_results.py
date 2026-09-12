"""What one run leaves behind: a single directory holding every cohort it covered.

    results/real_data/runs/<timestamp>-<name>/
        config.json         what was asked for: cohorts, gate, p-grid, reps, commit, time
        summary.json        headline numbers per cohort
        run.log             everything the run printed
        screening.csv       every tree of every cohort: eta and validity per operator
        per_tree.csv        every tree x every p: all metrics, both arms   (long format)
        curves.csv          per cohort x p: median, std, quartiles and a bootstrap
                            interval for the median, across trees          (plot-ready)
        recovery_curve.png  median NMI vs p, every cohort, both arms

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


def _sweep_arrays(cohort_name: str, tree_ids: Sequence[str] | None = None,
                  p_values: Sequence[float] | None = None):
    """``(p_values, {tree: {key: array}})`` for one p-grid.

    ``p_values`` is THE grid this run swept: the cache is shared across runs, so without
    it an export can pick up trees another run left behind on a different grid and plot
    them instead -- which is how a run asked for p>=0.01 produced a curve starting at
    1e-4. Only when no grid is given (a screening-only export) does it fall back to the
    grid most cached trees share.
    """
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
    if p_values is not None:
        want = tuple(np.round(np.asarray(p_values, float), 6))
        if want not in by_grid:
            return None, {}
        grid = want
    else:
        grid = max(by_grid, key=lambda k: len(by_grid[k]))
    out = {}
    for f in by_grid[grid]:
        z = np.load(f, allow_pickle=True)
        arrays = {}
        for k in z.files:
            if k in ("meta", "p_values"):
                continue
            a = z[k]
            # rule_L is a string; everything else is numeric. Keep it as it is rather
            # than forcing float on the whole archive.
            arrays[k] = a if a.dtype.kind in "OUS" else np.asarray(a, float)
        out[f.stem] = arrays
    return np.asarray(grid, float), out


def write_failures_csv(run_dir: Path, cohorts: Sequence[str]) -> Path | None:
    """Trees the screen could not process, with the error. None when there are none."""
    bad = []
    for c in cohorts:
        npz = screen_cache_path(c)
        if not npz.exists():
            continue
        bad += [(c, r) for r in np.load(npz, allow_pickle=True)["rows"] if "error" in r]
    if not bad:
        return None
    out = run_dir / "failures.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "tree", "error"])
        for c, r in bad:
            w.writerow([c, r["tree"], r["error"]])
    return out


def write_screening_csv(run_dir: Path, cohorts: Sequence[str]) -> Path | None:
    """Every tree of every cohort, one row: eta and validity per operator."""
    rows = [(c, r) for c in cohorts for r in _screen_rows(c)]
    if not rows:
        return None
    out = run_dir / "screening.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "m", "tree",
                    "eta_L", "valid_L", "rule_L",
                    "eta_L_kmeans", "valid_L_kmeans", "eta_L_sign", "valid_L_sign",
                    "eta_B", "valid_B"])
        for c, r in sorted(rows, key=lambda x: (x[0], x[1]["tree"])):
            w.writerow([c, r["m"], r["tree"],
                        f"{r['eta_S']:.4f}", int(r["valid_S"]), r.get("rule_S", ""),
                        f"{r.get('eta_S_kmeans', float('nan')):.4f}",
                        int(r.get("valid_S_kmeans", 0)),
                        f"{r.get('eta_S_sign', float('nan')):.4f}",
                        int(r.get("valid_S_sign", 0)),
                        f"{r['eta_B']:.4f}", int(r["valid_B"])])
    return out


def write_curve_csvs(run_dir: Path, cohort_ids: Dict[str, Sequence[str]],
                     p_values: Sequence[float] | None = None) -> List[Path]:
    """``per_tree.csv`` (raw) and ``curves.csv`` (median/std), all cohorts in one file."""
    written: List[Path] = []
    per_tree_rows: List[list] = []
    curve_rows: List[list] = []
    cols: List[str] = []

    for cohort, ids in cohort_ids.items():
        grid, trees = _sweep_arrays(cohort, ids, p_values)
        if grid is None:
            continue
        keys = [k for k in ALL_METRICS
                if any(k in v and v[k].dtype.kind == "f" and v[k].size == len(grid)
                       for v in trees.values())]
        cols = cols or keys
        for tree, arrays in sorted(trees.items()):
            rule = str(arrays.get("rule_L", [""])[0]) if "rule_L" in arrays else ""
            eta_ref_L = (float(arrays["eta_ref_L_" + rule][0])
                         if rule and f"eta_ref_L_{rule}" in arrays else float("nan"))
            eta_ref_B = (float(arrays["eta_ref_B"][0])
                         if "eta_ref_B" in arrays else float("nan"))
            for i, p in enumerate(grid):
                per_tree_rows.append(
                    [cohort, tree, f"{p:.6g}", rule,
                     f"{eta_ref_L:.4f}", f"{eta_ref_B:.4f}"]
                    + [f"{arrays[k][i]:.6f}" if k in arrays else "" for k in cols])
        rng = np.random.default_rng(0)      # fixed, so re-exporting a run is stable
        for i, p in enumerate(grid):
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
            w.writerow(["cohort", "tree", "p", "rule_L",
                        "eta_ref_L", "eta_ref_B"] + cols)
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
                  status: str = "completed", note: str = "",
                  p_values: Sequence[float] | None = None,
                  selection: Dict | None = None) -> Path:
    """Headline numbers per cohort, plus how the run ended.

    ``status`` matters: an interrupted run still exports every tree that finished, so its
    CSVs look exactly like a complete run's, only shorter. Without this the difference is
    invisible.
    """
    summary: dict = {"status": status, "cohorts": {}}
    if selection:
        # why trees_selected is what it is: the gate and the imbalance cap that produced it
        summary["selection"] = selection
    if note:
        summary["note"] = note
    summary["finished"] = datetime.now().isoformat(timespec="seconds")
    for cohort, ids in cohort_ids.items():
        rows = _screen_rows(cohort)
        grid, trees = _sweep_arrays(cohort, ids, p_values)
        entry: dict = {"trees_selected": len(list(ids))}
        if rows:
            entry["screening"] = dict(
                trees=len(rows), m=int(rows[0]["m"]),
                valid_L=sum(bool(r["valid_S"]) for r in rows),
                valid_B=sum(bool(r["valid_B"]) for r in rows),
                valid_both=sum(bool(r["valid_S"] and r["valid_B"]) for r in rows),
                median_eta_L=float(np.median([r["eta_S"] for r in rows])),
                median_eta_B=float(np.median([r["eta_B"] for r in rows])))
        if grid is not None and trees:
            entry["sweep"] = {"trees_requested": len(list(ids)),
                              "trees_done": len(trees),
                              "p_values": [float(x) for x in grid]}
            for key in ("nmi_L", "nmi_B", "nmi_gt_L", "nmi_gt_B"):
                arr = [v[key] for v in trees.values() if key in v]
                if arr:
                    entry["sweep"][f"median_{key}"] = [
                        float(x) for x in np.nanmedian(np.vstack(arr), 0)]
        summary["cohorts"][cohort] = entry
    out = run_dir / "summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    return out


def plot_recovery(run_dir: Path, cohort_ids: Dict[str, Sequence[str]],
                  p_values: Sequence[float] | None = None) -> Path | None:
    """One figure for the whole run: median NMI vs p, every cohort, both arms.

    Scored against each arm's full-matrix split only. The ``*_gt`` columns (vs the true
    tree's top bipartition) stay in the CSVs, but they are not plotted: that split is the
    root's two child subtrees, which on these trees is near-degenerate (998/2, 5944/56),
    so NMI against it is ~0 for every method at every p and the panel showed nothing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = []
    for cohort, ids in cohort_ids.items():
        grid, trees = _sweep_arrays(cohort, ids, p_values)
        if grid is not None and trees:
            series.append((cohort, grid, trees))
    if not series:
        return None

    # Colour distinguishes the OPERATORS -- the comparison the figure exists for -- and
    # line style distinguishes cohorts. The old figure did the opposite, so both arms of a
    # single cohort came out the same colour and only solid-vs-dashed told them apart.
    ARM = {"L": ("#4b5563", r"$L(S)$ Fiedler"),
           "B": ("#065f46", r"$B = H\mathcal{D}H$")}
    STYLES = ["-", "--", ":", "-."]
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for i, (cohort, p_values, trees) in enumerate(series):
        ls = STYLES[i % len(STYLES)]
        for arm, (color, arm_label) in ARM.items():
            key = f"nmi_{arm}"
            arr = [v[key] for v in trees.values() if key in v]
            if not arr:
                continue
            arr = np.vstack(arr)
            med = np.nanmedian(arr, 0)
            q25, q75 = (np.nanpercentile(arr, 25, axis=0),
                        np.nanpercentile(arr, 75, axis=0))
            rules = {str(v["rule_L"][0]) for v in trees.values() if "rule_L" in v}
            cut = (f", {rules.pop()} cut" if arm == "L" and len(rules) == 1
                   else (", per-tree cut" if arm == "L" and rules else ", sign cut"))
            ax.plot(p_values, med, color=color, ls=ls, lw=2, marker="o", ms=3,
                    label=f"{arm_label}{cut} - {cohort}, {arr.shape[0]} trees")
            # the middle half of the trees; a +-std band would leave [0, 1]
            ax.fill_between(p_values, q25, q75, color=color, alpha=0.13)
    ax.set_xscale("log")
    ax.set_xlabel(r"sub-sampling fraction $p$ (log)", fontsize=12)
    ax.set_ylabel("median NMI vs that operator's own full-matrix split", fontsize=11)
    ax.set_title("colour = operator, line style = cohort; band = middle half of trees",
                 fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.95)
    fig.tight_layout()
    out = run_dir / "recovery_curve.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def export_run(run_dir: Path, cohort_ids: Dict[str, Sequence[str]],
               status: str = "completed", note: str = "",
               p_values: Sequence[float] | None = None,
               selection: Dict | None = None) -> List[Path]:
    """Everything readable for a finished run. Returns the files written.

    ``p_values`` pins the export to the grid THIS run swept; without it the shared cache
    can contribute another run's trees.
    """
    written = [p for p in (write_screening_csv(run_dir, list(cohort_ids)),
                           write_failures_csv(run_dir, list(cohort_ids))) if p]
    written += write_curve_csvs(run_dir, cohort_ids, p_values)
    written.append(write_summary(run_dir, cohort_ids, status, note, p_values, selection))
    plot = plot_recovery(run_dir, cohort_ids, p_values)
    if plot:
        written.append(plot)
    return written
