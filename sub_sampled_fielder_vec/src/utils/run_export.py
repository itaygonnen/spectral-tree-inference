"""What one run leaves behind: a single directory holding every dataset it covered.

    results/real_data/runs/<timestamp>-<name>/
        config.json         what was asked for: datasets, gate, p-grid, reps, commit, time
        summary.json        headline numbers per dataset
        run.log             everything the run printed
        screening.csv       every tree of every dataset: eta and validity per operator
        per_tree.csv        every tree x every p: all metrics, every arm   (long format)
        curves.csv          per dataset x p: median, std, quartiles and a bootstrap
                            interval for the median, across trees          (plot-ready)
        recovery_curve.png  median NMI vs p, every dataset, every arm

``curves.csv`` is the file to plot from: one row per (dataset, p), one column per
metric x arm x {median, std, ...}. ``per_tree.csv`` is the same numbers before
aggregation, for a different cut (per-tree spread, tree-by-tree outliers).

The per-tree ``.npz`` files stay in ``_cache/<dataset>/`` so a killed run resumes; a run
directory is an export of them and never changes after the run.

Nothing here knows where the trees came from -- it reads the caches a run wrote, whether
those trees were FASTA alignments or simulated.
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

from ..plots.plot_operator_threshold import _median_ci
from ..runners.operator_sweep import REF_FULL, columns
from ..runners.operators import ARM_OF, OPERATORS, resolve
from .run_paths import screen_cache_path, sweep_cache_dir

# Spread across TREES, per p. std alone is a poor summary here: NMI is bounded and often
# bimodal (a tree either recovers its split or does not), so the mean +- std band leaves
# the data and hides the split. Quartiles say where the middle half of the trees sit;
# lo95/hi95 are a bootstrap interval for the median itself, i.e. how well this dataset
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


def write_config(run_dir: Path, cfg: dict, datasets: Sequence[str], stage: str) -> Path:
    """What was asked for, written before the work starts."""
    out = run_dir / "config.json"
    out.write_text(json.dumps(dict(cfg, datasets=list(datasets), stage=stage,
                                   commit=_commit(),
                                   started=datetime.now().isoformat(timespec="seconds")),
                              indent=2, default=str) + "\n")
    return out


def _screen_rows(dataset: str) -> List[dict]:
    npz = screen_cache_path(dataset)
    if not npz.exists():
        return []
    return [r for r in np.load(npz, allow_pickle=True)["rows"] if "error" not in r]


def _sweep_arrays(dataset: str, tree_ids: Sequence[str] | None = None,
                  p_values: Sequence[float] | None = None):
    """``(p_values, {tree: {key: array}})`` for one p-grid.

    ``p_values`` is THE grid this run swept: the cache is shared across runs, so without
    it an export can pick up trees another run left behind on a different grid and plot
    them instead -- which is how a run asked for p>=0.01 produced a curve starting at
    1e-4. Only when no grid is given (a screening-only export) does it fall back to the
    grid most cached trees share.
    """
    d = sweep_cache_dir(dataset)
    files = sorted(d.glob("*.npz")) if d.is_dir() else []
    if tree_ids is not None:
        keep = set(tree_ids)
        files = [f for f in files if f.stem in keep]
    if not files:
        return None, {}
    # one pass: each archive is opened once, and its arrays are kept keyed by the grid
    # they were swept on
    by_grid: Dict[tuple, Dict[str, dict]] = {}
    for f in files:
        z = np.load(f, allow_pickle=True)
        arrays = {}
        for k in z.files:
            if k in ("meta", "p_values"):
                continue
            a = z[k]
            # rule_<arm> is a string; everything else is numeric. Keep it as it is
            # rather than forcing float on the whole archive.
            arrays[k] = a if a.dtype.kind in "OUS" else np.asarray(a, float)
        by_grid.setdefault(tuple(np.round(np.asarray(z["p_values"], float), 6)),
                           {})[f.stem] = arrays
    if p_values is not None:
        want = tuple(np.round(np.asarray(p_values, float), 6))
        if want not in by_grid:
            return None, {}
        grid = want
    else:
        grid = max(by_grid, key=lambda k: len(by_grid[k]))
    return np.asarray(grid, float), by_grid[grid]


def write_failures_csv(run_dir: Path, datasets: Sequence[str]) -> Path | None:
    """Trees the screen could not process, with the error. None when there are none."""
    bad = []
    for c in datasets:
        npz = screen_cache_path(c)
        if not npz.exists():
            continue
        bad += [(c, r) for r in np.load(npz, allow_pickle=True)["rows"] if "error" in r]
    if not bad:
        return None
    out = run_dir / "failures.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["dataset", "tree", "error"])
        for c, r in bad:
            w.writerow([c, r["tree"], r["error"]])
    return out


def write_screening_csv(run_dir: Path, datasets: Sequence[str],
                        operators: Sequence[str] = ()) -> Path | None:
    """Every tree of every dataset, one row: eta and validity per operator.

    Columns are named by ARM (``L``, ``Lsym``, ``B``) while the cache stores them by
    operator key (``S``, ``Lsym``, ``B``) -- ``eta_S`` in a table beside ``eta_B`` read
    as two different things to everyone who saw it.
    """
    rows = [(c, r) for c in datasets for r in _screen_rows(c)]
    if not rows:
        return None
    ops = resolve(operators)
    # only write an operator's columns when some row actually carries them: a screen
    # from before L_sym existed would otherwise get a wall of blanks
    present = [k for k in ops if any(f"eta_{k}" in r for _, r in rows)]
    header = ["dataset", "m", "tree"]
    for k in present:
        a = ARM_OF[k]
        header += [f"eta_{a}", f"valid_{a}", f"rule_{a}"]
        header += [f"{x}_{a}_{rule}" for rule in OPERATORS[k].rules
                   for x in ("eta", "valid")]
    out = run_dir / "screening.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for c, r in sorted(rows, key=lambda x: (x[0], x[1]["tree"])):
            line = [c, r.get("m", ""), r["tree"]]
            for k in present:
                line += [_fmt(r.get(f"eta_{k}")), _flag(r.get(f"valid_{k}")),
                         r.get(f"rule_{k}", "")]
                for rule in OPERATORS[k].rules:
                    line += [_fmt(r.get(f"eta_{k}_{rule}")),
                             _flag(r.get(f"valid_{k}_{rule}"))]
            w.writerow(line)
    return out


def _fmt(v) -> str:
    return "" if v is None else f"{float(v):.4f}"


def _flag(v) -> str:
    return "" if v is None else str(int(bool(v)))


def write_curve_csvs(run_dir: Path, dataset_ids: Dict[str, Sequence[str]],
                     p_values: Sequence[float] | None = None,
                     operators: Sequence[str] = ()) -> List[Path]:
    """``per_tree.csv`` (raw) and ``curves.csv`` (median/std), all datasets in one file."""
    written: List[Path] = []
    per_tree_rows: List[list] = []
    curve_rows: List[list] = []
    arms = [ARM_OF[k] for k in resolve(operators)]
    all_cols = columns(resolve(operators))
    cols: List[str] = []
    per_tree_head: List[str] = []

    for dataset, ids in dataset_ids.items():
        grid, trees = _sweep_arrays(dataset, ids, p_values)
        if grid is None:
            continue
        keys = [k for k in all_cols
                if any(k in v and v[k].dtype.kind == "f" and v[k].size == len(grid)
                       for v in trees.values())]
        cols = cols or keys
        if not per_tree_head:
            per_tree_head = (["dataset", "tree", "p"]
                             + [f"rule_{a}" for a in arms]
                             + [f"eta_ref_{a}" for a in arms])
        for tree, arrays in sorted(trees.items()):
            rules = [str(arrays["rule_" + a][0]) if f"rule_{a}" in arrays else ""
                     for a in arms]
            etas = [_eta_ref(arrays, a) for a in arms]
            for i, p in enumerate(grid):
                per_tree_rows.append(
                    [dataset, tree, f"{p:.6g}"] + rules + etas
                    + [f"{arrays[k][i]:.6f}" if k in arrays else "" for k in cols])
        rng = np.random.default_rng(0)      # fixed, so re-exporting a run is stable
        for i, p in enumerate(grid):
            row: list = [dataset, f"{p:.6g}", len(trees)]
            for k in cols:
                vals = np.array([v[k][i] for v in trees.values() if k in v], float)
                st = _stats(vals, rng)
                row += [f"{st[stat]:.6f}" for stat in STATS]
            curve_rows.append(row)

    if per_tree_rows:
        out = run_dir / "per_tree.csv"
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(per_tree_head + cols)
            w.writerows(per_tree_rows)
        written.append(out)
    if curve_rows:
        out = run_dir / "curves.csv"
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["dataset", "p", "n_trees"]
                       + [f"{k}_{stat}" for k in cols for stat in STATS])
            w.writerows(curve_rows)
        written.append(out)
    return written


def _eta_ref(arrays: dict, arm: str) -> str:
    """The chosen rule's reference eta, from either the new column or the old pair."""
    if f"eta_ref_{arm}" in arrays:
        return f"{float(arrays[f'eta_ref_{arm}'][0]):.4f}"
    rule = str(arrays["rule_" + arm][0]) if f"rule_{arm}" in arrays else ""
    key = f"eta_ref_{arm}_{rule}"
    return f"{float(arrays[key][0]):.4f}" if key in arrays else ""


def write_summary(run_dir: Path, dataset_ids: Dict[str, Sequence[str]],
                  status: str = "completed", note: str = "",
                  p_values: Sequence[float] | None = None,
                  selection: Dict | None = None,
                  operators: Sequence[str] = ()) -> Path:
    """Headline numbers per dataset, plus how the run ended.

    ``status`` matters: an interrupted run still exports every tree that finished, so its
    CSVs look exactly like a complete run's, only shorter. Without this the difference is
    invisible.
    """
    ops = resolve(operators)
    # the commit too: summary.json is the file that gets mailed around on its own, and
    # "which version produced this" was not answerable from it
    summary: dict = {"status": status, "commit": _commit(), "datasets": {}}
    if selection:
        # why trees_selected is what it is: the gate and the cap that produced it
        summary["selection"] = selection
    if note:
        summary["note"] = note
    summary["finished"] = datetime.now().isoformat(timespec="seconds")
    for dataset, ids in dataset_ids.items():
        rows = _screen_rows(dataset)
        grid, trees = _sweep_arrays(dataset, ids, p_values)
        entry: dict = {"trees_selected": len(list(ids))}
        if rows:
            screen = {"trees": len(rows), "m": int(rows[0].get("m", 0))}
            for k in ops:
                if not any(f"eta_{k}" in r for r in rows):
                    continue
                a = ARM_OF[k]
                screen[f"valid_{a}"] = sum(bool(r.get(f"valid_{k}")) for r in rows)
                screen[f"median_eta_{a}"] = float(np.median(
                    [r[f"eta_{k}"] for r in rows if f"eta_{k}" in r]))
            screen["valid_L_and_B"] = sum(
                bool(r.get("valid_S") and r.get("valid_B")) for r in rows)
            entry["screening"] = screen
        if grid is not None and trees:
            entry["sweep"] = {"trees_requested": len(list(ids)),
                              "trees_done": len(trees),
                              "p_values": [float(x) for x in grid]}
            for k in ops:
                for ref in (REF_FULL, "vs_truetree"):
                    key = f"nmi_{ARM_OF[k]}_{ref}"
                    arr = [v[key] for v in trees.values() if key in v]
                    if arr:
                        entry["sweep"][f"median_{key}"] = [
                            float(x) for x in np.nanmedian(np.vstack(arr), 0)]
        summary["datasets"][dataset] = entry
    out = run_dir / "summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n")
    return out


# Colour distinguishes the OPERATORS -- the comparison the figure exists for -- and line
# style distinguishes datasets. The old figure did the opposite, so both arms of a single
# dataset came out the same colour and only solid-vs-dashed told them apart.
ARM_COLOUR = {"L": "#4b5563", "Lsym": "#b45309", "B": "#065f46"}
STYLES = ["-", "--", ":", "-."]


def plot_recovery(run_dir: Path, dataset_ids: Dict[str, Sequence[str]],
                  p_values: Sequence[float] | None = None,
                  operators: Sequence[str] = ()) -> Path | None:
    """One figure for the whole run: median NMI vs p, every dataset, every arm.

    Scored against each arm's full-matrix split only. The ``vs_truetree`` columns stay in
    the CSVs, but they are not plotted: that split is the root's two child subtrees,
    which on these trees is near-degenerate (998/2, 5944/56), so NMI against it is ~0 for
    every method at every p and the panel showed nothing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = []
    for dataset, ids in dataset_ids.items():
        grid, trees = _sweep_arrays(dataset, ids, p_values)
        if grid is not None and trees:
            series.append((dataset, grid, trees))
    if not series:
        return None

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for i, (dataset, grid, trees) in enumerate(series):
        ls = STYLES[i % len(STYLES)]
        for key in resolve(operators):
            arm = ARM_OF[key]
            col = f"nmi_{arm}_{REF_FULL}"
            arr = [v[col] for v in trees.values() if col in v]
            if not arr:
                continue
            arr = np.vstack(arr)
            med = np.nanmedian(arr, 0)
            q25, q75 = (np.nanpercentile(arr, 25, axis=0),
                        np.nanpercentile(arr, 75, axis=0))
            rules = {str(v[f"rule_{arm}"][0]) for v in trees.values()
                     if f"rule_{arm}" in v}
            cut = (f", {rules.pop()} cut" if len(rules) == 1
                   else (", per-tree cut" if rules else ""))
            ax.plot(grid, med, color=ARM_COLOUR.get(arm, "#111827"), ls=ls, lw=2,
                    marker="o", ms=3,
                    label=f"{OPERATORS[key].label}{cut} - {dataset}, "
                          f"{arr.shape[0]} trees")
            # the middle half of the trees; a +-std band would leave [0, 1]
            ax.fill_between(grid, q25, q75, color=ARM_COLOUR.get(arm, "#111827"),
                            alpha=0.13)
    ax.set_xscale("log")
    ax.set_xlabel(r"sub-sampling fraction $p$ (log)", fontsize=12)
    ax.set_ylabel("median NMI vs that operator's own full-matrix split", fontsize=11)
    ax.set_title("colour = operator, line style = dataset; band = middle half of trees",
                 fontsize=10)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.95)
    fig.tight_layout()
    out = run_dir / "recovery_curve.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def export_run(run_dir: Path, dataset_ids: Dict[str, Sequence[str]],
               status: str = "completed", note: str = "",
               p_values: Sequence[float] | None = None,
               selection: Dict | None = None,
               operators: Sequence[str] = ()) -> List[Path]:
    """Everything readable for a finished run. Returns the files written.

    ``p_values`` pins the export to the grid THIS run swept; without it the shared cache
    can contribute another run's trees.
    """
    names = list(dataset_ids)
    written = [p for p in (write_screening_csv(run_dir, names, operators),
                           write_failures_csv(run_dir, names)) if p]
    written += write_curve_csvs(run_dir, dataset_ids, p_values, operators)
    written.append(write_summary(run_dir, dataset_ids, status, note, p_values,
                                 selection, operators))
    plot = plot_recovery(run_dir, dataset_ids, p_values, operators)
    if plot:
        written.append(plot)
    return written
