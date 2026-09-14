"""One runner behind both real-data front-ends.

``scripts/run_real_sweep.py`` (command line, nohup'd on a cluster) and
``analysis/utils/real_interactive.py`` (the menu) used to carry a copy each of tree
selection, the eta gate, the p-grid, the run directory, the log wiring, the per-dataset
error containment and the export call. The copies had drifted: the two wrote different
vocabularies into the same ``summary.json`` field and different counts into the same
selection block. Everything that is not "how the knobs were collected" lives here now, so
a menu session and a batch job with the same parameters produce the same run directory.

The front-ends supply a :class:`RunSpec` and call :func:`execute`. Nothing else.
"""
from __future__ import annotations

import sys
import traceback
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.logging import (close_log_file, set_display_mode,   # noqa: E402
                               setup_log_file)

from .real_datasets import (get_dataset, new_run_dir,              # noqa: E402
                            screen_cache_path, sweep_cache_dir)
from .real_eta_screen import run_real_eta_screen                   # noqa: E402
from .real_recovery_sweep import run_sweep                         # noqa: E402
from .real_results import Tee, export_run, write_config            # noqa: E402

# Which screened trees a sweep may use. The gate is whether that operator's partition of
# the FULL matrix is a real single-edge split of the true tree -- recovery towards a
# reference that is not a tree edge measures stability, not correctness.
#
# The keys are the machine names (they go into summary.json and onto the command line);
# the labels are what the menu shows. The screen's row keys are still ``valid_S``/``eta_S``
# for the L arm -- that schema is frozen, because renaming it would invalidate every
# cached screen -- so the lambdas read _S and only the label says L(S).
GATES: Dict[str, tuple] = {
    "both":    ("L(S) and B both cut a real tree edge",
                lambda r: bool(r.get("valid_S")) and bool(r.get("valid_B"))),
    "valid_L": ("L(S) cuts a real tree edge", lambda r: bool(r.get("valid_S"))),
    "valid_B": ("B cuts a real tree edge", lambda r: bool(r.get("valid_B"))),
    "any":     ("L(S) or B cuts a real tree edge",
                lambda r: bool(r.get("valid_S")) or bool(r.get("valid_B"))),
    "all":     ("every tree, valid or not", lambda r: True),
}
# what the CLI called the L gate before the rename; still accepted so a command the
# professor already has in his shell history keeps working
GATE_ALIASES = {"valid_S": "valid_L"}

# Reference splits this lopsided are not a recovery question. At m=1000 the L(S) k-means
# cut routinely isolates ONE taxon (1/999, eta=999): a real pendant edge, so the validity
# gate passes it, but every sub-sample then picks a different singleton and NMI sits at 0
# for every p. 20 is well above the recoverability scale (m/log m)^(1/3) ~ 5.3 at m=1000.
MAX_ETA = 20.0

# p starts at 1e-4: the transition sits near log n / n, which is 1.4e-3 at m=6000, so the
# grid has to reach below that for the curve to show a floor rather than start on the ramp.
P_MIN, P_POINTS, REPS = 1e-4, 20, 10

# Per-tree cost, measured on this machine at m=6000 (see docs/RUNBOOK.md): a screen is one
# load plus two full-matrix eigensolves; a sweep is p x reps sub-sampled Fiedler solves at
# ~9 s each, and the solve is O(m^3).
SCREEN_SECS_LARGE, SCREEN_SECS_SMALL, SOLVE_SECS_AT_6000 = 75.0, 10.0, 9.0


def gate_key(name: str) -> str:
    """Normalise a gate name: accepts a key, a deprecated alias, or a menu label."""
    if name in GATES:
        return name
    if name in GATE_ALIASES:
        return GATE_ALIASES[name]
    for k, (label, _) in GATES.items():
        if label == name:
            return k
    raise ValueError(f"unknown gate {name!r}; expected one of {list(GATES)}")


def gate_label(name: str) -> str:
    return GATES[gate_key(name)][0]


@dataclass
class RunSpec:
    """Everything a run needs, however it was collected."""

    datasets: List[str]
    stage: str = "screen"                    # screen | sweep | both
    max_trees: int = 0                       # 0 = every tree in each dataset
    workers: int = 4
    gate: str = "both"
    max_eta: float = MAX_ETA
    p_min: float = P_MIN
    p_points: int = P_POINTS
    reps: int = REPS
    num_gaps: int = 10
    min_split: int = 5
    # On by default: measured +21% on the L arm at m=1000 (4 p x 5 reps, 0.48 s -> 0.58 s),
    # and the diagnostics are O(m^2) per replicate against the eigensolve's O(m^3), so the
    # share falls as 1/m -- about 4% at m=6000. Cheaper than ever re-running for them.
    extra_metrics: bool = True
    prefix: str = ""
    display_mode: str = "progress"

    def __post_init__(self):
        self.gate = gate_key(self.gate)
        if self.stage not in ("screen", "sweep", "both"):
            raise ValueError(f"stage must be screen|sweep|both, got {self.stage!r}")

    @property
    def runs_screen(self) -> bool:
        return self.stage in ("screen", "both")

    @property
    def runs_sweep(self) -> bool:
        return self.stage in ("sweep", "both")

    def p_values(self) -> Optional[np.ndarray]:
        """The one grid this run sweeps. None for a screening-only run."""
        if not self.runs_sweep:
            return None
        return np.logspace(np.log10(self.p_min), 0, self.p_points)


def verdicts(dataset: str) -> Dict[str, dict]:
    """Screen rows for a dataset, keyed by tree id. Empty when it has not been screened."""
    path = screen_cache_path(dataset)
    if not path.exists():
        return {}
    return {r["tree"]: r for r in np.load(path, allow_pickle=True)["rows"]
            if "error" not in r}


def select_ids(ids: Sequence[str], rows: Dict[str, dict], gate: str,
               max_eta: float = 0.0) -> List[str]:
    """Trees passing the validity gate, and (when ``max_eta`` > 0) not too lopsided."""
    if not rows:
        return list(ids)
    keep = GATES[gate_key(gate)][1]
    out = [t for t in ids if t in rows and keep(rows[t])]
    if max_eta and max_eta > 0:
        out = [t for t in out
               if max(rows[t].get("eta_S", 0.0), rows[t].get("eta_B", 0.0)) <= max_eta]
    return out


def default_gate(dataset_names: Sequence[str],
                 verdicts_by: Dict[str, dict]) -> str:
    """Strictest gate that still selects a tree in EVERY chosen dataset."""
    for key, (_, fn) in GATES.items():
        if all(any(fn(r) for r in verdicts_by.get(name, {}).values())
               for name in dataset_names):
            return key
    return "all"


@dataclass
class PlanRow:
    """What one dataset contributes to a run, before any of it has run."""

    dataset: object                          # real_datasets.Dataset
    m: int
    ids: List[str]                           # after --limit, before any gate
    selected: List[str]                      # what the sweep would take today
    n_screened: int
    n_gate: int
    hours: float
    warning: str = ""

    @property
    def name(self) -> str:
        return self.dataset.name


def plan(spec: RunSpec, verdicts_by: Optional[Dict[str, dict]] = None) -> List[PlanRow]:
    """Per-dataset tree lists and cost, under one shared configuration."""
    if verdicts_by is None:
        verdicts_by = {n: verdicts(n) for n in spec.datasets}
    rows: List[PlanRow] = []
    for name in spec.datasets:
        ds = get_dataset(name)
        m, _ = ds.shape()
        all_ids = ds.ids()
        ids = all_ids[:spec.max_trees] if spec.max_trees else list(all_ids)
        v = verdicts_by.get(name, {})
        n_screened = sum(1 for t in all_ids if t in v)
        secs = 0.0
        if spec.runs_screen:
            secs += (len(ids) * (SCREEN_SECS_LARGE if m >= 6000 else SCREEN_SECS_SMALL)
                     / max(1, spec.workers))
        n_gate, selected = len(ids), list(ids)
        if spec.runs_sweep:
            n_gate = len(select_ids(ids, v, spec.gate))
            selected = select_ids(ids, v, spec.gate, spec.max_eta)
            secs += (len(selected) * spec.p_points * spec.reps
                     * SOLVE_SECS_AT_6000 * (m / 6000.0) ** 3)
        warn = ""
        if spec.runs_sweep and not spec.runs_screen and n_screened < len(all_ids):
            warn = (f"screening covers only {n_screened} of {len(all_ids)} trees, and the "
                    f"sweep can only use screened trees -- screen this dataset first")
        rows.append(PlanRow(dataset=ds, m=m, ids=ids, selected=selected,
                            n_screened=n_screened, n_gate=n_gate,
                            hours=secs / 3600.0, warning=warn))
    return rows


def _run_dataset(row: PlanRow, spec: RunSpec) -> List[str]:
    """Both stages on one dataset. Returns the ids this run covered."""
    ds = row.dataset
    if spec.runs_screen:
        run_real_eta_screen(row.ids, screen_cache_path(ds.name),
                            dataset=ds.name, workers=spec.workers)
        rows = verdicts(ds.name)
        if not rows:
            raise RuntimeError(
                f"{ds.name}: screening produced nothing usable -- see the errors above "
                f"and failures.csv in the run directory")
    if not spec.runs_sweep:
        return list(row.ids)

    # re-read the verdicts: on a "both" run the screen above is what decides the selection
    ids = select_ids(row.ids, verdicts(ds.name), spec.gate, spec.max_eta)
    if not ids:
        print(f"  nothing to sweep -- no tree passes [{gate_label(spec.gate)}]"
              f"{f' with eta <= {spec.max_eta:g}' if spec.max_eta else ''}. "
              f"Widen the gate (--dataset-rule any, or all).", flush=True)
        return []
    run_sweep(ids, sweep_cache_dir(ds.name), spec.p_values(), reps=spec.reps,
              num_gaps=spec.num_gaps, min_split=spec.min_split,
              dataset=ds.name, m=row.m, extra_metrics=spec.extra_metrics)
    return ids


def execute(spec: RunSpec, rows: Optional[List[PlanRow]] = None,
            on_dataset: Optional[Callable[[int, int, PlanRow], None]] = None) -> int:
    """Run every dataset in ``spec`` into ONE run directory. Returns an exit code.

    Interruption and per-dataset failure are both survivable: the per-tree caches hold
    whatever finished, the remaining datasets still run, and the export happens either
    way with the outcome recorded in ``summary.json``.
    """
    rows = rows if rows is not None else plan(spec)
    run_dir = new_run_dir(spec.prefix)
    write_config(run_dir, asdict(spec), spec.datasets, spec.stage)
    set_display_mode(spec.display_mode)
    setup_log_file(str(run_dir))
    print(f"run directory: {run_dir}", flush=True)

    tee = Tee(run_dir / "run.log")
    sys.stdout = tee
    covered: Dict[str, List[str]] = {}
    failures: List[str] = []
    status = "completed"
    try:
        for k, row in enumerate(rows, 1):
            if on_dataset is not None:
                on_dataset(k, len(rows), row)
            elif len(rows) > 1:
                print(f"\n=== {row.name} ({k}/{len(rows)})", flush=True)
            try:
                covered[row.name] = _run_dataset(row, spec)
            except KeyboardInterrupt:
                print(f"{row.name}: interrupted -- finished trees are kept", flush=True)
                covered[row.name] = list(row.selected)
                status = "interrupted"
                break
            except Exception:
                failures.append(row.name)
                print(f"ERROR: {row.name} failed:", flush=True)
                traceback.print_exc(file=sys.stdout)
                covered[row.name] = list(row.selected)
    finally:
        tee.close()
        close_log_file()

    if failures:
        status = "failed"
    selection = None
    if spec.runs_sweep:
        selection = {
            "rule": spec.gate,
            "rule_label": gate_label(spec.gate),
            "max_eta": spec.max_eta,
            "per_dataset": {r.name: {"screened": r.n_screened,
                                     "after_rule": r.n_gate,
                                     "after_max_eta": len(covered.get(r.name, []))}
                            for r in rows},
        }
    for f in export_run(run_dir, covered or {r.name: r.ids for r in rows},
                        status, ", ".join(failures), spec.p_values(), selection):
        print(f"  {f.name}")
    if failures:
        print(f"ERROR: {len(failures)} dataset(s) failed: {', '.join(failures)} "
              f"-- traceback in {run_dir / 'run.log'}")
    if status != "completed":
        print(f"run marked {status!r} in summary.json -- re-run the same command to "
              "continue from the cache")
    print(f"everything for this run is in {run_dir}")
    return 1 if failures else 0
