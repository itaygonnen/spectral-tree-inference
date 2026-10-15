"""One run: several sources, two stages, one result directory.

A run is the same object whether its trees come from FASTA files or from a simulator.
It holds the knobs (:class:`RunSpec`), a list of :class:`Source` objects that say how to
enumerate and load trees, and nothing else -- no ``Dataset``, no ``.fasta``, no
simulator. Front-ends differ only in how they fill those in:

    scripts/run_sweep.py         argparse -> RunSpec + real sources
    analysis/utils/run_menu.py   a menu -> the same, for real OR simulated sources
    scripts/run_benchmark.py     a menu -> RunSpec + generated sources

Before this existed the menu and the command line carried a selection implementation
each and wrote different vocabularies into the same ``summary.json`` field, and the
real and generated pipelines were two separate stacks doing the same thing.
"""
from __future__ import annotations

import sys
import traceback
from dataclasses import asdict, dataclass, field
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from ..utils.logging import close_log_file, set_display_mode, setup_log_file
from ..utils.run_paths import new_run_dir, screen_cache_path, sweep_cache_dir
from .operator_screen import run_screen
from .operator_sweep import run_sweep
from .eta_buckets import ETA_TARGETS, ETA_TOL, PER_BUCKET, label, panel_lines
from .eta_buckets import select as bucket_select
from .operators import ALL_OPERATORS, ARM_OF, OPERATORS, resolve

# Which screened trees a sweep may use. The gate is whether that operator's partition of
# the FULL matrix is a real single-edge split of the true tree -- recovery towards a
# reference that is not a tree edge measures stability, not correctness.
#
# The keys are the machine names (they go into summary.json and onto the command line);
# the labels are what a menu shows. "both" means L(S) and B -- the pair the figures
# compare -- and keeps that meaning now that L_sym is screened too, so a selection
# recorded before L_sym existed still means what it said.
# Third element: the operators the imbalance cap weighs. It is the gate's operators, not
# the run's. A run sweeps L_sym as a third curve; capping on it too would drop a tree whose
# L and B splits are perfectly even, which is a loss nobody asked for and nobody could
# explain from the output. "any"/"all" name no operator, so they fall back to the pair the
# recovery figure compares -- which is also what the cap meant before L_sym was screened.
GATES: Dict[str, tuple] = {
    "both":       ("L(S) and B both cut a real tree edge",
                   lambda r: bool(r.get("valid_S")) and bool(r.get("valid_B")),
                   ("S", "B")),
    "valid_L":    ("L(S) cuts a real tree edge", lambda r: bool(r.get("valid_S")),
                   ("S",)),
    "valid_Lsym": ("L_sym cuts a real tree edge", lambda r: bool(r.get("valid_Lsym")),
                   ("Lsym",)),
    "valid_B":    ("B cuts a real tree edge", lambda r: bool(r.get("valid_B")),
                   ("B",)),
    "any":        ("any operator cuts a real tree edge",
                   lambda r: any(r.get(f"valid_{k}") for k in ALL_OPERATORS),
                   ("S", "B")),
    "all":        ("every tree, valid or not", lambda r: True, ("S", "B")),
}
# what the CLI called the L gate before the operators were named by arm
GATE_ALIASES = {"valid_S": "valid_L"}

# Reference splits this lopsided are not a recovery question. At m=1000 the L(S) k-means
# cut routinely isolates ONE taxon (1/999, eta=999): a real pendant edge, so the validity
# gate passes it, but every sub-sample then picks a different singleton and NMI sits at 0
# for every p. 20 is well above the recoverability scale (m/log m)^(1/3) ~ 5.3 at m=1000.
MAX_ETA = 20.0

# What a sweep reads by default: the similarity arm alone (2026-10-14 -- the distance
# route is no longer run). NOT ``ALL_OPERATORS``: L_sym is a full Fiedler arm costing
# about as much as L, and B is ``GRIFFING_COST_SHARE`` of one. Both stay in
# ``ALL_OPERATORS`` because the screen is cheap and still records their validity and eta,
# so the ``valid_Lsym`` / ``valid_B`` / ``both`` gates keep working, every screen cache
# keeps its meaning, and ``--operators S,Lsym,B`` still sweeps all three on request.
# ``DEFAULT_GATE`` moves with this: gating on B's validity while never sweeping B throws
# away trees that are perfectly good for L.
DEFAULT_OPERATORS: tuple = ("S",)

# Consecutive settled p's required before either end of the grid is extrapolated. The
# sweep visits the grid middle-out and stops each side on its own rule -- NMI == 1.0
# above, NMI <= LOW_FLOOR below -- so a curve that has said what it is going to say is
# not measured to both horizons. 0 measures every p. See
# ``src.runners.p_sweep_inner.MiddleOut`` for the traversal and the fill contract.
EARLY_STOP = 2

# The gate a run uses unless told otherwise. ``valid_L`` rather than ``both`` because the
# sweep default is S alone; ``gate_operators("valid_L")`` is ``("S",)``, so the imbalance
# cap and the eta buckets both weigh the similarity arm without anything else being said.
DEFAULT_GATE = "valid_L"

# p starts at 1e-4: the transition sits near log n / n, which is 1.4e-3 at m=6000, so the
# grid has to reach below that for the curve to show a floor rather than start on the ramp.
P_MIN, P_POINTS, REPS = 1e-4, 20, 10

# Per-tree cost, measured at m=6000: a screen is one load plus an eigensolve per operator
# (~75 s), a sweep is p x reps sub-sampled solves at ~9 s each. Both are dominated by the
# O(m^3) eigensolve, so both scale by (m/6000)^3 -- a flat estimate under m=6000 was ~100x
# over at n=128 and under at n=4096, and the "~X h" line is what people plan around.
# B is ~30x cheaper than a Fiedler arm (ARPACK for one eigenpair, no Laplacian).
SCREEN_SECS_AT_6000, SOLVE_SECS_AT_6000 = 75.0, 9.0
# ...and the screen splits into building the matrices, O(m^2 L), and the eigensolves,
# O(m^3): at m=6000/L=5000 that measured 60 s and 15 s. Keeping the two terms apart is
# what stops the estimate being ~30x over at n=128, where the old flat 10 s/tree sat.
SCREEN_BUILD_AT_6000, SCREEN_SOLVE_AT_6000 = 60.0, 15.0
REF_M, REF_L = 6000.0, 5000.0
GRIFFING_COST_SHARE = 0.03


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


def as_command(spec: "RunSpec", sources: Sequence["Source"],
               script: str = "python scripts/run_sweep.py",
               width: int = 0, indent: str = "    ") -> str:
    """The command line that reproduces this spec exactly.

    The menu used to suggest a bare ``--dataset X --stage sweep`` for long runs, which
    silently dropped every answer the user had just given -- at m=6000 that turned a
    38-tree ``valid_L`` sweep into a ``both`` sweep of nothing. Every field that differs
    from the CLI default is emitted, so the printed command IS the configured run.
    """
    parts = [script, f'--dataset "{",".join(s.name for s in sources)}"',
             f"--stage {spec.stage}"]
    if spec.max_trees:
        parts.append(f"--limit {spec.max_trees}")
    if spec.workers != 4:
        parts.append(f"--workers {spec.workers}")
    if spec.runs_sweep:
        parts.append(f"--dataset-rule {spec.gate}")
        parts.append(f"--max-eta {spec.max_eta:g}")
        # the bucket answers are the ones that decide how many trees run at all, so a
        # printed command that dropped them would reproduce a much longer run
        if not spec.bucket:
            parts.append("--no-buckets")
        else:
            if tuple(spec.eta_targets) != tuple(ETA_TARGETS):
                parts.append("--eta-targets " + ",".join(str(t) for t in spec.eta_targets))
            if tuple(spec.per_bucket) != tuple(PER_BUCKET):
                parts.append("--per-bucket " + ",".join(str(c) for c in spec.per_bucket))
            if float(spec.eta_tol) != float(ETA_TOL):
                parts.append(f"--eta-tol {spec.eta_tol:g}")
        parts.append(f"--p-min {spec.p_min:g}")
        parts.append(f"--p-points {spec.p_points}")
        parts.append(f"--reps {spec.reps}")
        if tuple(spec.operators) != tuple(DEFAULT_OPERATORS):
            parts.append(f"--operators {','.join(spec.operators)}")
        if spec.early_stop != EARLY_STOP:
            parts.append(f"--early-stop {spec.early_stop}")
        if not spec.extra_metrics:
            parts.append("--no-extra-metrics")
    if spec.prefix:
        parts.append(f'--prefix "{spec.prefix}"')
    if spec.display_mode != "progress":
        parts.append(f"--display-mode {spec.display_mode}")
    if width <= 0:
        return " ".join(parts)
    # backslash continuations at flag boundaries: a command that wraps mid-token in the
    # terminal is still copyable, but it is not readable, and this one is printed to be
    # read as much as pasted
    lines, cur = [], parts[0]
    for part in parts[1:]:
        if len(cur) + 1 + len(part) + 2 > width:
            lines.append(cur + " \\")
            cur = indent + part
        else:
            cur += " " + part
    lines.append(cur)
    return "\n".join(lines)


def screen_secs(m: int, seq_len: int = 0) -> float:
    """Rough seconds to screen one tree: matrix build O(m^2 L) plus eigensolves O(m^3).

    Calibrated on the m=6000 / L=5000 dataset. It is an estimate for the "~X h" line, not
    a measurement -- a real dataset's FASTA parse adds a term this does not model, so it
    reads low at small m.
    """
    scale_l = (seq_len / REF_L) if seq_len else 1.0
    return max(0.05,
               SCREEN_BUILD_AT_6000 * (m / REF_M) ** 2 * scale_l
               + SCREEN_SOLVE_AT_6000 * (m / REF_M) ** 3)


def gate_operators(name: str) -> tuple:
    """The operators a gate's imbalance cap weighs. See the note on ``GATES``."""
    return GATES[gate_key(name)][2]


@dataclass
class Source:
    """Where one group of trees comes from, and what to call it.

    ``name`` is the cache slug and the CSV column, ``loader`` turns a tree id into
    ``(S, labels, tree, D)``, ``ids`` is every tree available, ``m`` its taxon count and
    ``seq_len`` its alignment length (only the cost estimate uses that).
    Anything that satisfies those four is a source: a FASTA directory, a simulator, a
    future format nobody has written yet.
    """

    name: str
    loader: Callable
    ids: List[str] = field(default_factory=list)
    m: int = 0
    seq_len: int = 0


@dataclass
class RunSpec:
    """Everything a run needs, however it was collected. No data source in sight."""

    stage: str = "screen"                    # screen | sweep | both
    max_trees: int = 0                       # 0 = every tree in each source
    workers: int = 4
    gate: str = DEFAULT_GATE
    max_eta: float = MAX_ETA
    # Which trees are worth sweeping at all. Screened trees are bucketed by the eta the
    # gate weighs, against the paper's targets, and only ``per_bucket`` of each are swept
    # -- everything landing in no bucket is dropped, which is where the time goes. See
    # ``src.runners.eta_buckets``; parallel tuples rather than a dict so ``asdict`` into
    # config.json round-trips (JSON would key a dict by str).
    bucket: bool = True
    eta_targets: Sequence[int] = ETA_TARGETS
    per_bucket: Sequence[int] = PER_BUCKET
    eta_tol: float = ETA_TOL
    p_min: float = P_MIN
    p_points: int = P_POINTS
    reps: int = REPS
    num_gaps: int = 10
    min_split: int = 5
    operators: Sequence[str] = DEFAULT_OPERATORS
    early_stop: int = EARLY_STOP
    # On by default: measured +21% on the L arm at m=1000 (4 p x 5 reps, 0.48 s ->
    # 0.58 s), and the diagnostics are O(m^2) per replicate against the eigensolve's
    # O(m^3), so the share falls as 1/m -- about 4% at m=6000.
    extra_metrics: bool = True
    prefix: str = ""
    display_mode: str = "progress"

    def __post_init__(self):
        self.gate = gate_key(self.gate)
        self.operators = resolve(self.operators)
        if self.stage not in ("screen", "sweep", "both"):
            raise ValueError(f"stage must be screen|sweep|both, got {self.stage!r}")
        self.eta_targets = tuple(int(t) for t in self.eta_targets)
        self.per_bucket = tuple(int(c) for c in self.per_bucket)
        if len(self.per_bucket) != len(self.eta_targets):
            raise ValueError(
                f"per_bucket has {len(self.per_bucket)} entries for "
                f"{len(self.eta_targets)} eta targets; they are parallel lists")
        if self.eta_tol <= 0:
            raise ValueError(f"eta_tol must be > 0, got {self.eta_tol}")

    @property
    def runs_screen(self) -> bool:
        return self.stage in ("screen", "both")

    @property
    def runs_sweep(self) -> bool:
        return self.stage in ("sweep", "both")

    @property
    def arms(self) -> List[str]:
        return [ARM_OF[k] for k in self.operators]

    def p_values(self) -> Optional[np.ndarray]:
        """The one grid this run sweeps. None for a screening-only run."""
        if not self.runs_sweep:
            return None
        return np.logspace(np.log10(self.p_min), 0, self.p_points)


def verdicts(source: str) -> Dict[str, dict]:
    """Screen rows for a source, keyed by tree id. Empty when it has not been screened."""
    path = screen_cache_path(source)
    if not path.exists():
        return {}
    return {r["tree"]: r for r in np.load(path, allow_pickle=True)["rows"]
            if "error" not in r}


def select_ids(ids: Sequence[str], rows: Dict[str, dict], gate: str,
               max_eta: float = 0.0,
               operators: Sequence[str] = (),
               bucket: bool = False,
               eta_targets: Sequence[int] = ETA_TARGETS,
               per_bucket: Sequence[int] = PER_BUCKET,
               eta_tol: float = ETA_TOL) -> List[str]:
    """What a sweep would take: three filters, in order.

        screened -> [validity gate] -> [eta <= max_eta] -> [per-bucket cap] -> swept

    The imbalance cap weighs the operators the GATE names, not every operator the run
    happens to sweep -- see the note on :data:`GATES`. ``operators`` is accepted and
    ignored; it is kept so the older call signature does not break.

    ``bucket`` adds the third filter (:mod:`src.runners.eta_buckets`): keep only
    ``per_bucket`` trees per eta target, and drop everything that matches no target. It
    defaults OFF here and ON in :class:`RunSpec`, because this function is also called
    bare -- by the selftest and by anything asking "what does this gate select" -- and
    those callers mean the gate, not a run's budget. Note the bucket windows
    (1/5/10/15 +-1) are far tighter than ``MAX_ETA=20``, so when it is on the eta cap is
    subsumed and left in place only so turning buckets off restores the old behaviour.

    This is the one function ``plan()``, ``_run_source()`` and ``summary.json`` all agree
    on, which is why the bucket stage lives here and not in a front-end: real and
    generated sources are fixed by the same edit, and ``PlanRow.hours`` -- linear in the
    selected count -- corrects itself.
    """
    if not rows:
        return list(ids)
    _, keep, cap_ops = GATES[gate_key(gate)]
    out = [t for t in ids if t in rows and keep(rows[t])]
    if max_eta and max_eta > 0:
        keys = [f"eta_{k}" for k in cap_ops]
        out = [t for t in out
               if max((rows[t].get(k, 0.0) for k in keys), default=0.0) <= max_eta]
    if bucket:
        out, _buckets, _short = bucket_select(out, rows, cap_ops, eta_targets,
                                              per_bucket, eta_tol)
    return out


def select_buckets(ids: Sequence[str], rows: Dict[str, dict], spec: "RunSpec"):
    """``(kept, buckets, shortfall)`` for one source under ``spec``, gate and cap applied.

    What :func:`select_ids` does, with the occupancy it selected from handed back so a
    panel can be drawn from the same numbers the plan used rather than a second
    derivation that could disagree.
    """
    gated = select_ids(ids, rows, spec.gate, spec.max_eta, spec.operators)
    if not rows:
        return list(gated), {}, {}
    cap_ops = gate_operators(spec.gate)
    return bucket_select(gated, rows, cap_ops, spec.eta_targets, spec.per_bucket,
                         spec.eta_tol)


def default_gate(source_names: Sequence[str], verdicts_by: Dict[str, dict]) -> str:
    """Strictest gate that still selects a tree in EVERY chosen source."""
    for key, (_, fn, _ops) in GATES.items():
        if all(any(fn(r) for r in verdicts_by.get(name, {}).values())
               for name in source_names):
            return key
    return "all"


@dataclass
class PlanRow:
    """What one source contributes to a run, before any of it has run."""

    source: Source
    ids: List[str]                           # after --limit, before any gate
    selected: List[str]                      # what the sweep would take today
    n_screened: int
    n_gate: int
    hours: float
    warning: str = ""
    n_max_eta: int = 0                       # survivors of the eta cap, before buckets
    # ``{target: ids}`` plus a ``None`` key, the occupancy ``selected`` was drawn from.
    # Carried on the row so the panel, the plan line and summary.json are one derivation
    # -- three independent ones would eventually disagree.
    buckets: Dict = field(default_factory=dict)
    shortfall: Dict = field(default_factory=dict)
    n_done: int = 0                          # selected trees already swept on this grid

    @property
    def name(self) -> str:
        return self.source.name

    @property
    def m(self) -> int:
        return self.source.m


def plan(spec: RunSpec, sources: Sequence[Source],
         verdicts_by: Optional[Dict[str, dict]] = None) -> List[PlanRow]:
    """Per-source tree lists and cost, under one shared configuration."""
    if verdicts_by is None:
        verdicts_by = {s.name: verdicts(s.name) for s in sources}
    n_fiedler = sum(1 for k in spec.operators if OPERATORS[k].kind == "fiedler")
    n_griffing = len(spec.operators) - n_fiedler
    rows: List[PlanRow] = []
    for src in sources:
        ids = src.ids[:spec.max_trees] if spec.max_trees else list(src.ids)
        v = verdicts_by.get(src.name, {})
        n_screened = sum(1 for t in src.ids if t in v)
        secs = 0.0
        if spec.runs_screen:
            per_tree = screen_secs(src.m, src.seq_len)
            secs += len(ids) * per_tree / max(1, spec.workers)
        n_gate, selected = len(ids), list(ids)
        n_max_eta, buckets, shortfall, n_done = len(ids), {}, {}, 0
        if spec.runs_sweep:
            n_gate = len(select_ids(ids, v, spec.gate))
            capped = select_ids(ids, v, spec.gate, spec.max_eta, spec.operators)
            n_max_eta = len(capped)
            if spec.bucket and v:
                selected, buckets, shortfall = select_buckets(ids, v, spec)
            else:
                selected = capped
            # cheap: a header read per selected tree, and the selection is capped
            n_done = len(already_swept(src.name, src.m, selected, spec))
            # serial: run_sweep walks the trees one at a time, so this is wall time.
            # ``workers`` parallelises the screen only. An UPPER bound when
            # ``early_stop`` is on: how many p's each end skips depends on where that
            # tree's curve turns, which is not knowable before sweeping it.
            # Trees already swept on this grid are reused, so they cost nothing -- and
            # the bucket panel prints that count right beside this estimate, so an
            # estimate that ignored it would be visibly contradicted on screen.
            secs += (max(0, len(selected) - n_done) * spec.p_points * spec.reps
                     * SOLVE_SECS_AT_6000 * (src.m / 6000.0) ** 3
                     * (n_fiedler + n_griffing * GRIFFING_COST_SHARE))
        warn = ""
        if spec.runs_sweep and not spec.runs_screen and n_screened < len(src.ids):
            warn = (f"screening covers only {n_screened} of {len(src.ids)} trees, and "
                    f"the sweep can only use screened trees -- screen this source first")
        rows.append(PlanRow(source=src, ids=ids, selected=selected,
                            n_screened=n_screened, n_gate=n_gate,
                            hours=secs / 3600.0, warning=warn,
                            n_max_eta=n_max_eta, buckets=buckets,
                            shortfall=shortfall, n_done=n_done))
    return rows


def _run_source(row: PlanRow, spec: RunSpec) -> List[str]:
    """Both stages on one source. Returns the ids this run covered."""
    src = row.source
    if spec.runs_screen:
        run_screen(row.ids, screen_cache_path(src.name), src.loader,
                   operators=spec.operators, source=src.name, workers=spec.workers,
                   min_split=spec.min_split)
        if not verdicts(src.name):
            raise RuntimeError(
                f"{src.name}: screening produced nothing usable -- see the errors above "
                f"and failures.csv in the run directory")
    if not spec.runs_sweep:
        return list(row.ids)

    # re-read the verdicts: on a "both" run the screen above is what decides this
    ids = select_ids(row.ids, verdicts(src.name), spec.gate, spec.max_eta,
                     spec.operators, spec.bucket, spec.eta_targets, spec.per_bucket,
                     spec.eta_tol)
    if not ids:
        print(f"  nothing to sweep -- no tree passes [{gate_label(spec.gate)}]"
              f"{f' with eta <= {spec.max_eta:g}' if spec.max_eta else ''}. "
              f"Widen the gate (--dataset-rule any, or all).", flush=True)
        return []
    run_sweep(ids, sweep_cache_dir(src.name), src.loader, spec.p_values(),
              reps=spec.reps, num_gaps=spec.num_gaps, min_split=spec.min_split,
              operators=spec.operators, source=src.name, m=row.m,
              extra_metrics=spec.extra_metrics, early_stop=spec.early_stop)
    return ids


def already_swept(source: str, m: int, selected: Sequence[str],
                  spec: RunSpec) -> set:
    """Of this source's selected trees, those a sweep would REUSE rather than recompute.

    Not a glob of ``*.npz``: ``load_tree_sweep`` only reuses a cached tree whose p-grid
    and meta match this run's, so a tree swept on another grid is on disk and still has
    to be recomputed. Measured on 1000 taxa, two eta15 trees had a cache file and one of
    them was reusable -- counting the file would have promised half the work that
    remained.
    """
    from ..utils.run_paths import sweep_cache_dir
    from .operator_sweep import load_tree_sweep, sweep_meta

    pv = spec.p_values()
    if pv is None or not selected:
        return set()
    cache = sweep_cache_dir(source)
    if not cache.is_dir():
        return set()
    grid = np.round(np.asarray(pv, float), 6)
    meta = dict(sweep_meta(spec.reps, spec.num_gaps, spec.min_split, int(m),
                           spec.early_stop))
    return {t for t in selected
            if load_tree_sweep(cache, t, grid, meta) is not None}


def bucket_panel(row: "PlanRow", spec: RunSpec, legend: bool = True) -> List[str]:
    """The bucket occupancy for one planned source, or no lines when it does not apply.

    Lives here, beside the selection it describes, so the menu prompt and the
    non-interactive ``--dry-run`` print the same table from the same numbers.
    """
    if not (spec.bucket and spec.runs_sweep and row.buckets):
        return []
    return panel_lines(row.buckets, spec.eta_targets, spec.per_bucket, spec.eta_tol,
                       legend=legend,
                       done=already_swept(row.source.name, row.m, row.selected, spec))


def _selection_row(row: "PlanRow", spec: RunSpec, covered: Sequence[str]) -> dict:
    """One dataset's funnel for ``summary.json``: how many trees survived each filter,
    and -- when bucketing is on -- what each bucket held against what was asked for.

    Without the per-bucket got/want a summary can say "18 swept" while leaving a reader
    to assume the four buckets were even, which on real data they are not.
    """
    out = {"screened": row.n_screened, "after_rule": row.n_gate,
           "after_max_eta": row.n_max_eta, "swept": len(covered)}
    if spec.bucket:
        out["buckets"] = {label(t): {"available": len(ids),
                                     "want": int(c),
                                     "take": min(len(ids), int(c))}
                          for t, c, ids in zip(spec.eta_targets, spec.per_bucket,
                                               (row.buckets.get(int(t), [])
                                                for t in spec.eta_targets))}
        out["unbucketed"] = len(row.buckets.get(None, []))
    return out


def execute(spec: RunSpec, sources: Sequence[Source],
            rows: Optional[List[PlanRow]] = None,
            on_source: Optional[Callable[[int, int, PlanRow], None]] = None) -> int:
    """Run every source in ``sources`` into ONE run directory. Returns an exit code.

    Interruption and per-source failure are both survivable: the per-tree caches hold
    whatever finished, the remaining sources still run, and the export happens either
    way with the outcome recorded in ``summary.json``.
    """
    from ..utils.run_export import Tee, export_run, write_config

    rows = rows if rows is not None else plan(spec, sources)
    names = [s.name for s in sources]
    run_dir = new_run_dir(spec.prefix)
    write_config(run_dir, asdict(spec), names, spec.stage)
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
            if on_source is not None:
                on_source(k, len(rows), row)
            elif len(rows) > 1:
                print(f"\n=== {row.name} ({k}/{len(rows)})", flush=True)
            try:
                covered[row.name] = _run_source(row, spec)
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
            "max_eta_operators": list(gate_operators(spec.gate)),
            "operators": list(spec.operators),
            "buckets": {"enabled": bool(spec.bucket),
                        "eta_targets": [int(t) for t in spec.eta_targets],
                        "per_bucket": [int(c) for c in spec.per_bucket],
                        "eta_tol": float(spec.eta_tol)},
            "per_dataset": {r.name: _selection_row(r, spec, covered.get(r.name, []))
                            for r in rows},
        }
    for f in export_run(run_dir, covered or {r.name: r.ids for r in rows},
                        status, ", ".join(failures), spec.p_values(), selection,
                        spec.operators):
        print(f"  {f.name}")
    if failures:
        print(f"ERROR: {len(failures)} dataset(s) failed: {', '.join(failures)} "
              f"-- traceback in {run_dir / 'run.log'}")
    if status != "completed":
        print(f"run marked {status!r} in summary.json -- re-run the same command to "
              "continue from the cache")
    print(f"everything for this run is in {run_dir}")
    return 1 if failures else 0
