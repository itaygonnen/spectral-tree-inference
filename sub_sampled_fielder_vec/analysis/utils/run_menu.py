"""Interactive real-data menu for ``scripts/interactive_run.py``.

The launcher's other branch configures a *simulated* sweep (tree model, mutation rate,
sequence length). Real datasets have none of those knobs -- the alignment and the true
tree are given -- so this asks a different, shorter set of questions:

    which datasets  any <name>/{fasta,newick} under data/tree_sets/, several at a time
    which stage    screening (eta + validity per operator) or recovery sweep (vs p)
    parameters     asked ONCE and applied to every chosen dataset, so sizes stay
                   comparable -- the whole point of running m=1000 beside m=6000

This module is only the questions and the status panel. Selecting trees, planning,
running and exporting all live in :mod:`src.runners.experiment_run`, which
``scripts/run_sweep.py`` calls too -- so a menu session on a login node and a nohup'd
batch run produce the same run directory from the same cache. That runner knows nothing
about FASTA: this module's real job is turning the answers into a ``RunSpec`` and a list
of ``Source`` objects.

Lives beside the dataset helpers in ``analysis/utils/`` -- ``src/`` must never import
from ``analysis/``, so the dependency runs launcher -> analysis -> src, never back.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from dataclasses import replace

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, get_multi_choice, print_divider, print_error,
    print_header, print_option, print_success, print_warning)

from .real_datasets import describe_search, list_datasets, sweep_cache_dir  # noqa: E402
from src.runners.eta_buckets import (ETA_TARGETS, PER_BUCKET,       # noqa: E402
                                     panel_lines)
from src.runners.experiment_run import (DEFAULT_GATE, EARLY_STOP,   # noqa: E402
                                        GATES, MAX_ETA,
                                        P_MIN, P_POINTS, REPS, RunSpec, as_command,
                                        default_gate, execute, gate_label,
                                        gate_operators, plan, select_buckets, verdicts,
                                        already_swept)
from src.runners.operators import ALL_OPERATORS, ARM_OF, OPERATORS  # noqa: E402


def _cmd_width() -> int:
    """Width to wrap a printed command at, leaving room for the two-space indent."""
    import shutil
    return max(40, min(shutil.get_terminal_size(fallback=(80, 24)).columns, 100) - 8)


def _say(text: str, indent: str = "") -> None:
    """Print prose wrapped to the terminal.

    Long lines are not a cosmetic problem here: the line before a prompt used to be 160
    characters, so on an 80-column terminal it wrapped twice and the answer landed
    mid-sentence -- which reads as the program having crashed.
    """
    import shutil
    import textwrap

    width = max(40, min(shutil.get_terminal_size(fallback=(80, 24)).columns, 100))
    for line in textwrap.wrap(text, width=width, initial_indent=indent,
                              subsequent_indent=indent + "  ") or [indent]:
        print(line)


_SOURCES_BY_NAME: Dict[str, object] = {}     # set by _run so the panel can count
                                             # trees that are not screened yet


def _ask_config(names, is_screen: bool, verdicts_by: Dict[str, dict],
                pooled_as: Optional[tuple] = None) -> RunSpec:
    """One parameter set for every source: show the defaults, edit them only on request.

    ``pooled_as`` is ``(eta_targets, per_bucket)`` when the trees were drawn from the eta
    pool to those targets. The sweep then selects on the same terms they were built on,
    instead of overriding them with the paper's caps -- a run pooled at 10 per cell was
    otherwise cut to 10/10/3/3, silently discarding 14 of the 40 trees it had just
    simulated. Bucketing itself is the same code either way (``select_ids``); what the
    pool changes is only what the buckets are expected to hold.
    """
    spec = RunSpec(stage="screen" if is_screen else "sweep",
                   max_trees=0,                   # 0 = every tree in each dataset
                   workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
                   gate=default_gate(names, verdicts_by),
                   max_eta=MAX_ETA, p_min=P_MIN, p_points=P_POINTS, reps=REPS,
                   early_stop=EARLY_STOP,
                   eta_targets=(pooled_as[0] if pooled_as else ETA_TARGETS),
                   per_bucket=(pooled_as[1] if pooled_as else PER_BUCKET),
                   extra_metrics=True, prefix="", display_mode="progress")

    print_divider()
    print("configuration (every source):")
    if is_screen:
        print(f"  trees      all, {spec.workers} workers")
        print(f"  operators  {', '.join(spec.operators)}")
    else:
        # Only what varies and changes the run. The bucket settings used to be echoed
        # here as well, and they are the first line of the panel immediately below.
        _say(f"gate       [{gate_label(spec.gate)}], eta <= {spec.max_eta:g}",
             indent="  ")
        stop = (f"early stop {spec.early_stop}" if spec.early_stop else "no early stop")
        print(f"  p-grid     {spec.p_points} points, {spec.p_min:g} .. 1.0, "
              f"{spec.reps} reps, {stop}")
        print(f"  operators  {', '.join(spec.operators)}")

    print()
    if not confirm("Edit this configuration?", default=False):
        return spec

    cap = get_input("Max trees per source (blank = all)", default="")
    spec.max_trees = int(cap) if cap and cap.strip() else 0
    if is_screen:
        # screening has no free parameters, so it has one canonical output per dataset
        spec.workers = int(get_input("Workers", default=str(spec.workers)))
        return spec

    # Spell the counts out. A bare "2, 91" reads as a range, and even "2/66" hides that
    # the denominator a gate can act on is the SCREENED trees, not the dataset.
    # No counts on the options: the Screening status table above already gives coverage
    # and verdicts, and the plan printed after this shows what the choice actually selects.
    keys = list(GATES)
    labels = [gate_label(k) for k in keys]
    chosen = get_menu_choice("Reference partition must be a real tree edge under:",
                             labels, default_index=keys.index(spec.gate))
    spec.gate = keys[labels.index(chosen)]
    spec.max_eta = float(get_input(
        "Drop trees whose reference split is more lopsided than eta (0 = keep all)",
        default=str(spec.max_eta)))
    _ask_buckets(spec, names, verdicts_by)
    spec.p_min = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                 default=str(spec.p_min)))
    spec.p_points = int(get_input(
        f"p-grid points (log-spaced, {spec.p_min:g}..1)", default=str(spec.p_points)))
    spec.reps = int(get_input("Bootstrap reps per p", default=str(spec.reps)))
    spec.early_stop = int(get_input(
        "Settled points before an end of the grid is extrapolated (0 = measure every p)",
        default=str(spec.early_stop)))
    ops = get_input("Operators (comma-separated: " + ",".join(ALL_OPERATORS) + ")",
                    default=",".join(spec.operators))
    spec.operators = tuple(o.strip() for o in ops.split(",") if o.strip())
    spec.extra_metrics = confirm(
        "Also record sigma2, ||sub-full||_2 and the numerical ranks at every p?",
        default=True)
    name = get_input("Name for this run (blank = timestamp only)", default="")
    spec.prefix = (name or "").strip()
    spec.display_mode = get_menu_choice(
        "Output:", ["progress - bars, and every line in the run's experiment.log",
                    "debug - every line on the terminal too, no bars"],
        default_index=0).split(" ")[0]
    return spec


def _ask_buckets(spec: RunSpec, names, verdicts_by: Dict[str, dict]) -> None:
    """Show what each eta bucket holds, then take the caps. Mutates ``spec`` in place.

    This is where a run stops paying for trees it does not need, so it is a dialogue
    rather than a flag: on real data the paper's windows at +-1 keep a minority of the
    valid trees and some buckets come up empty, and the only way to find out is to look.
    The panel is redrawn after the tolerance so a choice can be made against the numbers
    it produces rather than remembered from before it.

    Sets ``spec.bucket``, ``spec.eta_tol``, ``spec.eta_targets`` and ``spec.per_bucket``.
    """
    spec.bucket = confirm(
        "Restrict the sweep to the paper's eta buckets?  (n sweeps every gated tree)",
        default=spec.bucket)
    if not spec.bucket:
        return
    raw = get_input("Eta buckets - comma-separated targets",
                    default=",".join(str(t) for t in spec.eta_targets))
    spec.eta_targets = tuple(int(x) for x in raw.split(",") if x.strip())
    if len(spec.per_bucket) != len(spec.eta_targets):
        # keep the parallel lists parallel: reuse the last cap for any new target
        last = spec.per_bucket[-1] if spec.per_bucket else 3
        spec.per_bucket = tuple(
            (spec.per_bucket[i] if i < len(spec.per_bucket) else last)
            for i in range(len(spec.eta_targets)))

    while True:
        tol = get_input("Bucket half-width (a tree joins the nearest target within it)",
                        default=f"{spec.eta_tol:g}")
        new_tol = float(tol)
        if new_tol == spec.eta_tol:
            break
        spec.eta_tol = new_tol          # redraw once at the new tolerance, then move on
        _bucket_panels(spec, names, verdicts_by, legend=False)
        if confirm("Keep this half-width?", default=True):
            break

    raw = get_input("Trees per bucket - one count per target, in order",
                    default=",".join(str(c) for c in spec.per_bucket))
    caps = [int(x) for x in raw.split(",") if x.strip()]
    if len(caps) == len(spec.eta_targets):
        spec.per_bucket = tuple(caps)
    else:
        print_warning(f"need {len(spec.eta_targets)} counts, got {len(caps)} -- keeping "
                      + ",".join(str(c) for c in spec.per_bucket))
    _bucket_panels(spec, names, verdicts_by, legend=False)


def _bucket_panels(spec: RunSpec, names, verdicts_by: Dict[str, dict],
                   legend: bool = True, is_screen: bool = False) -> None:
    """One occupancy panel per source, under the spec as it stands.

    Shown on every run. A tree's bucket is its MEASURED eta, so on a screening run the
    trees not yet screened have no bucket yet -- they are counted on their own line
    rather than guessed at, and the table then reads as progress towards the buckets
    instead of a selection.

    ``legend=False`` for a redraw: the dialogue shows each panel up to three times while
    the tolerance and the caps are settled, and three copies of the same footer pushes
    the numbers off the screen.
    """
    for i, name in enumerate(names):
        rows = verdicts_by.get(name) or {}
        n_all = len(_ids_of(name, names, verdicts_by))
        kept, buckets, _short = (select_buckets(sorted(rows), rows, spec)
                                 if rows else ([], {}, {}))
        # the screen row carries the taxon count, which the sweep-cache meta is keyed on
        m = next((r.get("m", 0) for r in rows.values() if r.get("m")), 0)
        done = already_swept(name, m, kept, spec) if kept else set()
        todo = len(kept) - len(done)
        head = (f"{name}  ->  {len(rows)} of {n_all} screened"
                if is_screen else
                f"{name}  ->  {todo} tree(s) to sweep"
                + (f"  ({len(done)} already done)" if done else ""))
        body = panel_lines(buckets, spec.eta_targets, spec.per_bucket, spec.eta_tol,
                           legend=(legend and i == 0), done=done)
        if is_screen and n_all > len(rows):
            # inserted after the data rows, before TOTAL/legend, so the column that
            # cannot be known yet is visible rather than silently missing
            at = next(j for j, ln in enumerate(body) if ln.startswith("TOTAL"))
            body.insert(at, f"{'unscreened':<11}{n_all - len(rows):>7}{'-':>7}{'-':>7}"
                            f"  bucket unknown until screened")
        _framed([head] + body, rule_after=0)


def _ids_of(name, names, verdicts_by) -> list:
    """Every tree id the source offers, screened or not. Falls back to the screened ones
    when the caller did not carry the source objects (the panel is also drawn from
    ``--dry-run``, which has them, and from the menu, which has the names)."""
    src = _SOURCES_BY_NAME.get(name)
    return list(src.ids) if src is not None else list(verdicts_by.get(name) or {})


def _framed(lines: List[str], rule_after: int = -1) -> None:
    """Draw ``lines`` inside a box, so a status panel is not mistaken for a menu.

    A box wider than the terminal is worse than no box: every line wraps and the borders
    land mid-sentence, which is what an 80-column terminal did to this panel once it
    carried a column pair per operator. When it will not fit, rule lines are used
    instead -- the panel still reads as a panel, and nothing is mangled.
    """
    import shutil

    width = max(len(x) for x in lines)
    cols = shutil.get_terminal_size(fallback=(80, 24)).columns
    if width + 4 > cols:
        rule = "─" * min(cols, width)
        print(rule)
        for i, line in enumerate(lines):
            print(line)
            if i == rule_after:
                print(rule)
        print(rule)
        print()
        return
    print("┌" + "─" * (width + 2) + "┐")
    for i, line in enumerate(lines):
        print(f"│ {line.ljust(width)} │")
        if i == rule_after:
            print("├" + "─" * (width + 2) + "┤")
    print("└" + "─" * (width + 2) + "┘")
    print()


def _print_status(sources, verdicts_by: Dict[str, dict]) -> None:
    """One line per source: how much is screened, how much is usable, how much is swept.

    This used to be a framed table of per-operator columns plus eight lines of legend --
    19 lines for two sources, most of it restating what the bucket panel below it says
    better. What survives is what the table alone could tell you: coverage, whether the
    trees are usable at all, and how lopsided they are.
    """
    print_header("Screening status")
    for c in sources:
        v, ids = verdicts_by[c.name], c.ids
        done = [t for t in ids if t in v]
        swept = (len(list(sweep_cache_dir(c.name).glob("*.npz")))
                 if sweep_cache_dir(c.name).is_dir() else 0)
        if not done:
            _say(f"{c.name}: {len(ids)} trees, none screened yet", indent="  ")
            continue
        key = f"eta_{gate_operators(DEFAULT_GATE)[0]}"
        etas = [v[t][key] for t in done if key in v[t]]
        n_ok = sum(1 for t in done if GATES[DEFAULT_GATE][1](v[t]))
        med = f", median eta {float(np.median(etas)):.1f}" if etas else ""
        _say(f"{c.name}: {len(done)}/{len(ids)} screened, {n_ok} usable"
             f"{med}, {swept} swept", indent="  ")
    print()


def _run(sources, *, batch_hint: str, pooled_as: Optional[tuple] = None,
         reconfigure=None) -> None:
    """Status, stage, one configuration, the plan, then the run. Source-agnostic.

    ``reconfigure`` is offered as a third choice on the stage menu and returns a fresh
    ``(sources, pooled_as)`` or None to quit. It is a menu entry rather than a prompt so
    that NOTHING is asked before the status is on screen: changing the settings is a
    thing you do after seeing that the defaults do not suit, not before.
    """
    while True:
        verdicts_by = {c.name: verdicts(c.name) for c in sources}
        _SOURCES_BY_NAME.clear()
        _SOURCES_BY_NAME.update({c.name: c for c in sources})
        _print_status(sources, verdicts_by)

        # The bucket table FIRST -- straight after the source is picked, before any
        # question at all. It is the status this launcher exists to show, and burying it
        # behind prompts meant answering several to find out whether there was anything
        # worth running. Occupancy does not depend on the stage, so there is nothing to
        # wait for.
        _preview = RunSpec(stage="sweep", gate=default_gate(
            [c.name for c in sources], verdicts_by),
            **(dict(eta_targets=pooled_as[0], per_bucket=pooled_as[1])
               if pooled_as else {}))
        _bucket_panels(_preview, [c.name for c in sources], verdicts_by)

        options = ["screening - split the full matrix, check eta and validity",
                   "recovery sweep - sweep trees over p for NMI curve "
                   "(needs screening)"]
        if reconfigure is not None:
            options.append("change the data settings - model, size, pooling")
        stage = get_menu_choice("Run:", options, default_index=0)
        if reconfigure is not None and stage.startswith("change"):
            got = reconfigure()
            if not got or not got[0]:
                return
            sources, pooled_as = got
            continue
        break
    is_screen = stage.startswith("screening")

    spec = _ask_config([c.name for c in sources], is_screen, verdicts_by, pooled_as)
    if spec.max_trees:
        sources = [replace(c, ids=c.ids[:spec.max_trees]) for c in sources]
    rows = plan(spec, sources, verdicts_by)

    print_divider()
    if not is_screen:
        pv = spec.p_values()
        print(f"  p-grid     {len(pv)} log-spaced points, {pv[0]:.4g} .. {pv[-1]:.4g}")
    width = max([12] + [len(r.name) + 1 for r in rows])
    total = sum(r.hours for r in rows)
    for r in rows:
        print(f"  {r.name:<{width}} m={r.m:<5} {len(r.selected):>4} trees"
              f"  ~{r.hours:.1f} h")
        if not is_screen:
            # where the trees went: the two filters, in the order they are applied
            _say(f"{r.n_screened} screened -> {r.n_gate} "
                 f"[{gate_label(spec.gate)}] -> {len(r.selected)} with eta <= "
                 f"{spec.max_eta:g}", indent="      ")
        if r.warning:
            print_warning(f"{r.name}: {r.warning}")
    print()
    _say(f"total ~{total:.1f} h"
         + ("" if is_screen else " (the sweep runs one tree at a time; workers "
                                "parallelise the screen only)")
         + ". Every tree is cached on its own, so this is safe to interrupt and "
           "resume.")
    if not confirm("Run it here now?  (n prints the command instead)",
                   default=total <= 1.0):
        # "n" is not a cancellation, it is a request: the command carrying the answers
        # just given is the whole point of having answered them. Printing a warning here
        # made choosing "run it elsewhere" look like a failure.
        print()
        print_success("Not running here. This is the run you configured:")
        print()
        print("  " + as_command(spec, sources, width=_cmd_width(), indent="      "))
        print()
        print("  # detached, with a log:")
        print("  nohup " + as_command(spec, sources, width=_cmd_width(),
                                      indent="      ")
              + f" \\\n      > logs/{spec.stage}.log 2>&1 &")
        print()
        _say("Add --dry-run to either one to see the selection without running it. "
             "Nothing has been written.")
        return

    def _announce(k, n, row):
        print_divider()
        print_header(f"[{k}/{n}] {row.name}")

    execute(spec, sources, rows, on_source=_announce)
    print_success("done")


def run_real_data_menu() -> None:
    """Real datasets: pick the FASTA/newick directories, then run."""
    datasets = list_datasets()
    if not datasets:
        print_error("No real datasets found.")
        print(describe_search())
        return

    print_header("Real data")
    _say("Tip: select several with commas (e.g. '1,2') to run every size in turn",
         indent="  ")
    print()
    for i, c in enumerate(datasets, 1):
        m, seq_len = c.shape()
        print_option(str(i), f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})",
                     highlight=(i == len(datasets)))
    print()
    picks = get_multi_choice("Dataset(s)", [str(i) for i in range(1, len(datasets) + 1)])
    chosen = [datasets[int(i) - 1] for i in picks]
    _run([c.source() for c in chosen], batch_hint="python scripts/run_sweep.py")


def run_generated_menu() -> None:
    """Simulated trees: ask the model, sizes and sequence length, then run.

    The same experiment as the real branch -- a population of trees, screened, gated and
    swept, with the curve a median over trees. The simulated branch used to run one tree
    per size and take its spread from bootstrap replicates, which is a different
    quantity and was not comparable with a real-data curve.
    """
    from src.runners.generated_source import sources_from_plan
    from src.utils.generated_prompts import (baseline_plan, describe_plan,
                                             prompt_generated_plan, resolve_pool_ids)

    def _build(plan_g):
        """A plan -> (sources, pooled_as), or (None, None) when it yields no trees."""
        _say(describe_plan(plan_g), indent="  ")
        if plan_g.pooled:
            # the eta pool is the only place a rejection-sampled tree exists; top it
            # up before anything asks for ids, because what comes back is what is on
            # disk
            plan_g.tree_ids = resolve_pool_ids(plan_g)
            if not plan_g.tree_ids:
                print_error("the eta pool is empty for these parameters and nothing "
                            "could be built -- nothing to run")
                return None, None
        srcs = sources_from_plan(plan_g)
        if not srcs:
            print_error("no trees to run")
            return None, None
        # selected on the terms they were pooled on, not the paper's caps
        pooled = ((tuple(plan_g.etas), (plan_g.per_cell,) * len(plan_g.etas))
                  if plan_g.pooled and plan_g.etas else None)
        return srcs, pooled

    print_header("Generated data")
    # Start from the model's baseline -- no questions -- so the status is the first
    # thing on screen. The six-step dialogue is the "change the data settings" entry on
    # the menu below.
    sources, pooled_as = _build(baseline_plan())
    if not sources:
        return
    _run(sources, batch_hint="", pooled_as=pooled_as,
         reconfigure=lambda: _build(prompt_generated_plan()))
