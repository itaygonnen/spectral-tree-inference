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
from typing import Dict, List

import numpy as np

from dataclasses import replace

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, get_multi_choice, print_divider, print_error,
    print_header, print_option, print_success, print_warning)

from .real_datasets import describe_search, list_datasets, sweep_cache_dir  # noqa: E402
from src.runners.experiment_run import (GATES, MAX_ETA, P_MIN,  # noqa: E402
                                        P_POINTS, REPS, RunSpec, as_command,
                                        default_gate, execute, gate_label,
                                        gate_operators, plan, verdicts)
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


def _ask_config(names, is_screen: bool, verdicts_by: Dict[str, dict]) -> RunSpec:
    """One parameter set for every source: show the defaults, edit them only on request."""
    spec = RunSpec(stage="screen" if is_screen else "sweep",
                   max_trees=0,                   # 0 = every tree in each dataset
                   workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
                   gate=default_gate(names, verdicts_by),
                   max_eta=MAX_ETA, p_min=P_MIN, p_points=P_POINTS, reps=REPS,
                   extra_metrics=True, prefix="", display_mode="progress")

    print_divider()
    print("configuration (applies to every source chosen):")
    if is_screen:
        print("  trees      all")
        print(f"  workers    {spec.workers}")
        print(f"  operators  {', '.join(OPERATORS[k].key for k in spec.operators)}"
              "  (each Fiedler arm cut by k-means or sign, whichever is more even)")
    else:
        _say(f"trees      all, gated on [{gate_label(spec.gate)}], "
             f"eta <= {spec.max_eta:g}", indent="  ")
        print(f"  p-grid     {spec.p_points} log-spaced points, {spec.p_min:g} .. 1.0")
        print(f"  reps       {spec.reps} bootstrap replicates per p")
        print(f"  operators  {', '.join(spec.operators)}"
              "   (NMI is what the figure plots)")
        print("  metrics    NMI, ARI, agreement, sign agreement, dot")
        print("  extras     on - sigma2, ||sub-full||_2, numerical ranks per p")
        print("  output     one run directory with every source in it")
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
    spec.p_min = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                 default=str(spec.p_min)))
    spec.p_points = int(get_input(
        f"p-grid points (log-spaced, {spec.p_min:g}..1)", default=str(spec.p_points)))
    spec.reps = int(get_input("Bootstrap reps per p", default=str(spec.reps)))
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
    """Screening coverage, verdicts and median imbalance, before anything is picked."""
    print_header("Screening status")
    cap = f"eta<={MAX_ETA:g}"
    if not sources:
        return
    width = max(8, max(len(s.name) for s in sources) + 1)
    # only operators some chosen source actually has verdicts for: a column pair of "-"
    # per unscreened operator is what pushed this panel past an 80-column terminal
    shown = [k for k in ALL_OPERATORS
             if any(any(f"valid_{k}" in v for v in verdicts_by[c.name].values())
                    for c in sources)] or ["S", "B"]
    head = f"{'source':<{width}}{'trees':>6}{'scrn':>6}"
    for k in shown:
        head += f"{ARM_OF[k] + ' edge':>9}{'eta':>7}"
    head += f"{'L+B':>5}{cap:>9}{'swept':>7}"
    lines = [head]
    for c in sources:
        v, ids = verdicts_by[c.name], c.ids
        done = [t for t in ids if t in v]
        line = f"{c.name:<{width}}{len(ids):>6}{len(done):>6}"
        for k in shown:
            n_ok = sum(1 for t in done if v[t].get(f"valid_{k}"))
            etas = [v[t][f"eta_{k}"] for t in done if f"eta_{k}" in v[t]]
            med = float(np.median(etas)) if etas else float("nan")
            line += (f"{n_ok:>9}{med:>7.1f}" if etas else f"{'-':>9}{'-':>7}")
        n_both = sum(1 for t in done if v[t].get("valid_S") and v[t].get("valid_B"))
        # the same operators the cap will actually weigh, so the panel cannot promise
        # a count the plan then contradicts
        cap_ops = gate_operators("both")
        n_cap = sum(1 for t in done
                    if max((v[t].get(f"eta_{k}", 0.0) for k in cap_ops),
                           default=0.0) <= MAX_ETA)
        swept = (len(list(sweep_cache_dir(c.name).glob("*.npz")))
                 if sweep_cache_dir(c.name).is_dir() else 0)
        lines.append(line + f"{n_both:>5}{n_cap:>9}{swept:>7}")
    lines += [
        "• <arm> edge  - trees whose split of the full matrix is a real edge of",
        "                the true tree; eta is the median imbalance of that split",
        "                (larger clan / smaller clan, so 1 is a perfectly even cut)",
        "• L+B         - valid under both L(S) and B, the pair the figure compares",
        f"• eta<={MAX_ETA:<7g}- even enough on L(S) and B to carry a recovery signal;",
        "                what a 'both' sweep starts from. A gate naming one",
        "                operator weighs only that one.",
        "• swept       - trees the recovery sweep has already covered",
        "• operators with no verdicts yet are not shown; screen to add them",
    ]
    _framed(lines, rule_after=len(sources))


def _run(sources, *, batch_hint: str) -> None:
    """Status, stage, one configuration, the plan, then the run. Source-agnostic."""
    verdicts_by = {c.name: verdicts(c.name) for c in sources}
    _print_status(sources, verdicts_by)

    stage = get_menu_choice(
        "Run:", ["screening - split the full matrix, check eta and validity",
                 "recovery sweep - sweep trees over p for NMI curve "
                 "(needs screening)"], default_index=0)
    is_screen = stage.startswith("screening")

    spec = _ask_config([c.name for c in sources], is_screen, verdicts_by)
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
    from src.utils.generated_prompts import (prompt_generated_plan, print_pool_table,
                                             resolve_pool_ids)

    print_header("Generated data")
    plan_g = prompt_generated_plan()
    if plan_g.pooled:
        # the eta pool is the only place a rejection-sampled tree exists; top it up
        # before anything asks for ids, because what comes back is what is on disk
        print_pool_table(plan_g)
        plan_g.tree_ids = resolve_pool_ids(plan_g)
        if not plan_g.tree_ids:
            print_error("the eta pool is empty for these parameters and nothing could "
                        "be built -- nothing to run")
            return
    sources = sources_from_plan(plan_g)
    if not sources:
        print_error("no trees to run")
        return
    _run(sources, batch_hint="")
