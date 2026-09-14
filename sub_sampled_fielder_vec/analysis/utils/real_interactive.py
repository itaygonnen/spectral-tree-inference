"""Interactive real-data menu for ``scripts/interactive_run.py``.

The launcher's other branch configures a *simulated* sweep (tree model, mutation rate,
sequence length). Real datasets have none of those knobs -- the alignment and the true
tree are given -- so this asks a different, shorter set of questions:

    which datasets  any <name>/{fasta,newick} under data/tree_sets/, several at a time
    which stage    screening (eta + validity per operator) or recovery sweep (vs p)
    parameters     asked ONCE and applied to every chosen dataset, so sizes stay
                   comparable -- the whole point of running m=1000 beside m=6000

This module is only the questions and the status panel. Selecting trees, planning,
running and exporting all live in :mod:`analysis.utils.real_run`, which
``scripts/run_real_sweep.py`` calls too -- so a menu session on a login node and a
nohup'd batch run produce the same run directory from the same cache.

Lives beside the dataset helpers in ``analysis/utils/`` -- ``src/`` must never import
from ``analysis/``, so the dependency runs launcher -> analysis -> src, never back.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, get_multi_choice, print_divider, print_error,
    print_header, print_option, print_success, print_warning)

from .real_datasets import describe_search, list_datasets, sweep_cache_dir  # noqa: E402
from .real_run import (GATES, MAX_ETA, P_MIN, P_POINTS, REPS,  # noqa: E402
                       RunSpec, default_gate, execute, gate_label, plan, verdicts)


def _ask_config(datasets, is_screen: bool, verdicts_by: Dict[str, dict]) -> RunSpec:
    """One parameter set for all datasets: show the defaults, edit them only on request."""
    names = [c.name for c in datasets]
    spec = RunSpec(datasets=names,
                   stage="screen" if is_screen else "sweep",
                   max_trees=0,                   # 0 = every tree in each dataset
                   workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
                   gate=default_gate(names, verdicts_by),
                   max_eta=MAX_ETA, p_min=P_MIN, p_points=P_POINTS, reps=REPS,
                   extra_metrics=True, prefix="", display_mode="progress")

    print_divider()
    print("configuration (applies to every dataset chosen):")
    if is_screen:
        print("  trees      all")
        print(f"  workers    {spec.workers}")
        print("  arms       L(S)+k-means and B=HDH+sign, min_split=5")
    else:
        print(f"  trees      all, gated on [{gate_label(spec.gate)}], "
              f"eta <= {spec.max_eta:g}")
        print(f"  p-grid     {spec.p_points} log-spaced points, {spec.p_min:g} .. 1.0")
        print(f"  reps       {spec.reps} bootstrap replicates per p")
        print("  metrics    NMI, ARI, agreement, sign agreement, dot "
              "(NMI is what the figure plots)")
        print("  extras     on  - sigma2, ||sub-full||_2 and the numerical ranks,\n                       recorded at every p (~4% of the run at m=6000)")
        print("  output     one run directory with every dataset in it")
    print()
    if not confirm("Edit this configuration?", default=False):
        return spec

    cap = get_input("Max trees per dataset (blank = all)", default="")
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
    """Draw ``lines`` inside a box, so a status panel is not mistaken for a menu."""
    width = max(len(x) for x in lines)
    print("┌" + "─" * (width + 2) + "┐")
    for i, line in enumerate(lines):
        print(f"│ {line.ljust(width)} │")
        if i == rule_after:
            print("├" + "─" * (width + 2) + "┤")
    print("└" + "─" * (width + 2) + "┘")
    print()


def _print_status(datasets, verdicts_by: Dict[str, dict]) -> None:
    """Screening coverage, verdicts and median imbalance, before anything is picked."""
    print_header("Screening status")
    cap = f"eta<={MAX_ETA:g}"
    lines = [f"{'dataset':<12}{'trees':>7}{'screened':>10}{'L(S) edge':>11}"
             f"{'B edge':>8}{'both':>6}{'med eta_L':>11}{'med eta_B':>11}"
             f"{cap:>10}{'swept':>7}"]
    for c in datasets:
        v, ids = verdicts_by[c.name], c.ids()
        done = [t for t in ids if t in v]
        n_s = sum(1 for t in done if v[t].get("valid_S"))
        n_b = sum(1 for t in done if v[t].get("valid_B"))
        n_both = sum(1 for t in done if v[t].get("valid_S") and v[t].get("valid_B"))
        med_l = (float(np.median([v[t].get("eta_S", np.nan) for t in done]))
                 if done else float("nan"))
        med_b = (float(np.median([v[t].get("eta_B", np.nan) for t in done]))
                 if done else float("nan"))
        n_cap = sum(1 for t in done
                    if max(v[t].get("eta_S", 0.0), v[t].get("eta_B", 0.0)) <= MAX_ETA)
        swept = (len(list(sweep_cache_dir(c.name).glob("*.npz")))
                 if sweep_cache_dir(c.name).is_dir() else 0)
        lines.append(f"{c.name:<12}{len(ids):>7}{len(done):>10}{n_s:>11}"
                     f"{n_b:>8}{n_both:>6}{med_l:>11.1f}{med_b:>11.1f}"
                     f"{n_cap:>10}{swept:>7}")
    lines += [
        "• L(S) edge / B edge - trees whose split of the full matrix is a real edge",
        "                       of the true tree",
        "• med eta            - median imbalance of that split: larger clan / smaller",
        "                       clan, so 1 is a perfectly even cut",
        f"• eta<={MAX_ETA:<14g}- trees even enough on BOTH operators to carry a",
        "                       recovery signal; this is what the sweep starts from",
        "• swept              - trees the recovery sweep has already covered",
    ]
    _framed(lines, rule_after=len(datasets))


def run_real_data_menu() -> None:
    """Pick datasets and a stage, configure once, then run the stage on each."""
    datasets = list_datasets()
    if not datasets:
        print_error("No real datasets found.")
        print(describe_search())
        return

    print_header("Real data")
    print("  Tip: select several with commas (e.g. '1,2') to run every size in turn")
    print()
    for i, c in enumerate(datasets, 1):
        m, seq_len = c.shape()
        print_option(str(i), f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})",
                     highlight=(i == len(datasets)))
    print()
    picks = get_multi_choice("Dataset(s)", [str(i) for i in range(1, len(datasets) + 1)])
    chosen = [datasets[int(i) - 1] for i in picks]

    verdicts_by = {c.name: verdicts(c.name) for c in chosen}
    _print_status(chosen, verdicts_by)

    stage = get_menu_choice(
        "Run:", ["screening - split the full matrix, check eta and validity",
                 "recovery sweep - sweep trees over p for NMI curve "
                 "(needs screening)"], default_index=0)
    is_screen = stage.startswith("screening")

    spec = _ask_config(chosen, is_screen, verdicts_by)
    rows = plan(spec, verdicts_by)

    print_divider()
    total = sum(r.hours for r in rows)
    for r in rows:
        print(f"  {r.name:<12} m={r.m:<5} {len(r.selected):>4} trees  ~{r.hours:.1f} h")
        if not is_screen:
            # where the trees went: the two filters, in the order they are applied
            print(f"  {'':<12} {r.n_screened} screened -> {r.n_gate} "
                  f"[{gate_label(spec.gate)}] -> {len(r.selected)} with eta <= "
                  f"{spec.max_eta:g}")
        if r.warning:
            print_warning(f"{r.name}: {r.warning}")
    print(f"\ntotal ~{total:.1f} h. Every tree is cached on its own, so this is safe to "
          "interrupt and resume.")
    if total > 1.0:
        print("for anything this long prefer:\n  nohup python scripts/run_real_sweep.py "
              f"--dataset \"{','.join(c.name for c in chosen)}\" "
              f"--stage {spec.stage} "
              f"> logs/real_{spec.stage}.log 2>&1 &")
    if not confirm("Run it here?", default=total <= 1.0):
        print_warning("Cancelled")
        return

    def _announce(k, n, row):
        print_divider()
        print_header(f"[{k}/{n}] {row.name}")

    execute(spec, rows, on_dataset=_announce)
    print_success("done")
