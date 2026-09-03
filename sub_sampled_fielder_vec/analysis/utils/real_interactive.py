"""Interactive real-data menu for ``scripts/interactive_run.py``.

The launcher's other branches configure a *simulated* sweep (tree model, mutation rate,
sequence length). Real cohorts have none of those knobs -- the alignment and the true
tree are given -- so this asks a different, shorter set of questions:

    which cohort   any data/real_datasets/Datasets/<name>/{fasta,newick} on this machine
    which stage    screen (eta + validity per operator) or sweep (recovery vs p)
    how many trees, how many workers, what p-grid

Everything it runs is the same code path as ``scripts/run_real_sweep.py``, so an
interactive session on a login node and a nohup'd batch run share one cache.

Lives beside the cohort helpers in ``analysis/utils/`` -- ``src/`` must never import
from ``analysis/``, so the dependency runs launcher -> analysis -> src, never back.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, print_divider, print_error, print_header,
    print_success, print_warning)

from .real_cohorts import (list_cohorts, screen_cache_path,  # noqa: E402
                           sweep_cache_dir)
from .real_eta_screen import run_real_eta_screen             # noqa: E402
from .real_recovery_sweep import run_sweep                   # noqa: E402


def _screen_verdicts(screen_cache: Path) -> dict:
    if not screen_cache.exists():
        return {}
    return {r["tree"]: r for r in np.load(screen_cache, allow_pickle=True)["rows"]
            if "error" not in r}


def _select_ids(ids, verdicts: dict, rule: str) -> list:
    if rule == "all" or not verdicts:
        return list(ids)
    keep = {
        "both valid": lambda r: r.get("valid_S") and r.get("valid_B"),
        "L(S) valid": lambda r: r.get("valid_S"),
        "B valid": lambda r: r.get("valid_B"),
    }[rule]
    return [t for t in ids if t in verdicts and keep(verdicts[t])]


def run_real_data_menu() -> None:
    """Ask for a cohort and a stage, then run it. Returns when the stage finishes."""
    cohorts = list_cohorts()
    if not cohorts:
        print_error("No real cohorts found under data/real_datasets/Datasets/")
        print_warning("Expected <name>/fasta/*.fasta beside <name>/newick/*.nwk")
        return

    print_header("Real-data cohorts")
    labels = []
    for c in cohorts:
        m, seq_len = c.shape()
        labels.append(f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})")
        print(f"  • {labels[-1]}")
    print()

    choice = get_menu_choice("Cohort:", labels, default_index=len(labels) - 1)
    cohort = cohorts[labels.index(choice)]
    m, _seq_len = cohort.shape()
    all_ids = cohort.ids()

    stage = get_menu_choice(
        "Stage:", ["screen (eta + validity per operator)",
                   "sweep (recovery NMI vs p)"], default_index=0)
    n_trees = int(get_input(f"How many trees (max {len(all_ids)})",
                            default=str(len(all_ids))))
    ids = all_ids[:max(1, min(n_trees, len(all_ids)))]

    screen_cache = screen_cache_path(cohort.name)
    print_divider()

    if stage.startswith("screen"):
        default_workers = str(max(1, min(8, (os.cpu_count() or 4) // 2)))
        workers = int(get_input("Workers", default=default_workers))
        print(f"~{'75 s' if m >= 6000 else '10 s'} per tree; cached and resumable.")
        if not confirm(f"Screen {len(ids)} trees of {cohort.name!r}?", default=True):
            print_warning("Cancelled")
            return
        run_real_eta_screen(ids, screen_cache, cohort_name=cohort.name,
                                 workers=workers)
        rows = _screen_verdicts(screen_cache)
        n_s = sum(bool(r.get("valid_S")) for r in rows.values())
        n_b = sum(bool(r.get("valid_B")) for r in rows.values())
        print_success(f"screened {len(rows)} trees: L(S) valid {n_s}, B valid {n_b} "
                      f"-> {screen_cache}")
        return

    # ---- sweep -------------------------------------------------------------
    verdicts = _screen_verdicts(screen_cache)
    if not verdicts:
        print_warning("No screen cache for this cohort -- run the screen stage first "
                      "to gate the sweep on valid partitions. Sweeping all trees.")
        rule = "all"
    else:
        rule = get_menu_choice(
            "Which trees:", ["both valid", "L(S) valid", "B valid", "all"],
            default_index=0)
    ids = _select_ids(ids, verdicts, rule)
    if not ids:
        print_error(f"No tree passes {rule!r}. Pick a wider rule.")
        return

    p_points = int(get_input("p-grid points (log-spaced, 0.01..1)", default="20"))
    reps = int(get_input("Bootstrap reps per p", default="10"))
    p_values = np.logspace(-2, 0, p_points)

    # one sub-sampled Fiedler solve is ~9 s at m=6000 and scales as O(m^3)
    per_solve = 9.0 * (m / 6000.0) ** 3
    hours = len(ids) * p_points * reps * per_solve / 3600.0
    print(f"\nestimated {hours:.1f} h for {len(ids)} trees "
          f"({p_points} p x {reps} reps, ~{per_solve:.1f} s per L solve)")
    print("one .npz per tree -- safe to interrupt and resume; for anything over an "
          "hour prefer:\n  nohup python scripts/run_real_sweep.py --cohort "
          f"{cohort.name!r} --stage sweep > logs/real_sweep.log 2>&1 &")
    if not confirm("Run it here anyway?", default=hours < 1.0):
        print_warning("Cancelled")
        return

    run_sweep(ids, sweep_cache_dir(cohort.name), p_values, reps=reps,
                   cohort_name=cohort.name, m=m)
    print_success(f"sweep cached -> {sweep_cache_dir(cohort.name)}")
