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


# Which screened trees a sweep may use. The gate is whether that operator's partition of
# the FULL matrix is a real single-edge split of the true tree -- sub-sampling recovery
# towards a reference that is not a tree edge measures stability, not correctness.
RULES = {
    "L(S) and B both cut a real tree edge": lambda r: r.get("valid_S") and r.get("valid_B"),
    "L(S) cuts a real tree edge": lambda r: r.get("valid_S"),
    "B cuts a real tree edge": lambda r: r.get("valid_B"),
    "every tree, valid or not": lambda r: True,
}


def _select_ids(ids, verdicts: dict, rule: str) -> list:
    if not verdicts:
        return list(ids)
    keep = RULES[rule]
    return [t for t in ids if t in verdicts and keep(verdicts[t])]


ALL_TREES = "every tree, valid or not"

# Defaults for a real-cohort run. The p-grid and reps match the notebook's figure, so a
# run left on defaults extends the same caches the notebook plots from.
P_MIN, P_POINTS, REPS = 0.01, 20, 10


def _default_rule(ids, verdicts: dict) -> str:
    """Strictest gate that still selects a tree: both arms, else one, else no gate."""
    for name, fn in RULES.items():
        if any(t in verdicts and fn(verdicts[t]) for t in ids):
            return name
    return ALL_TREES


def _defaults(cohort, ids, verdicts: dict, is_screen: bool) -> dict:
    d = dict(n_trees=len(ids), workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
             rule=_default_rule(ids, verdicts) if verdicts else ALL_TREES,
             p_min=P_MIN, p_points=P_POINTS, reps=REPS)
    if is_screen:
        d["summary"] = [f"trees      {d['n_trees']} (all)",
                        f"workers    {d['workers']}",
                        "arms       L(S)+k-means and B=HDH+sign, min_split=5"]
    else:
        # A tree can only be gated on if it was screened. Say so out loud: a cohort with
        # 66 trees and 2 screened rows otherwise reports "1 of 66" with no hint why.
        n_screened = sum(1 for t in ids if t in verdicts)
        n_sel = len(_select_ids(ids, verdicts, d["rule"])) if verdicts else len(ids)
        d["n_screened"] = n_screened
        d["summary"] = [f"trees      {n_sel} of {n_screened} screened "
                        f"({d['n_trees']} in cohort)  [{d['rule']}]",
                        f"p-grid     {d['p_points']} log-spaced points, "
                        f"{d['p_min']:g} .. 1.0",
                        f"reps       {d['reps']} bootstrap replicates per p",
                        "metrics    NMI, ARI, agreement, sign agreement, dot "
                        "(NMI is what the figure plots)"]
    return d


def run_real_data_menu() -> None:
    """Ask for a cohort and a stage, then run it. Returns when the stage finishes."""
    cohorts = list_cohorts()
    if not cohorts:
        print_error("No real cohorts found under data/real_datasets/Datasets/")
        print_warning("Expected <name>/fasta/*.fasta beside <name>/newick/*.nwk")
        return

    print_header("Real data")
    labels = [f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})"
              for c, (m, seq_len) in ((c, c.shape()) for c in cohorts)]
    choice = get_menu_choice("Cohort:", labels, default_index=len(labels) - 1)
    cohort = cohorts[labels.index(choice)]
    m, _seq_len = cohort.shape()
    all_ids = cohort.ids()

    stage = get_menu_choice(
        "Stage:", ["screen (eta + validity per operator)",
                   "sweep (recovery NMI vs p)"], default_index=0)

    screen_cache = screen_cache_path(cohort.name)
    verdicts = _screen_verdicts(screen_cache)
    is_screen = stage.startswith("screen")
    d = _defaults(cohort, all_ids, verdicts, is_screen)

    print_divider()
    if not is_screen and d.get("n_screened", 0) < len(all_ids):
        print_warning(f"only {d.get('n_screened', 0)} of {len(all_ids)} trees are "
                      "screened, and the sweep can only gate on screened trees -- "
                      "run the screen stage first to use the whole cohort")
    print("defaults:")
    for line in d["summary"]:
        print(f"  {line}")
    print()
    tune = confirm("Edit these defaults?", default=False)

    n_trees = int(get_input(f"How many trees (max {len(all_ids)})",
                            default=str(d["n_trees"]))) if tune else d["n_trees"]
    ids = all_ids[:max(1, min(n_trees, len(all_ids)))]

    if is_screen:
        workers = int(get_input("Workers", default=str(d["workers"]))) if tune \
            else d["workers"]
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
    if not verdicts:
        print_warning("No screen cache for this cohort -- run the screen stage first "
                      "to gate the sweep on valid partitions. Sweeping all trees.")
        rule = ALL_TREES
    elif tune:
        # show the size of each option, so a rule that selects nothing is visible up front
        n_screened = sum(1 for t in ids if t in verdicts)
        labels = [f"{name}  ({sum(1 for t in ids if t in verdicts and fn(verdicts[t]))} "
                  f"of {n_screened} screened)" for name, fn in RULES.items()]
        chosen = get_menu_choice("Reference partition must be a real tree edge under:",
                                 labels, default_index=list(RULES).index(d["rule"]))
        rule = list(RULES)[labels.index(chosen)]
    else:
        rule = d["rule"]
    ids = _select_ids(ids, verdicts, rule)
    if not ids:
        print_error(f"No tree passes {rule!r}. Pick a wider rule.")
        return

    if tune:
        p_min = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                default=str(d["p_min"])))
        p_points = int(get_input(f"p-grid points (log-spaced, {p_min:g}..1)",
                                 default=str(d["p_points"])))
        reps = int(get_input("Bootstrap reps per p", default=str(d["reps"])))
    else:
        p_min, p_points, reps = d["p_min"], d["p_points"], d["reps"]
    p_values = np.logspace(np.log10(p_min), 0, p_points)

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
