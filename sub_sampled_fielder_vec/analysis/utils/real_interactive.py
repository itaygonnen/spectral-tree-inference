"""Interactive real-data menu for ``scripts/interactive_run.py``.

The launcher's other branch configures a *simulated* sweep (tree model, mutation rate,
sequence length). Real cohorts have none of those knobs -- the alignment and the true
tree are given -- so this asks a different, shorter set of questions:

    which cohorts  any <name>/{fasta,newick} under data/cohorts/, several at a time
    which stage    screen (eta + validity per operator) or sweep (recovery vs p)
    parameters     asked ONCE and applied to every chosen cohort, so sizes stay
                   comparable -- the whole point of running m=1000 beside m=6000

Everything it runs is the same code path as ``scripts/run_real_sweep.py``, so an
interactive session on a login node and a nohup'd batch run share one cache.

Lives beside the cohort helpers in ``analysis/utils/`` -- ``src/`` must never import
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

from .real_cohorts import (cohort_results_dir, list_cohorts,  # noqa: E402
                           log_path, screen_cache_path, sweep_cache_dir)
from .real_eta_screen import run_real_eta_screen             # noqa: E402
from .real_recovery_sweep import run_sweep                   # noqa: E402
from .real_results import Tee, export_all                    # noqa: E402

# Defaults, shared by every cohort in one run. The p-grid and reps match the notebook's
# figure, so a run left on defaults extends the caches the notebook plots from. p starts
# at 1e-4: the transition sits near log n / n, which is 1.4e-3 at m=6000, so the grid has
# to reach below that for the curve to show a floor rather than start on the ramp.
P_MIN, P_POINTS, REPS = 1e-4, 20, 10
ALL_TREES = "every tree, valid or not"

# Which screened trees a sweep may use. The gate is whether that operator's partition of
# the FULL matrix is a real single-edge split of the true tree -- sub-sampling recovery
# towards a reference that is not a tree edge measures stability, not correctness.
RULES = {
    "L(S) and B both cut a real tree edge": lambda r: r.get("valid_S") and r.get("valid_B"),
    "L(S) cuts a real tree edge": lambda r: r.get("valid_S"),
    "B cuts a real tree edge": lambda r: r.get("valid_B"),
    ALL_TREES: lambda r: True,
}


def _verdicts(cohort) -> dict:
    """Screen rows for a cohort, keyed by tree id. Empty when it has not been screened."""
    path = screen_cache_path(cohort.name)
    if not path.exists():
        return {}
    return {r["tree"]: r for r in np.load(path, allow_pickle=True)["rows"]
            if "error" not in r}


def _select_ids(ids, verdicts: dict, rule: str) -> list:
    if not verdicts:
        return list(ids)
    keep = RULES[rule]
    return [t for t in ids if t in verdicts and keep(verdicts[t])]


def _default_rule(cohorts, verdicts_by: Dict[str, dict]) -> str:
    """Strictest gate that still selects a tree in EVERY chosen cohort."""
    for name, fn in RULES.items():
        if all(any(t in verdicts_by[c.name] and fn(verdicts_by[c.name][t])
                   for t in c.ids()) for c in cohorts):
            return name
    return ALL_TREES


def _ask_config(cohorts, is_screen: bool, verdicts_by: Dict[str, dict]) -> dict:
    """One parameter set for all cohorts: show the defaults, edit them only on request."""
    cfg = dict(max_trees=0,                       # 0 = every tree in each cohort
               workers=max(1, min(8, (os.cpu_count() or 4) // 2)),
               rule=_default_rule(cohorts, verdicts_by),
               p_min=P_MIN, p_points=P_POINTS, reps=REPS)

    print_divider()
    print("configuration (applies to every cohort chosen):")
    if is_screen:
        print("  trees      all")
        print(f"  workers    {cfg['workers']}")
        print("  arms       L(S)+k-means and B=HDH+sign, min_split=5")
    else:
        print(f"  trees      all, gated on [{cfg['rule']}]")
        print(f"  p-grid     {cfg['p_points']} log-spaced points, {cfg['p_min']:g} .. 1.0")
        print(f"  reps       {cfg['reps']} bootstrap replicates per p")
        print("  metrics    NMI, ARI, agreement, sign agreement, dot "
              "(NMI is what the figure plots)")
    print()
    if not confirm("Edit this configuration?", default=False):
        return cfg

    cap = get_input("Max trees per cohort (blank = all)", default="")
    cfg["max_trees"] = int(cap) if cap and cap.strip() else 0
    if is_screen:
        cfg["workers"] = int(get_input("Workers", default=str(cfg["workers"])))
        return cfg

    # Spell the counts out. A bare "2, 91" reads as a range, and even "2/66" hides that
    # the denominator a gate can act on is the SCREENED trees, not the cohort.
    print()
    print("  Step 1 has to run before the sweep: it is what tells us, for each tree,")
    print("  whether an operator's split of the full matrix is a real edge of the")
    print("  true tree. Step 1 has been run on:")
    for c in cohorts:
        n_scr = sum(1 for t in c.ids() if t in verdicts_by[c.name])
        print(f"    {c.name}: {n_scr} of {len(c.ids())} trees")
    print("  The sweep can only use those trees; each option keeps this many of them:")
    labels = [
        f"{name}  ->  "
        + ", ".join(
            f"{c.name}: {len(_select_ids(c.ids(), verdicts_by[c.name], name))} tree(s)"
            for c in cohorts)
        for name in RULES]
    chosen = get_menu_choice("Reference partition must be a real tree edge under:",
                             labels, default_index=list(RULES).index(cfg["rule"]))
    cfg["rule"] = list(RULES)[labels.index(chosen)]
    cfg["p_min"] = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                   default=str(cfg["p_min"])))
    cfg["p_points"] = int(get_input(
        f"p-grid points (log-spaced, {cfg['p_min']:g}..1)", default=str(cfg["p_points"])))
    cfg["reps"] = int(get_input("Bootstrap reps per p", default=str(cfg["reps"])))
    return cfg


def _plan(cohorts, cfg: dict, is_screen: bool, verdicts_by: Dict[str, dict]) -> List[dict]:
    """Per-cohort tree lists and cost, under the shared configuration."""
    rows = []
    for c in cohorts:
        m, _ = c.shape()
        ids = c.ids(cfg["max_trees"] or None)
        if is_screen:
            secs = len(ids) * (75.0 if m >= 6000 else 10.0) / max(1, cfg["workers"])
        else:
            ids = _select_ids(ids, verdicts_by[c.name], cfg["rule"])
            # one sub-sampled Fiedler solve is ~9 s at m=6000 and scales as O(m^3)
            secs = len(ids) * cfg["p_points"] * cfg["reps"] * 9.0 * (m / 6000.0) ** 3
        rows.append(dict(cohort=c, m=m, ids=ids, hours=secs / 3600.0,
                         screened=sum(1 for t in c.ids() if t in verdicts_by[c.name])))
    return rows


def _print_status(cohorts, verdicts_by: Dict[str, dict]) -> None:
    """Where step 1 stands on each chosen cohort, before the stage is picked."""
    print_header("Step 1 status")
    for c in cohorts:
        v = verdicts_by[c.name]
        ids = c.ids()
        done = [t for t in ids if t in v]
        if not done:
            print(f"  {c.name}: not started ({len(ids)} trees)")
            continue
        n_s = sum(1 for t in done if v[t].get("valid_S"))
        n_b = sum(1 for t in done if v[t].get("valid_B"))
        n_both = sum(1 for t in done
                     if v[t].get("valid_S") and v[t].get("valid_B"))
        state = "complete" if len(done) == len(ids) else f"{len(done)} of {len(ids)} trees"
        print(f"  {c.name}: {state} - real tree edge under L(S) on {n_s}, "
              f"B on {n_b}, both on {n_both}")
        n_swept = len(list(sweep_cache_dir(c.name).glob("*.npz"))) \
            if sweep_cache_dir(c.name).is_dir() else 0
        if n_swept:
            print(f"  {' ' * len(c.name)}  step 2 done on {n_swept} tree(s)")
    print()


def run_real_data_menu() -> None:
    """Pick cohorts and a stage, configure once, then run the stage on each."""
    cohorts = list_cohorts()
    if not cohorts:
        print_error("No real cohorts found under data/cohorts/")
        print_warning("Expected <name>/fasta/*.fasta beside <name>/newick/*.nwk")
        return

    print_header("Real data")
    print("  Tip: select several with commas (e.g. '1,2') to run every size in turn")
    print()
    for i, c in enumerate(cohorts, 1):
        m, seq_len = c.shape()
        print_option(str(i), f"{c.name}  ({len(c.ids())} trees, m={m}, L={seq_len})",
                     highlight=(i == len(cohorts)))
    print()
    picks = get_multi_choice("Cohort(s)", [str(i) for i in range(1, len(cohorts) + 1)])
    chosen = [cohorts[int(i) - 1] for i in picks]

    verdicts_by = {c.name: _verdicts(c) for c in chosen}
    _print_status(chosen, verdicts_by)

    stage = get_menu_choice(
        "Stage:", ["step 1 - screen: which split each operator reads off the full "
                   "matrix, and is it a real tree edge",
                   "step 2 - sweep: how much of that split survives sub-sampling "
                   "(needs step 1)"], default_index=0)
    is_screen = stage.startswith("step 1")

    for c in chosen:
        n_screened = sum(1 for t in c.ids() if t in verdicts_by[c.name])
        if not is_screen and n_screened < len(c.ids()):
            print_warning(f"{c.name}: step 1 has only run on {n_screened} of "
                          f"{len(c.ids())} trees, and the sweep can only use trees "
                          "step 1 has covered -- run step 1 on this cohort first")

    cfg = _ask_config(chosen, is_screen, verdicts_by)
    plan = _plan(chosen, cfg, is_screen, verdicts_by)

    print_divider()
    total = sum(r["hours"] for r in plan)
    for r in plan:
        gate = "" if is_screen else f", gated on [{cfg['rule']}]"
        print(f"  {r['cohort'].name:<12} m={r['m']:<5} {len(r['ids']):>4} trees"
              f"  ~{r['hours']:.1f} h{gate}")
    print(f"\ntotal ~{total:.1f} h. Every tree is cached on its own, so this is safe to "
          "interrupt and resume.")
    if total > 1.0:
        print("for anything this long prefer:\n  nohup python scripts/run_real_sweep.py "
              f"--cohort \"{','.join(c.name for c in chosen)}\" "
              f"--stage {'screen' if is_screen else 'sweep'} "
              f"> logs/real_{'screen' if is_screen else 'sweep'}.log 2>&1 &")
    if not confirm("Run it here?", default=total <= 1.0):
        print_warning("Cancelled")
        return

    stage_name = "screen" if is_screen else "sweep"
    for k, r in enumerate(plan, 1):
        cohort, ids = r["cohort"], r["ids"]
        print_divider()
        print_header(f"[{k}/{len(plan)}] {cohort.name}")
        if not ids:
            print_warning(f"no tree passes [{cfg['rule']}] -- skipped")
            continue
        # an interactive run should leave the same record a nohup'd one does
        tee = Tee(log_path(cohort.name, stage_name))
        sys.stdout = tee
        try:
            _run_stage(cohort, ids, cfg, is_screen, r["m"])
        finally:
            tee.close()
        for f in export_all(cohort.name):
            print(f"  wrote {f}")
        print_success(f"results -> {cohort_results_dir(cohort.name)}")


def _run_stage(cohort, ids, cfg: dict, is_screen: bool, m: int) -> None:
    """One stage on one cohort, under the shared configuration."""
    if is_screen:
        run_real_eta_screen(ids, screen_cache_path(cohort.name),
                            cohort_name=cohort.name, workers=cfg["workers"])
        rows = _verdicts(cohort)
        n_s = sum(bool(v.get("valid_S")) for v in rows.values())
        n_b = sum(bool(v.get("valid_B")) for v in rows.values())
        print_success(f"{len(rows)} trees: L(S) cuts a real edge on {n_s}, B on {n_b}")
    else:
        p_values = np.logspace(np.log10(cfg["p_min"]), 0, cfg["p_points"])
        run_sweep(ids, sweep_cache_dir(cohort.name), p_values, reps=cfg["reps"],
                  cohort_name=cohort.name, m=m)
