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
import traceback
from pathlib import Path
from typing import Dict, List

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]      # sub_sampled_fielder_vec
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from src.utils.interactive_ui import (                       # noqa: E402
    confirm, get_input, get_menu_choice, get_multi_choice, print_divider, print_error,
    print_header, print_option, print_success, print_warning)

from .real_cohorts import (describe_search, list_cohorts,     # noqa: E402
                           new_run_dir, screen_cache_path, sweep_cache_dir)
from .real_eta_screen import run_real_eta_screen             # noqa: E402
from .real_recovery_sweep import run_sweep                   # noqa: E402
from .real_results import Tee, export_run, write_config      # noqa: E402
from src.utils.logging import (close_log_file, set_display_mode,  # noqa: E402
                               setup_log_file)

# Defaults, shared by every cohort in one run. The p-grid and reps match the notebook's
# figure, so a run left on defaults extends the caches the notebook plots from. p starts
# at 1e-4: the transition sits near log n / n, which is 1.4e-3 at m=6000, so the grid has
# to reach below that for the curve to show a floor rather than start on the ramp.
P_MIN, P_POINTS, REPS = 1e-4, 20, 10
ALL_TREES = "every tree, valid or not"

# Reference splits this lopsided are not a recovery question. At m=1000 the L(S) k-means
# cut routinely isolates ONE taxon (1/999, eta=999): a real pendant edge, so the validity
# gate passes it, but every sub-sample then picks a different singleton and NMI sits at 0
# for every p. 20 is well above the recoverability scale (m/log m)^(1/3) ~ 5.3 at m=1000,
# so it removes the degenerate cases without touching the ones the theory speaks about.
MAX_ETA = 20.0

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


def _select_ids(ids, verdicts: dict, rule: str, max_eta: float = 0.0) -> list:
    """Trees passing the validity gate, and (when ``max_eta`` > 0) not too lopsided."""
    if not verdicts:
        return list(ids)
    keep = RULES[rule]
    out = [t for t in ids if t in verdicts and keep(verdicts[t])]
    if max_eta and max_eta > 0:
        out = [t for t in out
               if max(verdicts[t].get("eta_S", 0.0),
                      verdicts[t].get("eta_B", 0.0)) <= max_eta]
    return out


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
               p_min=P_MIN, p_points=P_POINTS, reps=REPS, prefix="",
               display_mode="progress", max_eta=MAX_ETA)

    print_divider()
    print("configuration (applies to every cohort chosen):")
    if is_screen:
        print("  trees      all")
        print(f"  workers    {cfg['workers']}")
        print("  arms       L(S)+k-means and B=HDH+sign, min_split=5")
    else:
        print(f"  trees      all, gated on [{cfg['rule']}], "
              f"eta <= {cfg['max_eta']:g}")
        print(f"  p-grid     {cfg['p_points']} log-spaced points, {cfg['p_min']:g} .. 1.0")
        print(f"  reps       {cfg['reps']} bootstrap replicates per p")
        print("  metrics    NMI, ARI, agreement, sign agreement, dot "
              "(NMI is what the figure plots)")
        print("  output     <cohort>/sweep/  (name it below to keep runs side by side)")
    print()
    if not confirm("Edit this configuration?", default=False):
        return cfg

    cap = get_input("Max trees per cohort (blank = all)", default="")
    cfg["max_trees"] = int(cap) if cap and cap.strip() else 0
    if is_screen:
        # screening has no free parameters, so it has one canonical output per cohort
        cfg["workers"] = int(get_input("Workers", default=str(cfg["workers"])))
        return cfg

    # Spell the counts out. A bare "2, 91" reads as a range, and even "2/66" hides that
    # the denominator a gate can act on is the SCREENED trees, not the cohort.
    # No counts on the options: the Screening status table above already gives coverage
    # and verdicts, and the plan printed after this shows what the choice actually selects
    # (these labels were also counting the whole cohort, ignoring the cap asked for above).
    labels = list(RULES)
    chosen = get_menu_choice("Reference partition must be a real tree edge under:",
                             labels, default_index=list(RULES).index(cfg["rule"]))
    cfg["rule"] = list(RULES)[labels.index(chosen)]
    cfg["max_eta"] = float(get_input(
        "Drop trees whose reference split is more lopsided than eta (0 = keep all)",
        default=str(cfg["max_eta"])))
    cfg["p_min"] = float(get_input("Smallest sub-sampling rate p (largest is always 1.0)",
                                   default=str(cfg["p_min"])))
    cfg["p_points"] = int(get_input(
        f"p-grid points (log-spaced, {cfg['p_min']:g}..1)", default=str(cfg["p_points"])))
    cfg["reps"] = int(get_input("Bootstrap reps per p", default=str(cfg["reps"])))
    name = get_input("Name for this run (blank = timestamp only)", default="")
    cfg["prefix"] = (name or "").strip()
    cfg["display_mode"] = get_menu_choice(
        "Output:", ["progress - bars, and every line in the run's experiment.log",
                    "debug - every line on the terminal too, no bars"],
        default_index=0).split(" ")[0]
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
            after_rule = _select_ids(ids, verdicts_by[c.name], cfg["rule"])
            ids = _select_ids(ids, verdicts_by[c.name], cfg["rule"],
                              cfg.get("max_eta", 0.0))
            n_rule = len(after_rule)
            # one sub-sampled Fiedler solve is ~9 s at m=6000 and scales as O(m^3)
            secs = len(ids) * cfg["p_points"] * cfg["reps"] * 9.0 * (m / 6000.0) ** 3
        rows.append(dict(cohort=c, m=m, ids=ids, hours=secs / 3600.0,
                         n_rule=None if is_screen else n_rule,
                         screened=sum(1 for t in c.ids() if t in verdicts_by[c.name])))
    return rows


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


def _print_status(cohorts, verdicts_by: Dict[str, dict]) -> None:
    """Screening coverage, verdicts and median imbalance, before anything is picked."""
    print_header("Screening status")
    lines = [f"{'cohort':<12}{'trees':>7}{'screened':>10}{'L(S) edge':>11}"
             f"{'B edge':>8}{'both':>6}{'med eta_L':>11}{'med eta_B':>11}{'swept':>7}"]
    for c in cohorts:
        v, ids = verdicts_by[c.name], c.ids()
        done = [t for t in ids if t in v]
        n_s = sum(1 for t in done if v[t].get("valid_S"))
        n_b = sum(1 for t in done if v[t].get("valid_B"))
        n_both = sum(1 for t in done if v[t].get("valid_S") and v[t].get("valid_B"))
        med_l = (float(np.median([v[t].get("eta_S", np.nan) for t in done]))
                 if done else float("nan"))
        med_b = (float(np.median([v[t].get("eta_B", np.nan) for t in done]))
                 if done else float("nan"))
        swept = (len(list(sweep_cache_dir(c.name).glob("*.npz")))
                 if sweep_cache_dir(c.name).is_dir() else 0)
        lines.append(f"{c.name:<12}{len(ids):>7}{len(done):>10}{n_s:>11}"
                     f"{n_b:>8}{n_both:>6}{med_l:>11.1f}{med_b:>11.1f}{swept:>7}")
    lines += [
        "• L(S) edge / B edge - trees whose split of the full matrix is a real edge",
        "                       of the true tree",
        "• med eta            - median imbalance of that split: larger clan / smaller",
        "                       clan, so 1 is a perfectly even cut",
        "• swept              - trees the recovery sweep has already covered",
    ]
    _framed(lines, rule_after=len(cohorts))


def run_real_data_menu() -> None:
    """Pick cohorts and a stage, configure once, then run the stage on each."""
    cohorts = list_cohorts()
    if not cohorts:
        print_error("No real cohorts found.")
        print(describe_search())
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
        "Run:", ["screening - split the full matrix, check eta and validity",
                 "recovery sweep - sweep trees over p for NMI curve "
                 "(needs screening)"], default_index=0)
    is_screen = stage.startswith("screening")

    for c in chosen:
        n_screened = sum(1 for t in c.ids() if t in verdicts_by[c.name])
        if not is_screen and n_screened < len(c.ids()):
            print_warning(f"{c.name}: screening covers only {n_screened} of "
                          f"{len(c.ids())} trees, and the sweep can only use screened "
                          "trees -- screen this cohort first")

    cfg = _ask_config(chosen, is_screen, verdicts_by)
    plan = _plan(chosen, cfg, is_screen, verdicts_by)

    print_divider()
    total = sum(r["hours"] for r in plan)
    for r in plan:
        print(f"  {r['cohort'].name:<12} m={r['m']:<5} {len(r['ids']):>4} trees"
              f"  ~{r['hours']:.1f} h")
        if not is_screen:
            # where the trees went: the two filters, in the order they are applied
            print(f"  {'':<12} {r['screened']} screened -> {r['n_rule']} "
                  f"[{cfg['rule']}] -> {len(r['ids'])} with eta <= "
                  f"{cfg.get('max_eta', 0.0):g}")
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
    run_dir = new_run_dir(cfg["prefix"])
    write_config(run_dir, cfg, [c.name for c in chosen], stage_name)
    # progress bars on the terminal, every line in the run's own experiment.log
    set_display_mode(cfg.get("display_mode", "progress"))
    setup_log_file(str(run_dir))
    print(f"run directory: {run_dir}")

    # one log for the whole run, beside its results
    tee = Tee(run_dir / "run.log")
    sys.stdout = tee
    selected, failures, status = {}, [], "completed"
    try:
        for k, r in enumerate(plan, 1):
            cohort, ids = r["cohort"], r["ids"]
            print_divider()
            print_header(f"[{k}/{len(plan)}] {cohort.name}")
            if not ids:
                print_warning(f"no tree passes [{cfg['rule']}] -- skipped")
                continue
            try:
                _run_stage(cohort, ids, cfg, is_screen, r["m"])
                selected[cohort.name] = list(ids)
            except KeyboardInterrupt:
                print_warning(f"{cohort.name}: interrupted -- finished trees are kept")
                selected[cohort.name] = list(ids)
                status = "interrupted"
                break
            except Exception:
                # the per-tree cache holds whatever finished, and the other cohorts are
                # still worth running, so record the failure and carry on
                failures.append(cohort.name)
                print_error(f"{cohort.name} failed:")
                traceback.print_exc(file=sys.stdout)
                selected[cohort.name] = list(ids)
    finally:
        tee.close()
        close_log_file()

    if failures:
        status = "failed"
    run_grid = (None if is_screen
                else np.logspace(np.log10(cfg["p_min"]), 0, cfg["p_points"]))
    selection = (None if is_screen else
                 {"rule": cfg["rule"], "max_eta": cfg.get("max_eta", 0.0),
                  "per_cohort": {r["cohort"].name:
                                 {"screened": r["screened"], "after_rule": r["n_rule"],
                                  "after_max_eta": len(r["ids"])} for r in plan}})
    for f in export_run(run_dir, selected or {c.name: c.ids() for c in chosen}, status,
                        ", ".join(failures), run_grid, selection):
        print(f"  {f.name}")
    if failures:
        print_error(f"{len(failures)} cohort(s) failed: {', '.join(failures)} "
                    f"-- traceback in {run_dir / 'run.log'}")
    if status != "completed":
        print_warning(f"run marked {status!r} in summary.json -- re-run the same command "
                      "to continue from the cache")
    print_success(f"everything for this run is in {run_dir}")


def _run_stage(cohort, ids, cfg: dict, is_screen: bool, m: int) -> None:
    """One stage on one cohort, under the shared configuration."""
    if is_screen:
        run_real_eta_screen(ids, screen_cache_path(cohort.name),
                            cohort_name=cohort.name, workers=cfg["workers"])
        rows = _verdicts(cohort)
        if not rows:
            raise RuntimeError(
                f"{cohort.name}: screening produced nothing usable -- see the errors "
                f"above and failures.csv in the run directory")
        n_s = sum(bool(v.get("valid_S")) for v in rows.values())
        n_b = sum(bool(v.get("valid_B")) for v in rows.values())
        print_success(f"{len(rows)} trees: L(S) cuts a real edge on {n_s}, B on {n_b}")
    else:
        p_values = np.logspace(np.log10(cfg["p_min"]), 0, cfg["p_points"])
        run_sweep(ids, sweep_cache_dir(cohort.name), p_values,
                  reps=cfg["reps"], cohort_name=cohort.name, m=m)
