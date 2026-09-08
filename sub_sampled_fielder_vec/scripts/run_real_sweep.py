#!/usr/bin/env python3
"""Non-interactive runner for the real-cohort screen and recovery sweep.

Two stages, both resumable, both driven by :mod:`analysis.utils.real_cohorts`:

    screen   eta + validity per tree for L(S)+k-means and B=HDH+sign  (one .npz of rows)
    sweep    NMI vs sub-sampling rate p for the same two arms         (one .npz per tree)

This is what runs on a cluster: no prompts, everything on the command line, safe to
nohup and safe to kill. The interactive launcher (``scripts/interactive_run.py`` ->
"real data") calls the same functions.

    python scripts/run_real_sweep.py --list
    python scripts/run_real_sweep.py --cohort "1000 taxa,6000 taxa" --stage screen --workers 8
    nohup python scripts/run_real_sweep.py --cohort "6000 taxa" --stage sweep \
        > logs/real_sweep.log 2>&1 &

The sweep cohort defaults to the trees whose partition is a real tree edge under BOTH
operators; ``--cohort-rule`` widens that when one arm never clears the gate.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]          # sub_sampled_fielder_vec
_REPO = _ROOT.parent
for _p in (str(_ROOT), str(_REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.utils.real_cohorts import (                       # noqa: E402
    cohort_results_dir, get_cohort, list_cohorts, log_path, screen_cache_path,
    sweep_cache_dir)
from analysis.utils.real_results import Tee, export_all         # noqa: E402
from analysis.utils.real_eta_screen import run_real_eta_screen  # noqa: E402
from analysis.utils.real_recovery_sweep import run_sweep        # noqa: E402

COHORT_RULES = ("both", "valid_S", "valid_B", "any", "all")


def screen_rows(cohort_name: str) -> dict:
    path = screen_cache_path(cohort_name)
    if not path.exists():
        return {}
    return {r["tree"]: r for r in np.load(path, allow_pickle=True)["rows"]
            if "error" not in r}


def select_ids(cohort_name: str, ids, rule: str) -> list:
    """Filter ``ids`` by the screen's validity verdicts under ``rule``."""
    if rule == "all":
        return list(ids)
    rows = screen_rows(cohort_name)
    if not rows:
        print(f"WARNING: no screen cache for {cohort_name!r} -- run --stage screen first; "
              f"using all {len(ids)} ids")
        return list(ids)
    keep = {
        "both": lambda r: r.get("valid_S") and r.get("valid_B"),
        "valid_S": lambda r: r.get("valid_S"),
        "valid_B": lambda r: r.get("valid_B"),
        "any": lambda r: r.get("valid_S") or r.get("valid_B"),
    }[rule]
    out = [t for t in ids if t in rows and keep(rows[t])]
    print(f"cohort rule {rule!r}: {len(out)}/{len(ids)} trees")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show cohorts on disk and exit")
    ap.add_argument("--cohort", default="6000 taxa",
                    help="dataset directory name; comma-separate to run several sizes "
                         "in turn, or 'all' for every cohort on disk")
    ap.add_argument("--stage", choices=("screen", "sweep", "both"), default="both")
    ap.add_argument("--limit", type=int, default=0, help="first N trees only")
    ap.add_argument("--workers", type=int, default=4, help="screen workers")
    ap.add_argument("--p-min", type=float, default=1e-4,
                    help="smallest sub-sampling rate; the grid always ends at 1.0")
    ap.add_argument("--p-points", type=int, default=20, help="log-spaced p in [p-min, 1]")
    ap.add_argument("--reps", type=int, default=10, help="bootstrap reps per p")
    ap.add_argument("--num-gaps", type=int, default=10)
    ap.add_argument("--min-split", type=int, default=5)
    ap.add_argument("--cohort-rule", choices=COHORT_RULES, default="both",
                    help="which screened trees enter the sweep")
    ap.add_argument("--prefix", default="",
                    help="name this run: the sweep goes to <cohort>/sweep_<prefix>/ so "
                         "two grids or gates can sit side by side")
    args = ap.parse_args()

    if args.list:
        for c in list_cohorts():
            m, L = c.shape()
            print(f"  {c.name!r}: {len(c.ids())} trees, m={m}, L={L}  -> {c.dir}")
        return

    names = ([c.name for c in list_cohorts()] if args.cohort.strip() == "all"
             else [n.strip() for n in args.cohort.split(",") if n.strip()])
    for name in names:
        if len(names) > 1:
            print(f"\n=== {name} ({names.index(name) + 1}/{len(names)})", flush=True)
        _run_cohort(get_cohort(name), args)
    print("done.")


def _run_cohort(cohort, args) -> None:
    prefix = "" if args.stage == "screen" else args.prefix
    stage_log = log_path(cohort.name, "screen" if args.stage == "screen" else "sweep",
                         prefix)
    tee = Tee(stage_log)
    sys.stdout = tee
    try:
        _run_stages(cohort, args)
    finally:
        tee.close()
    for f in export_all(cohort.name, prefix):
        print(f"  wrote {f}")
    print(f"results -> {cohort_results_dir(cohort.name)}  (log: {stage_log.name})")


def _run_stages(cohort, args) -> None:
    ids = cohort.ids(args.limit or None)
    m, seq_len = cohort.shape()
    print(f"cohort {cohort.name!r}: {len(ids)} trees, m={m}, L={seq_len}", flush=True)

    if args.stage in ("screen", "both"):
        run_real_eta_screen(ids, screen_cache_path(cohort.name),
                            cohort_name=cohort.name, workers=args.workers)

    if args.stage in ("sweep", "both"):
        sweep_ids = select_ids(cohort.name, ids, args.cohort_rule)
        if not sweep_ids:
            print("nothing to sweep -- no tree passes the cohort rule "
                  f"{args.cohort_rule!r}. Try --cohort-rule valid_S or all.")
            return
        cache_dir = sweep_cache_dir(cohort.name, args.prefix)
        p_values = np.logspace(np.log10(args.p_min), 0, args.p_points)
        print(f"sweep: {len(p_values)} p x {args.reps} reps over {len(sweep_ids)} trees "
              f"-> {cache_dir}", flush=True)
        run_sweep(sweep_ids, cache_dir, p_values, reps=args.reps,
                  num_gaps=args.num_gaps, min_split=args.min_split,
                  cohort_name=cohort.name, m=m)


if __name__ == "__main__":
    main()
