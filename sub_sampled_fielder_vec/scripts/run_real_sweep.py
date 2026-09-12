#!/usr/bin/env python3
"""Non-interactive runner for the real-dataset screen and recovery sweep.

Two stages, both resumable, both driven by :mod:`analysis.utils.real_datasets`:

    screen   eta + validity per tree for L(S)+k-means and B=HDH+sign  (one .npz of rows)
    sweep    NMI vs sub-sampling rate p for the same two arms         (one .npz per tree)

This is what runs on a cluster: no prompts, everything on the command line, safe to
nohup and safe to kill. The interactive launcher (``scripts/interactive_run.py`` ->
"real data") calls the same functions.

    python scripts/run_real_sweep.py --list
    python scripts/run_real_sweep.py --dataset "1000 taxa,6000 taxa" --stage screen --workers 8
    nohup python scripts/run_real_sweep.py --dataset "6000 taxa" --stage sweep \
        > logs/real_sweep.log 2>&1 &

The sweep dataset defaults to the trees whose partition is a real tree edge under BOTH
operators; ``--dataset-rule`` widens that when one arm never clears the gate.
"""
from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]          # sub_sampled_fielder_vec
_REPO = _ROOT.parent
for _p in (str(_ROOT), str(_REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.utils.real_datasets import (                       # noqa: E402
    describe_search, get_dataset, list_datasets, new_run_dir, screen_cache_path,
    sweep_cache_dir)
from analysis.utils.real_results import (Tee, export_run,        # noqa: E402
                                         write_config)
from src.utils.logging import (close_log_file, set_display_mode,  # noqa: E402
                               setup_log_file)
from analysis.utils.real_eta_screen import run_real_eta_screen  # noqa: E402
from analysis.utils.real_recovery_sweep import run_sweep        # noqa: E402

SELECTION_RULES = ("both", "valid_S", "valid_B", "any", "all")


def screen_rows(dataset: str) -> dict:
    path = screen_cache_path(dataset)
    if not path.exists():
        return {}
    return {r["tree"]: r for r in np.load(path, allow_pickle=True)["rows"]
            if "error" not in r}


def select_ids(dataset: str, ids, rule: str, max_eta: float = 0.0) -> list:
    """Filter ``ids`` by the screen's validity verdicts under ``rule``."""
    if rule == "all":
        return list(ids)
    rows = screen_rows(dataset)
    if not rows:
        print(f"WARNING: no screen cache for {dataset!r} -- run --stage screen first; "
              f"using all {len(ids)} ids")
        return list(ids)
    keep = {
        "both": lambda r: r.get("valid_S") and r.get("valid_B"),
        "valid_S": lambda r: r.get("valid_S"),
        "valid_B": lambda r: r.get("valid_B"),
        "any": lambda r: r.get("valid_S") or r.get("valid_B"),
    }[rule]
    out = [t for t in ids if t in rows and keep(rows[t])]
    n_gate = len(out)
    if max_eta and max_eta > 0:
        out = [t for t in out
               if max(rows[t].get("eta_S", 0.0), rows[t].get("eta_B", 0.0)) <= max_eta]
    print(f"dataset rule {rule!r}: {n_gate}/{len(ids)} trees"
          + (f", {len(out)} after eta<={max_eta:g}" if max_eta else ""))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show datasets on disk and exit")
    ap.add_argument("--dataset", "--cohort", dest="dataset", default="6000 taxa",
                    help="dataset directory name; comma-separate to run several sizes "
                         "in turn, or 'all' for every dataset on disk")
    ap.add_argument("--stage", choices=("screen", "sweep", "both"), default="both")
    ap.add_argument("--limit", type=int, default=0, help="first N trees only")
    ap.add_argument("--workers", type=int, default=4, help="screen workers")
    ap.add_argument("--p-min", type=float, default=1e-4,
                    help="smallest sub-sampling rate; the grid always ends at 1.0")
    ap.add_argument("--p-points", type=int, default=20, help="log-spaced p in [p-min, 1]")
    ap.add_argument("--reps", type=int, default=10, help="bootstrap reps per p")
    ap.add_argument("--num-gaps", type=int, default=10)
    ap.add_argument("--min-split", type=int, default=5)
    ap.add_argument("--dataset-rule", "--cohort-rule", dest="dataset_rule",
                    choices=SELECTION_RULES, default="both",
                    help="which screened trees enter the sweep")
    ap.add_argument("--max-eta", type=float, default=20.0,
                    help="drop trees whose reference split is more lopsided than this "
                         "(0 keeps all); a 1/999 split is a real edge but has no "
                         "recovery signal")
    ap.add_argument("--prefix", default="",
                    help="name this run: results land in runs/<timestamp>-<name>/")
    ap.add_argument("--display-mode", choices=("progress", "debug"), default="progress",
                    help="progress: bars on the terminal, full detail in experiment.log; "
                         "debug: every line on the terminal, no bars")
    args = ap.parse_args()

    if args.list:
        datasets = list_datasets()
        for c in datasets:
            pr = c.probe()
            L = (f"L={pr['L_min']}" if pr["L_min"] == pr["L_max"]
                 else f"L={pr['L_min']}..{pr['L_max']}  NOT ALIGNED")
            print(f"  {c.name!r}: {len(c.ids())} trees, m={pr['m']}, {L}  -> {c.dir}")
        if not datasets:
            print("no datasets found.\n" + describe_search())
        return

    names = ([c.name for c in list_datasets()] if args.dataset.strip() == "all"
             else [n.strip() for n in args.dataset.split(",") if n.strip()])

    # one directory for the whole run, every dataset inside it
    run_dir = new_run_dir(args.prefix)
    write_config(run_dir, vars(args), names, args.stage)
    set_display_mode(args.display_mode)
    setup_log_file(str(run_dir))
    print(f"run directory: {run_dir}", flush=True)

    tee = Tee(run_dir / "run.log")
    sys.stdout = tee
    selected, failures, status = {}, [], "completed"
    try:
        for i, name in enumerate(names, 1):
            if len(names) > 1:
                print(f"\n=== {name} ({i}/{len(names)})", flush=True)
            try:
                selected[name] = _run_cohort(get_dataset(name), args)
            except KeyboardInterrupt:
                print(f"{name}: interrupted -- finished trees are kept", flush=True)
                status = "interrupted"
                break
            except Exception:
                # keep the other datasets and the export; the cache holds what finished
                failures.append(name)
                print(f"ERROR: {name} failed:", flush=True)
                traceback.print_exc(file=sys.stdout)
    finally:
        tee.close()
        close_log_file()

    if failures:
        status = "failed"
    run_grid = (None if args.stage == "screen"
                else np.logspace(np.log10(args.p_min), 0, args.p_points))
    selection = (None if args.stage == "screen"
                 else {"rule": args.dataset_rule, "max_eta": args.max_eta,
                       "per_dataset": {n: {"selected": len(v)}
                                      for n, v in selected.items()}})
    for f in export_run(run_dir, selected, status, ", ".join(failures), run_grid,
                        selection):
        print(f"  {f.name}")
    if failures:
        print(f"ERROR: {len(failures)} dataset(s) failed: {', '.join(failures)} "
              f"-- traceback in {run_dir / 'run.log'}")
    if status != "completed":
        print(f"run marked {status!r} in summary.json -- re-run the same command to "
              "continue from the cache")
    print(f"everything for this run is in {run_dir}")
    if failures:
        sys.exit(1)


def _run_cohort(dataset, args) -> list:
    ids = dataset.ids(args.limit or None)
    m, seq_len = dataset.shape()
    print(f"dataset {dataset.name!r}: {len(ids)} trees, m={m}, L={seq_len}", flush=True)

    if args.stage in ("screen", "both"):
        run_real_eta_screen(ids, screen_cache_path(dataset.name),
                            dataset=dataset.name, workers=args.workers)

    if args.stage in ("sweep", "both"):
        sweep_ids = select_ids(dataset.name, ids, args.dataset_rule, args.max_eta)
        if not sweep_ids:
            print("nothing to sweep -- no tree passes the dataset rule "
                  f"{args.dataset_rule!r}. Try --dataset-rule valid_S or all.")
            return []
        cache_dir = sweep_cache_dir(dataset.name)
        p_values = np.logspace(np.log10(args.p_min), 0, args.p_points)
        print(f"sweep: {len(p_values)} p x {args.reps} reps over {len(sweep_ids)} trees "
              f"-> {cache_dir}", flush=True)
        run_sweep(sweep_ids, cache_dir, p_values, reps=args.reps,
                  num_gaps=args.num_gaps, min_split=args.min_split,
                  dataset=dataset.name, m=m)
        return sweep_ids
    return list(ids)


if __name__ == "__main__":
    main()
