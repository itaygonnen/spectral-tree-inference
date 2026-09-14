#!/usr/bin/env python3
"""Non-interactive runner for the real-dataset screen and recovery sweep.

Two stages, both resumable, both driven by :mod:`analysis.utils.real_datasets`:

    screen   eta + validity per tree for L(S) (k-means or sign) and B=HDH+sign
             -> one .npz of rows per dataset
    sweep    NMI vs sub-sampling rate p for the same two arms
             -> one .npz per tree

This is what runs on a cluster: no prompts, everything on the command line, safe to
nohup and safe to kill. The interactive launcher (``scripts/interactive_run.py`` ->
"real data") asks the same questions in a menu and then calls the same
:func:`analysis.utils.real_run.execute`, so the two produce identical run directories.

    python scripts/run_real_sweep.py --list
    python scripts/run_real_sweep.py --dataset "1000 taxa,6000 taxa" --stage screen --workers 8
    nohup python scripts/run_real_sweep.py --dataset "6000 taxa" --stage sweep \
        > logs/real_sweep.log 2>&1 &

The sweep defaults to the trees whose partition is a real tree edge under BOTH operators
and whose reference split is no more lopsided than eta=20; ``--dataset-rule`` widens the
first gate when one arm never clears it, ``--max-eta 0`` removes the second.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]          # sub_sampled_fielder_vec
_REPO = _ROOT.parent
for _p in (str(_ROOT), str(_REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.utils.real_datasets import describe_search, list_datasets   # noqa: E402
from analysis.utils.real_run import (GATE_ALIASES, GATES, MAX_ETA,        # noqa: E402
                                     P_MIN, P_POINTS, REPS, RunSpec,
                                     execute, gate_label, plan)


def _list() -> None:
    datasets = list_datasets()
    for c in datasets:
        pr = c.probe()
        L = (f"L={pr['L_min']}" if pr["L_min"] == pr["L_max"]
             else f"L={pr['L_min']}..{pr['L_max']}  NOT ALIGNED")
        print(f"  {c.name!r}: {len(c.ids())} trees, m={pr['m']}, {L}  -> {c.dir}")
    if not datasets:
        print("no datasets found.\n" + describe_search())


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="show datasets on disk and exit")
    ap.add_argument("--dataset", "--cohort", dest="dataset", default="6000 taxa",
                    help="dataset directory name; comma-separate to run several sizes "
                         "in turn, or 'all' for every dataset on disk")
    ap.add_argument("--stage", choices=("screen", "sweep", "both"), default="both")
    ap.add_argument("--limit", type=int, default=0, help="first N trees only")
    ap.add_argument("--workers", type=int, default=4, help="screen workers")
    ap.add_argument("--p-min", type=float, default=P_MIN,
                    help="smallest sub-sampling rate; the grid always ends at 1.0")
    ap.add_argument("--p-points", type=int, default=P_POINTS,
                    help="log-spaced p in [p-min, 1]")
    ap.add_argument("--reps", type=int, default=REPS, help="bootstrap reps per p")
    ap.add_argument("--num-gaps", type=int, default=10)
    ap.add_argument("--min-split", type=int, default=5)
    ap.add_argument("--dataset-rule", "--cohort-rule", dest="dataset_rule",
                    choices=tuple(GATES) + tuple(GATE_ALIASES), default="both",
                    help="which screened trees enter the sweep: "
                         + "; ".join(f"{k} = {v[0]}" for k, v in GATES.items()))
    ap.add_argument("--max-eta", type=float, default=MAX_ETA,
                    help="drop trees whose reference split is more lopsided than this "
                         "(0 keeps all); a 1/999 split is a real edge but has no "
                         "recovery signal")
    ap.add_argument("--extra-metrics", dest="extra_metrics", action="store_true",
                    default=True, help=argparse.SUPPRESS)   # the default; kept explicit
    ap.add_argument("--no-extra-metrics", dest="extra_metrics", action="store_false",
                    help="skip sigma2, ||sub-full||_2 and the numerical ranks at every p "
                         "(the diagnostics the original experiment kept). They cost ~4%% "
                         "of a sweep at m=6000, so this is rarely worth it")
    ap.add_argument("--prefix", default="",
                    help="name this run: results land in runs/<timestamp>-<name>/")
    ap.add_argument("--display-mode", choices=("progress", "debug"), default="progress",
                    help="progress: bars on the terminal, full detail in experiment.log; "
                         "debug: every line on the terminal, no bars")
    return ap


def spec_from_args(args) -> RunSpec:
    """The parsed command line as a :class:`RunSpec` -- the same object the menu builds."""
    names = ([c.name for c in list_datasets()] if args.dataset.strip() == "all"
             else [n.strip() for n in args.dataset.split(",") if n.strip()])
    return RunSpec(datasets=names, stage=args.stage, max_trees=args.limit,
                   workers=args.workers, gate=args.dataset_rule, max_eta=args.max_eta,
                   p_min=args.p_min, p_points=args.p_points, reps=args.reps,
                   num_gaps=args.num_gaps, min_split=args.min_split,
                   extra_metrics=args.extra_metrics, prefix=args.prefix,
                   display_mode=args.display_mode)


def main() -> None:
    args = build_parser().parse_args()

    if args.list:
        _list()
        return

    spec = spec_from_args(args)
    rows = plan(spec)
    for r in rows:
        print(f"dataset {r.name!r}: {len(r.ids)} trees, m={r.m}  ~{r.hours:.1f} h",
              flush=True)
        if spec.runs_sweep:
            print(f"  {r.n_screened} screened -> {r.n_gate} [{gate_label(spec.gate)}] "
                  f"-> {len(r.selected)} with eta <= {spec.max_eta:g}", flush=True)
        if r.warning:
            print(f"  WARNING: {r.warning}", flush=True)
    sys.exit(execute(spec, rows))


if __name__ == "__main__":
    main()
