#!/usr/bin/env python3
"""Bucket occupancy per dataset, per tolerance -- and the invariants the capping rests on.

    python scripts/check_eta_buckets.py                  # every screened dataset
    python scripts/check_eta_buckets.py --tol 1 2 3      # compare tolerances
    python scripts/check_eta_buckets.py --arm B          # bucket on eta_B instead

Two jobs in one script, because they want the same table.

**The tool.** "How many trees do I actually get at tol=1, and which buckets come up
empty?" is the question that decides whether a run is worth launching, and it is
answerable from the screen caches alone -- no sweep, no run directory. Widening the
tolerance is the knob, so the table is printed for each tolerance asked for.

**The test.** Exits non-zero unless, for every dataset and tolerance: every kept tree's
eta lies inside its bucket's window, no bucket exceeds its cap, no unbucketed tree is
kept, and the kept list is exactly the concatenated prefixes of the buckets. Those four
are what the selection claims; a figure built on a bucket holding an out-of-window tree
would be mislabelled rather than merely noisy.

Occupancy is a property of the data, so the numbers move when a dataset is rescreened.
The invariants are properties of the code and must hold always.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.runners.eta_buckets import (                      # noqa: E402
    ETA_TARGETS, PER_BUCKET, UNBUCKETED, apply_caps, bucket_counts, label,
    panel_lines, tree_eta, window,
)
from src.runners.experiment_run import GATES, gate_key, gate_operators  # noqa: E402
from src.utils.run_paths import results_root, screen_cache_path         # noqa: E402


def datasets() -> List[str]:
    """Every source with a screen cache on disk."""
    root = results_root() / "_cache"
    if not root.is_dir():
        return []
    return sorted(d.name for d in root.iterdir()
                  if d.is_dir() and screen_cache_path(d.name).exists())


def screen_rows(name: str) -> Dict[str, dict]:
    path = screen_cache_path(name)
    rows = np.load(path, allow_pickle=True)["rows"]
    return {r["tree"]: r for r in rows if "error" not in r}


def check(ids, rows, cap_ops, tol, targets, caps) -> List[str]:
    """The four invariants. Returns a list of failures, empty when sound."""
    buckets = bucket_counts(ids, rows, cap_ops, targets, tol)
    kept, _short = apply_caps(buckets, targets, caps)
    cap_of = {int(t): int(c) for t, c in zip(targets, caps)}
    bad: List[str] = []

    unbucketed = set(buckets[UNBUCKETED])
    in_bucket = {t: set(v) for t, v in buckets.items() if t is not UNBUCKETED}

    for tree_id in kept:
        owner = [t for t, v in in_bucket.items() if tree_id in v]
        if not owner:
            bad.append(f"kept {tree_id} belongs to no bucket")
            continue
        lo, hi = window(owner[0], tol)
        e = tree_eta(rows[tree_id], cap_ops)
        if not (lo <= e <= hi):
            bad.append(f"kept {tree_id} has eta={e:.3f} outside {label(owner[0])} "
                       f"window [{lo:g}, {hi:g}]")
    if set(kept) & unbucketed:
        bad.append(f"{len(set(kept) & unbucketed)} unbucketed tree(s) were kept")
    for t, v in in_bucket.items():
        n = len(set(kept) & v)
        if n > cap_of.get(t, 0):
            bad.append(f"{label(t)} kept {n} > cap {cap_of.get(t, 0)}")
    expected = sorted(x for t in (int(t) for t in targets)
                      for x in buckets.get(t, [])[:cap_of.get(t, 0)])
    if kept != expected:
        bad.append("kept list is not the concatenated per-bucket prefixes")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tol", type=float, nargs="+", default=[1.0, 2.0, 3.0],
                    help="tolerances to tabulate (default 1 2 3)")
    ap.add_argument("--targets", type=int, nargs="+", default=list(ETA_TARGETS))
    ap.add_argument("--per-bucket", type=int, nargs="+", default=list(PER_BUCKET))
    ap.add_argument("--gate", default="valid_L", choices=sorted(GATES),
                    help="which operators' eta to bucket on (default valid_L = S only)")
    ap.add_argument("--arm", default="", help="override the gate and bucket on this "
                                              "operator key alone, e.g. S or B")
    args = ap.parse_args()

    if len(args.per_bucket) != len(args.targets):
        print(f"--per-bucket needs {len(args.targets)} values, got {len(args.per_bucket)}")
        return 2
    cap_ops = (args.arm,) if args.arm else gate_operators(gate_key(args.gate))

    names = datasets()
    print(f"bucketing on eta of {list(cap_ops)}  targets={args.targets}  "
          f"caps={args.per_bucket}")
    if not names:
        print(f"no screened dataset under {results_root() / '_cache'}")
    failures = 0
    for name in names:
        rows = screen_rows(name)
        # the gate's own validity predicate, so the table counts what a run would select
        keep = GATES[gate_key(args.gate)][1]
        ids = sorted(t for t, r in rows.items() if keep(r))
        etas = [tree_eta(rows[t], cap_ops) for t in ids]
        print(f"\n=== {name}: {len(rows)} screened, {len(ids)} pass [{args.gate}]")
        if etas:
            print(f"    eta  median={np.median(etas):.2f}  "
                  + "  ".join(f"p{q}={np.percentile(etas, q):.1f}"
                              for q in (10, 25, 50, 75, 90)))
        for tol in args.tol:
            buckets = bucket_counts(ids, rows, cap_ops, args.targets, tol)
            print(f"\n  tol = {tol:g}")
            for line in panel_lines(buckets, args.targets, args.per_bucket, tol):
                print("    " + line)
            for msg in check(ids, rows, cap_ops, tol, args.targets, args.per_bucket):
                print(f"    FAIL {msg}")
                failures += 1

    print(f"\n{'FAILED' if failures else 'OK'}: {failures} invariant violation(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
