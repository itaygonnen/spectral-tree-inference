#!/usr/bin/env python3
"""Headless bootstrap p-sweep builder (populates the sweep cache).

Thin CLI over :mod:`src.utils.eta_pool_sweep`, which holds the single
source of truth for the "discover samples -> run bootstrap p-sweep ->
collect results" loop shared with the ``02_real_data_sweeps`` notebooks
(same p_values / bootstrap_reps / seed / num_gaps / early-stop / per-eta
caps), so results land on the same cache keys the notebook reads.

Default method is ``kmeans`` only -- the paper's main-text operator
(``PAPER_METHOD`` in ``src/utils/eta_pool_sweep.py``). ``sign`` and ``sigma2``
are omitted by default (``sigma2`` is the expensive operator; skip it for new
tree sizes) but remain available via ``--methods`` (``choices=METHOD_SPECS``)
for rebuilding the appendix three-operator comparison.

Usage:
    python scripts/build_sweeps.py --ns 3000 6000 --methods kmeans
    python scripts/build_sweeps.py --ns 3000 6000 --methods sign sigma2 kmeans
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# --- project root on sys.path (matches notebook bootstrap) ---
ROOT = Path(__file__).resolve().parent
while ROOT.parent != ROOT and not (ROOT / "setup.py").exists():
    ROOT = ROOT.parent
PROJECT_ROOT = ROOT / "sub_sampled_fielder_vec"
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.utils.eta_pool_sweep import (
    METHOD_SPECS, ETA_TARGETS, run_and_collect_sweeps,
)


def build_sweeps(ns, methods, etas=ETA_TARGETS, max_per_eta=None):
    counts = {"computed": 0, "cached": 0, "skipped": 0}

    def _on_sweep(name, n, eta, idx, was_cached):
        if was_cached is None:
            counts["skipped"] += 1
            print(f"    SKIP    {name:6s} idx={idx:04d} (no pool entry)", flush=True)
        elif was_cached:
            counts["cached"] += 1
        else:
            counts["computed"] += 1
            print(f"    COMPUTE {name:6s} idx={idx:04d}", flush=True)

    def _on_group(n, eta, n_idxs):
        print(f"[n={n:5d} eta={eta:2d}] {n_idxs} sample(s)", flush=True)

    run_and_collect_sweeps(
        ns=ns, eta_targets=etas, methods=methods,
        max_per_eta=max_per_eta, use_cache=True, verbose=False,
        on_sweep=_on_sweep, on_group=_on_group,
    )
    print(f"\ndone: computed={counts['computed']} cached={counts['cached']} "
          f"skipped={counts['skipped']}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ns", type=int, nargs="+", required=True)
    p.add_argument("--methods", nargs="+", default=["kmeans"],
                   choices=list(METHOD_SPECS))
    p.add_argument("--etas", type=int, nargs="+", default=ETA_TARGETS)
    p.add_argument("--max-per-eta", type=int, default=None,
                   help="override per-eta sample cap uniformly (default: ETA_SAMPLE_CAP)")
    args = p.parse_args()
    print(f"ns={args.ns} methods={args.methods} etas={args.etas} max_per_eta={args.max_per_eta}", flush=True)
    build_sweeps(args.ns, args.methods, args.etas, max_per_eta=args.max_per_eta)


if __name__ == "__main__":
    main()
