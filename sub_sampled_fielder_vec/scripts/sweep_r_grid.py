#!/usr/bin/env python3
"""Sweep the synthetic flat-CBM cells of the iso-``r`` grid.

One cell per ``(r_target, eta)`` from :mod:`analysis.utils.r_grid`, swept over the
same ``p`` grid and trial count the paper's synthetic figure uses, scoring NMI --
the metric ``p_hat_star`` is read off.

The point of the grid is that cells sharing an ``r`` should share a ``p_hat_star``,
so the cells are swept **smallest first**: the small-``m`` end of every band is
cheap, and a partial run already shows whether the bands are flat.

Resumable. Every trial is cached individually under ``cache/sweep_trial``, so a
killed run is continued by re-running the identical command. Output is flushed per
line, since a backgrounded run's stdout is block-buffered otherwise.

    python scripts/sweep_r_grid.py --dry-run
    python scripts/sweep_r_grid.py --max-m 8500
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis.utils.block_model import build_flat_cbm_S           # noqa: E402
from analysis.utils.cache import get_or_compute_full, make_full_key  # noqa: E402
from analysis.utils.recovery import compute_nmi                   # noqa: E402
from analysis.utils.spectral import compute_fiedler_of_S          # noqa: E402
from analysis.utils.sweep import run_sweep                        # noqa: E402
from analysis.utils.r_grid import EXPONENT, R_TARGETS, synthetic_grid  # noqa: E402
from src.cache_io import CACHE_ROOT                               # noqa: E402

# Matched to analysis/paper/fig02_pstar_synth_cbm.ipynb so a cell of this grid and
# a cell of that figure are the same measurement and share the cache.
S_IN, S_OUT = 0.9, 0.05
P_GRID = np.geomspace(1e-4, 1.0, 25)
N_TRIALS = 8

OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid"


def log(msg: str) -> None:
    print(msg, flush=True)


def cached_trials(key: str) -> int:
    d = CACHE_ROOT / "sweep_trial" / key
    return len(list(d.glob("p*_seed*/.complete"))) if d.is_dir() else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-m", type=int, default=8500,
                    help="skip cells above this size (cost is cubic in m)")
    ap.add_argument("--r-targets", type=float, nargs="+", default=list(R_TARGETS))
    ap.add_argument("--dry-run", action="store_true",
                    help="list the cells and their cache state, sweep nothing")
    args = ap.parse_args()

    cells = [c for c in synthetic_grid(args.r_targets, etas=range(1, 9), exponent=EXPONENT)
             if c.m <= args.max_m]
    cells.sort(key=lambda c: c.m)

    log(f"{len(cells)} cells, m <= {args.max_m}, smallest first")
    log(f"p grid: {len(P_GRID)} points, {N_TRIALS} trials -> "
        f"{len(P_GRID) * N_TRIALS} solves per cell\n")

    for c in cells:
        key = make_full_key("flat_cbm", eta=c.eta, n=c.m, s_in=S_IN, s_out=S_OUT)
        log(f"  r={c.r_target:.2f}  eta={c.eta}  m={c.m:5d}  n1={c.n1:5d}  "
            f"cached trials={cached_trials(key):3d}  {key}")
    if args.dry_run:
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    t_all = time.time()

    for i, c in enumerate(cells, 1):
        key = make_full_key("flat_cbm", eta=c.eta, n=c.m, s_in=S_IN, s_out=S_OUT)
        t0 = time.time()

        def builder(c=c, key=key):
            S = build_flat_cbm_S(c.n1, c.m - c.n1, S_IN, S_OUT)
            return S, compute_fiedler_of_S(S), {
                "tree_model": "flat_cbm", "eta": c.eta, "n1": c.n1,
                "n2": c.m - c.n1, "m": c.m, "s_in": S_IN, "s_out": S_OUT,
                "full_key": key, "r_target": c.r_target, "r_actual": c.r_actual,
            }

        S, fiedler_full = get_or_compute_full(None, key, builder)

        # Guard the thing the grid depends on: eta realised == eta nominal. The
        # builder splits by integer division, so a bad m would silently move the
        # cell to another r, and the band would be comparing unlike cells.
        realised = max(c.n1, c.m - c.n1) / min(c.n1, c.m - c.n1)
        assert abs(realised - c.eta) < 1e-9, (
            f"eta mismatch at m={c.m}: nominal {c.eta}, realised {realised}"
        )

        df = run_sweep(None, key, S, fiedler_full, P_GRID, N_TRIALS,
                       metrics={"nmi": compute_nmi}, progress=False)
        df["source"] = "synthetic"
        df["eta"] = c.eta
        df["m"] = c.m
        df["r_target"] = c.r_target
        df["r_actual"] = c.r_actual
        frames.append(df)

        del S
        log(f"[{i}/{len(cells)}] r={c.r_target:.2f} eta={c.eta} m={c.m:5d} "
            f"done in {time.time() - t0:7.1f}s "
            f"(elapsed {(time.time() - t_all)/60:.1f} min)")

        # Written every cell, not at the end: the run is long and a killed job
        # should leave everything it finished behind it. CSV, not parquet: the
        # project's .venv has no parquet engine, and a format the reader cannot
        # open is not a checkpoint.
        # Written through a temp file and renamed: the plot scripts read this
        # path while the sweep is still going, and a half-written CSV is a
        # parse error rather than a short table.
        tmp = OUT / "sweep_nmi.csv.tmp"
        pd.concat(frames, ignore_index=True).to_csv(tmp, index=False)
        tmp.replace(OUT / "sweep_nmi.csv")

    log(f"\nwrote {OUT / 'sweep_nmi.csv'}  "
        f"({sum(len(f) for f in frames)} rows, {(time.time() - t_all)/60:.1f} min)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
