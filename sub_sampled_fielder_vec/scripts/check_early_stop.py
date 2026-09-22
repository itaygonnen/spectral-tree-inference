#!/usr/bin/env python3
"""Test the two assumptions the middle-out early stop rests on, against cached curves.

    python scripts/check_early_stop.py                # non-zero on any violation
    python scripts/check_early_stop.py --early-stop 3 --low 0.05

The stop skips p's it has not measured, so it is only sound if the curve really is
settled where it claims to be. Both claims are checked by replaying
``MiddleOut`` over every curve already on disk -- data that predates the code:

    HIGH   after ``early_stop`` consecutive NMI == 1.0, every larger p is also 1.0
    LOW    after ``early_stop`` consecutive NMI <= low, every smaller p is also <= low

A third check is structural rather than empirical: over a grid of ``(n_p, early_stop)``
and every possible pair of stop points, the visited indices and the two fill ranges must
partition ``range(n_p)`` exactly -- no index computed twice, none left neither computed
nor filled. That is the failure mode ``middle_out_runner`` has (see the note on
``MiddleOut``), and it is why this check exists rather than a comment.

Saving is reported two ways and neither invents a constant: the share of grid POINTS
skipped, and the share of DENSE points skipped (p > 0.1, where ``FiedlerVectorComputer``
switches from sparse ``eigsh`` to a dense ``eigh`` and the cost becomes O(m^3)). Arms are
weighted the way ``experiment_run.plan`` weights them -- a Griffing arm is
``GRIFFING_COST_SHARE`` of a Fiedler one.
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.runners.p_sweep_inner import HIGH_PERFECT, LOW_FLOOR, MiddleOut  # noqa: E402
from src.runners.experiment_run import GRIFFING_COST_SHARE  # noqa: E402
from src.runners.operators import ARM_OF, OPERATORS  # noqa: E402
from src.utils.run_paths import results_root  # noqa: E402

DENSE_ABOVE = 0.1        # src/core/fiedler_computer.py:65 -- dense eigh above this density
ARM_WEIGHT = {op.arm: (GRIFFING_COST_SHARE if op.kind == "griffing" else 1.0)
              for op in OPERATORS.values()}


def curves(path: Path) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """``(p_values, {arm: nmi curve})`` for one cached tree, full-matrix reference.

    Accepts both column spellings: the current ``nmi_<arm>_vs_fullmatrix`` and the
    pre-two-reference ``nmi_<arm>``, which the oldest caches still use.
    """
    z = np.load(path, allow_pickle=True)
    pv = np.asarray(z["p_values"], float)
    out: Dict[str, np.ndarray] = {}
    for arm in ARM_OF.values():
        for key in (f"nmi_{arm}_vs_fullmatrix", f"nmi_{arm}"):
            if key in z.files:
                a = np.asarray(z[key], float)
                if a.size == pv.size:
                    out[arm] = a
                break
    return pv, out


def replay(nmi: np.ndarray, early_stop: int, low: float,
           traversal: str = "middle_out") -> MiddleOut:
    """Drive ``MiddleOut`` with a curve that is already known, as an oracle."""
    trav = MiddleOut(len(nmi), early_stop=early_stop, high=HIGH_PERFECT, low=low,
                     traversal=traversal)
    for idx in trav.indices():
        trav.observe(idx, nmi[idx])
    return trav


def violations(nmi: np.ndarray, trav: MiddleOut, low: float) -> List[str]:
    """Where the stop claimed something the stored curve contradicts."""
    bad = []
    for i in trav.high_fill:
        if not nmi[i] >= HIGH_PERFECT:
            bad.append(f"high fill at p index {i}: stored NMI {nmi[i]:.4f} < 1.0")
    for i in trav.low_fill:
        if not nmi[i] <= low:
            bad.append(f"low fill at p index {i}: stored NMI {nmi[i]:.4f} > {low:g}")
    return bad


def structural_gaps(max_n_p: int = 40, max_stop: int = 4) -> List[str]:
    """Every (n_p, early_stop, where each side stops) must partition range(n_p) exactly.

    The oracle is a synthetic curve: 1.0 from ``h`` upwards, 0.0 from ``l`` downwards and
    0.5 in between, which fires the high side at ``h`` and the low side at ``l``
    independently. Sweeping h and l over the whole grid covers every stop position,
    including both sides firing in the same step and neither firing at all.
    """
    bad = []
    for n_p in range(1, max_n_p + 1):
        for stop in range(0, max_stop + 1):
            for h, l in itertools.product(range(n_p + 1), repeat=2):
                nmi = np.full(n_p, 0.5)
                nmi[h:] = 1.0
                nmi[:l] = 0.0
                for how in MiddleOut.TRAVERSALS:
                    trav = replay(nmi, stop, LOW_FLOOR, traversal=how)
                    tag = f"{how} n_p={n_p} stop={stop} h={h} l={l}"
                    seen = trav.measured
                    filled = set(trav.high_fill) | set(trav.low_fill)
                    if len(seen) != len(set(seen)):
                        bad.append(f"{tag}: index computed twice")
                    covered = set(seen) | filled
                    if sorted(covered) != list(range(n_p)):
                        bad.append(f"{tag}: neither computed nor filled: "
                                   f"{sorted(set(range(n_p)) - covered)}")
                    if set(seen) & filled:
                        bad.append(f"{tag}: filled AND computed: "
                                   f"{sorted(set(seen) & filled)}")
                    if 0 not in seen:
                        bad.append(f"{tag}: p_min not measured")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--early-stop", type=int, default=2,
                    help="consecutive settled points required per side (default 2)")
    ap.add_argument("--low", type=float, default=LOW_FLOOR,
                    help=f"low-side NMI floor (default {LOW_FLOOR:g})")
    ap.add_argument("--cache", default="", help="override the sweep cache root")
    args = ap.parse_args()

    root = Path(args.cache) if args.cache else results_root() / "_cache"
    files = sorted(root.glob("*/sweep/*.npz"))
    print(f"structural check: traversal covers every index exactly once ...", flush=True)
    gaps = structural_gaps()
    if gaps:
        for g in gaps[:20]:
            print(f"  FAIL {g}")
        print(f"  {len(gaps)} structural failure(s)")
    else:
        print("  ok -- no index computed twice, none left uncovered, p_min always measured")

    print(f"\nreplay over {len(files)} cached tree(s) under {root}"
          f"  (early_stop={args.early_stop}, low={args.low:g})")
    if not files:
        print("  no cached sweeps found -- nothing to replay")

    bad_total = 0
    saved_pts: List[Tuple[float, float]] = []      # (weight, fraction of points skipped)
    saved_dense: List[Tuple[float, float]] = []    # (weight, fraction of dense skipped)
    for f in files:
        pv, per_arm = curves(f)
        if not per_arm:
            continue
        n_dense = int(np.sum(pv > DENSE_ABOVE))
        print(f"\n  {f.parent.parent.name}/{f.stem}  n_p={len(pv)} "
              f"p={pv[0]:.3g}..{pv[-1]:.3g}")
        for arm, nmi in sorted(per_arm.items()):
            trav = replay(nmi, args.early_stop, args.low)
            bad = violations(nmi, trav, args.low)
            bad_total += len(bad)
            skipped = sorted(set(trav.high_fill) | set(trav.low_fill))
            dense_skipped = [i for i in skipped if pv[i] > DENSE_ABOVE]
            w = ARM_WEIGHT.get(arm, 1.0)
            saved_pts.append((w, len(skipped) / len(pv)))
            if n_dense:
                saved_dense.append((w, len(dense_skipped) / n_dense))
            hi = (f"p>={pv[trav.high_fill[0]]:.3g}" if trav.high_fill else "-")
            lo = (f"p<={pv[trav.low_fill[-1]]:.3g}" if trav.low_fill else "-")
            print(f"    {arm:<4} skips {len(skipped):>2}/{len(pv)} points "
                  f"({len(dense_skipped)}/{n_dense} dense)   high {hi:<12} low {lo}")
            for b in bad:
                print(f"      FAIL {arm}: {b}")

    def weighted(rows):
        if not rows:
            return float("nan")
        w = sum(r[0] for r in rows)
        return sum(r[0] * r[1] for r in rows) / w if w else float("nan")

    if saved_pts:
        print(f"\n  cost-weighted saving (arm weight: Fiedler 1.0, Griffing "
              f"{GRIFFING_COST_SHARE:g})")
        print(f"    grid points skipped : {100 * weighted(saved_pts):.1f}%")
        print(f"    dense points skipped: {100 * weighted(saved_dense):.1f}%   "
              f"(p > {DENSE_ABOVE:g}; this is where the O(m^3) solves are)")

    total = len(gaps) + bad_total
    print(f"\n{'FAILED' if total else 'OK'}: {total} violation(s)")
    return 1 if total else 0


if __name__ == "__main__":
    raise SystemExit(main())
