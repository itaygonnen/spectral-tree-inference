#!/usr/bin/env python3
"""Read p_hat_star out of the cached IQ-TREE + Sailfish (Yule) tree sweeps.

The twin of ``harvest_pool_sweeps.py`` for the simulated datasets: every
``sailfish_yule_*_buckets`` source the sweep has touched holds one ``.npz`` per tree
with one NMI curve per swept arm and the eta the screen measured. Both sizes and the
alignment length come from the source name, so a row of the collapse figure can
hold L fixed.

    python scripts/harvest_sailfish_sweeps.py
    python scripts/harvest_sailfish_sweeps.py --arm Lsym
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.utils.threshold_utils import find_discrete_threshold  # noqa: E402

CACHE = ROOT.parent / "results" / "real_data" / "_cache"
OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid" / "sailfish_nmi.csv"
LSYM_OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid" / "lsym_nmi.csv"

SOURCE = re.compile(r"^sailfish_yule_n(?P<m>\d+)_l(?P<L>\d+)_g.+_buckets$")
TREE = re.compile(r"^(?P<bin>eta[0-9p]+)_(?P<idx>\d+)$")


def harvest_tree(path: Path, threshold: float, arm: str = "L") -> dict | None:
    """One tree's p_hat_star on the measured points, or None if it never crosses.

    ``arm`` picks the operator, and every field it reads is keyed by it: a sweep
    run with ``--operators Lsym`` carries no L columns at all, so a function that
    hardcodes the arm silently returns the wrong population or nothing.
    """
    z = np.load(path, allow_pickle=True)
    key = f"nmi_{arm}_vs_fullmatrix"
    if key not in z.files:
        return None
    p = np.asarray(z["p_values"], float)
    nmi = np.asarray(z[key], float)
    ok = np.isfinite(nmi)
    if f"measured_{arm}" in z.files:
        ok &= np.asarray(z[f"measured_{arm}"], float) > 0
    ps = find_discrete_threshold(p[ok], nmi[ok], threshold)
    if ps is None:
        return None
    return {
        "eta": float(np.asarray(z[f"eta_ref_{arm}"]).ravel()[0]),
        "p_star": float(ps),
        # Pinned at an end of the grid is a bound, not a threshold.
        "censored": bool(ps <= p[ok].min() * 1.001 or ps >= p[ok].max() * 0.999),
    }


def _L_from(name: str) -> int:
    """Alignment length out of the dataset directory name, e.g. ``..._l1000_...``."""
    m = re.search(r"_l(\d+)_", name)
    return int(m.group(1)) if m else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", type=float, default=0.90)
    ap.add_argument("--arm", default="L", choices=("L", "Lsym"),
                    help="which operator's curve to read off (default: L)")
    args = ap.parse_args()

    rows, missed = [], 0
    for src in sorted(CACHE.glob("sailfish_yule_*_buckets")):
        s = SOURCE.match(src.name)
        if not s:
            continue
        for f in sorted((src / "sweep").glob("*.npz")):
            t = TREE.match(f.stem)
            row = harvest_tree(f, args.threshold, args.arm) if t else None
            if row is None:
                missed += 1
                continue
            rows.append({"source": "sailfish" if args.arm == "L" else "lsym",
                         "m": int(s["m"]), "L": int(s["L"]),
                         "eta_bin": t["bin"], "sample": int(t["idx"]), **row})

    df = pd.DataFrame(rows)
    out = OUT if args.arm == "L" else LSYM_OUT
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"{len(df)} trees ({missed} never crossed {args.threshold} or unreadable)")
    print(df.groupby(["L", "m"]).size().rename("trees").to_string())
    print(f"  eta      : {df.eta.min():.2f} .. {df.eta.max():.2f}")
    print(f"  censored : {int(df.censored.sum())}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
