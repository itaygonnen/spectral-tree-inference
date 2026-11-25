#!/usr/bin/env python3
"""Read p_hat_star out of the cached Kingman eta-pool sweeps.

The iso-``r`` grid asks whether cells sharing an ``r`` recover alike. The eta-pool
sweeps already on disk answer that over a far wider range than any grid we would
build now: eta in {1,5,10,15} crossed with n in {500..8000} spans
``r = (1+eta)^3 log m/m`` from 0.009 to about 51 -- four decades, already computed.

Each cached sweep gives the curve; the pool's own metadata gives the tree's TRUE
eta, which is what ``r`` must be computed from. The bin name is a target, not a
measurement, and at eta15 the realised value is nearer 14.4.

    python scripts/harvest_pool_sweeps.py
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.utils.threshold_utils import find_discrete_threshold  # noqa: E402

SWEEPS = ROOT / "cache" / "bootstrap_sweep"
POOL = ROOT / "cache" / "pool_sample"
OUT = ROOT / "analysis" / "notebooks_cache" / "r_grid" / "pool_nmi.csv"

#: The one configuration comparable to the CBM grid and the Sailfish
#: sweeps: sign cut, unnormalized Laplacian.
WANT = ("sign", None)

# <param_key>__eta<TT>__sample_<NNNN>__<hash>
NAME = re.compile(r"^(?P<param>.+?)__(?P<bin>eta\d+)__sample_(?P<idx>\d+)__[0-9a-f]+$")


def true_eta(param: str, bin_name: str, idx: str) -> float | None:
    """The tree's measured eta, from the pool metadata it was sampled into.

    Falls back to ``None`` rather than the bin's nominal value: a sweep whose pool
    entry is gone cannot be placed on the r axis, and guessing from the bin would
    put it at an r it was never measured at.
    """
    meta = POOL / param / bin_name / f"sample_{idx}" / "metadata.json"
    if not meta.exists():
        return None
    return float(json.loads(meta.read_text())["eta"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", type=float, default=0.90)
    args = ap.parse_args()

    rows, skipped = [], 0
    for d in sorted(SWEEPS.iterdir()):
        m = NAME.match(d.name)
        res = d / "result.json"
        if not m or not res.exists():
            skipped += 1
            continue
        eta = true_eta(m["param"], m["bin"], m["idx"])
        n = int(re.search(r"^n(\d+)_", m["param"]).group(1))
        if eta is None:
            skipped += 1
            continue

        blob = json.loads(res.read_text())
        cfg = blob.get("config", {})
        # The cache holds three different measurements side by side: kmeans on the
        # normalized Laplacian, the sign rule, and a sigma2 criterion. They are not
        # the same quantity, and the other two populations here are scored by the
        # sign rule on the unnormalized L = D - S, so only that subset is comparable.
        if (cfg.get("partition_method"), cfg.get("laplacian")) != WANT:
            skipped += 1
            continue
        j = blob["result"]
        p = np.asarray(j["p_values"], float)
        nmi = np.asarray(j["partition_nmi_M"], float)
        ok = np.isfinite(nmi)
        ps = find_discrete_threshold(p[ok], nmi[ok], args.threshold)
        if ps is None:
            continue
        rows.append({
            "source": "kingman_pool", "eta": eta, "m": n,
            "eta_bin": m["bin"], "sample": int(m["idx"]),
            "p_star": float(ps),
            # Pinned at an end of the grid is a bound, not a threshold.
            "censored": bool(ps <= p[ok].min() * 1.001 or ps >= p[ok].max() * 0.999),
        })

    df = pd.DataFrame(rows)
    # One sweep per (tree, p-grid); the same tree can appear under several hashes
    # (different p grids or stop settings), so keep the one that measured most.
    df = (df.sort_values("censored")
            .drop_duplicates(subset=["m", "eta_bin", "sample"], keep="first"))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    r = ((1.0 + df.eta) ** 3) * np.log(df.m) / df.m
    print(f"{len(df)} trees from {len(rows)} sweeps ({skipped} dirs skipped)")
    print(f"  m        : {sorted(df.m.unique())}")
    print(f"  eta      : {df.eta.min():.2f} .. {df.eta.max():.2f}")
    print(f"  r        : {r.min():.4f} .. {r.max():.1f}  ({np.log10(r.max()/r.min()):.1f} decades)")
    print(f"  p*       : {df.p_star.min():.5f} .. {df.p_star.max():.4f}")
    print(f"  censored : {int(df.censored.sum())}")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
