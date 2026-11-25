#!/usr/bin/env python3
"""Where the L_sym grid stands: trees per size, per iso-``r`` cell.

Reads the screen caches rather than the datasets, so it reports what is actually
measurable: a tree only enters a cell once its eta is known, and only counts once
its L_sym split is valid. Batches of the same size (``..._b2``) are pooled, since
they are i.i.d. draws from one prior and differ only in when they were generated.

    python scripts/grid_status.py
    python scripts/grid_status.py --tol 0.25
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis.utils.r_grid import LSYM_TARGETS, bin_by_r, r_of  # noqa: E402

CACHE = ROOT.parent / "results" / "real_data" / "_cache"
#: i.i.d. Sailfish sets only: a ``_buckets`` set is rejection-sampled, so counting
#: it here would report the targets it was built to as if they were the prior.
IID = re.compile(r"^sailfish_yule_n(?P<m>\d+)_l(?P<L>\d+)_g[^/]*?(?<!_buckets)$")


def gather(seq_len: int = 1000) -> dict:
    """``{m: (screened, valid_eta, swept_stems)}`` pooled over batches."""
    out: dict = {}
    for d in sorted(CACHE.iterdir()):
        m_ = IID.match(d.name)
        if not m_ or int(m_["L"]) != seq_len:
            continue
        f = d / "screen.npz"
        if not f.exists():
            continue
        rows = np.load(f, allow_pickle=True)["rows"]
        eta = np.array([q.get("eta_Lsym", np.nan) for q in rows], float)
        ok = np.array([bool(q.get("valid_Lsym", False)) for q in rows]) & np.isfinite(eta)
        swept = {p.stem for p in (d / "sweep").glob("*.npz")} if (d / "sweep").is_dir() else set()
        n, v, s = out.get(int(m_["m"]), (0, [], 0))
        out[int(m_["m"])] = (n + len(rows), v + list(eta[ok]), s + len(swept))
    return {k: (n, np.array(v), s) for k, (n, v, s) in sorted(out.items())}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tol", type=float, default=0.18,
                    help="relative half-width of an r cell (default 0.18)")
    ap.add_argument("--seq-len", type=int, default=1000)
    args = ap.parse_args()

    data = gather(args.seq_len)
    if not data:
        print("no i.i.d. Sailfish screens found")
        return 1
    sizes = list(data)
    w = 8

    def row(label, cells, note=""):
        print(f"  {label:<14}" + "".join(f"{c:>{w}}" for c in cells) + f"  {note}")

    print(f"\nL_sym grid status   (cell = r within {args.tol:.0%}, L={args.seq_len})\n")
    print("  " + " " * 14 + "".join(f"{'m=' + str(m):>{w}}" for m in sizes))
    row("screened", [data[m][0] for m in sizes])
    row("valid L_sym", [len(data[m][1]) for m in sizes],
        f"{100 * sum(len(data[m][1]) for m in sizes) / max(1, sum(data[m][0] for m in sizes)):.0f}% overall")
    print()

    totals = {t: 0 for t in LSYM_TARGETS}
    unb = 0
    for t in LSYM_TARGETS:
        cells = []
        for m in sizes:
            eta = data[m][1]
            if not len(eta):
                cells.append(0)
                continue
            lab = bin_by_r(eta, np.full(len(eta), float(m)), tol=args.tol)
            c = int((lab == t).sum())
            cells.append(c)
            totals[t] += c
        # The eta a cell asks for, at the smallest and largest size present.
        lo, hi = (float(np.exp((np.log(t) - np.log(np.log(m) / m)) / 3) - 1)
                  for m in (sizes[0], sizes[-1]))
        row(f"r={t:<10.3f}", cells, f"eta {lo:.2f}..{hi:.2f}")
    print()
    for m in sizes:
        eta = data[m][1]
        if len(eta):
            lab = bin_by_r(eta, np.full(len(eta), float(m)), tol=args.tol)
            unb += int(np.isnan(lab).sum())
    row("unbinned", [int(np.isnan(bin_by_r(data[m][1],
                                           np.full(len(data[m][1]), float(m)),
                                           tol=args.tol)).sum())
                     if len(data[m][1]) else 0 for m in sizes],
        "matched no cell")
    row("swept", [data[m][2] for m in sizes], "L_sym curves on disk")
    print()
    print(f"  in a cell: {sum(totals.values())}    unbinned: {unb}    "
          f"empty cells: {sum(1 for t in LSYM_TARGETS for m in sizes if not _c(data, m, t, args.tol))}"
          f" of {len(LSYM_TARGETS) * len(sizes)}")
    return 0


def _c(data, m, t, tol) -> int:
    eta = data[m][1]
    if not len(eta):
        return 0
    lab = bin_by_r(eta, np.full(len(eta), float(m)), tol=tol)
    return int((lab == t).sum())


if __name__ == "__main__":
    raise SystemExit(main())
