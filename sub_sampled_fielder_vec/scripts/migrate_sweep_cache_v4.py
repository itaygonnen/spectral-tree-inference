#!/usr/bin/env python3
"""One-off: re-key the eta-pool sweep cache from schema v3 to v4, correcting the fill.

    python scripts/migrate_sweep_cache_v4.py --dry-run     # report, touch nothing
    python scripts/migrate_sweep_cache_v4.py               # rewrite alongside v3
    python scripts/migrate_sweep_cache_v4.py --drop-v3     # and delete the v3 entries

Why a migration and not a recompute
-----------------------------------
v3 entries were produced by an ascending scan that stopped after
``early_stop_consecutive_100`` consecutive p's at exactly 100% partition agreement and
then FABRICATED the tail: ``dot_product = 1.0`` and ``sign_agreement = 100.0``, neither
of which follows from a settled partition. v4 leaves both NaN.

The corrupt indices are not a guess. The old rule is deterministic in the stored
``partition_agreement_M`` curve and the stored ``config``, so this script reproduces the
exact fill boundary and nulls from it. The alternative -- bumping the schema and letting
503 entries miss -- is on the order of 60 hours of dense normalized-Laplacian solves
(52 entries at n=6000, 12 at n=8000) to arrive at values this computes exactly.

``partition_nmi_M`` / ``partition_ari_M`` / ``partition_agreement_M`` are untouched, and
those are the only fields ``eta_pool_sweep`` reads, so no figure moves.

v3 entries are left in place by default: keep them until ``--verify`` passes and the
figure has been rebuilt once.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from src.cache_io import SENTINEL, bootstrap_sweep as _scope  # noqa: E402
from src.utils.sweep_cache import compute_sweep_key           # noqa: E402

# The v3 fill trigger, verbatim: ``partition_agreement_M == 100.0`` exactly, N in a row,
# scanning ascending; the fill starts at the index AFTER the Nth.
V3_TRIGGER = 100.0
UNTRUSTED = ("dot_product", "sign_agreement")


def fill_start(agreement: List[float], consecutive: int) -> Optional[int]:
    """First index the v3 early stop filled, or None if it never fired."""
    if consecutive <= 0:
        return None
    run = 0
    for i, a in enumerate(agreement):
        if run >= consecutive:
            return i
        run = run + 1 if a == V3_TRIGGER else 0
    return None


def v4_key(config: Dict) -> Tuple[str, Dict]:
    """The v4 sweep key for a v3 config, with the field rename applied."""
    return compute_sweep_key(
        p_values=config["p_values"],
        bootstrap_reps=config["bootstrap_reps"],
        seed=config["seed"],
        num_gaps=config["num_gaps"],
        min_split=config["min_split"],
        early_stop=config["early_stop_consecutive_100"],
        partition_method=config["partition_method"],
        matrix_kind=config.get("matrix_kind", "similarity"),
        distance_alpha=config.get("distance_alpha", 1.0),
        laplacian=config.get("laplacian", "unnormalized"),
        # v3 had no traversal field because it had only one: walk the grid ascending
        # and stop at the top. Naming it is what keeps these entries distinct from a
        # middle-out sweep of the same grid.
        traversal="ascending",
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="report only")
    ap.add_argument("--drop-v3", action="store_true",
                    help="delete each v3 entry once its v4 twin is written")
    ap.add_argument("--verify", action="store_true",
                    help="check the v4 entries instead of writing them")
    args = ap.parse_args()

    root = _scope.root
    entries = sorted(d for d in root.iterdir()
                     if d.is_dir() and (d / "result.json").exists()) if root.exists() else []
    print(f"{len(entries)} cache entr{'y' if len(entries) == 1 else 'ies'} under {root}")

    n_v3 = n_v4 = n_written = n_nulled = n_skipped = n_bad = 0
    for d in entries:
        payload = json.loads((d / "result.json").read_text())
        config = payload.get("config", {})
        version = config.get("schema_version")
        if version == "v4":
            n_v4 += 1
            continue
        if version != "v3" or "early_stop_consecutive_100" not in config:
            n_skipped += 1
            continue
        n_v3 += 1

        result = payload["result"]
        start = fill_start(result["partition_agreement_M"],
                           int(config["early_stop_consecutive_100"]))
        key, new_config = v4_key(config)
        # the flat key is <param>__<eta bin>__sample_NNNN__<sweep key>; only the last
        # segment changes, so the entry keeps naming the same matrix
        target = root / (d.name.rsplit("__", 1)[0] + f"__{key}")

        if args.verify:
            if not (target / SENTINEL).exists():
                print(f"  MISSING v4 twin for {d.name}")
                n_bad += 1
                continue
            got = json.loads((target / "result.json").read_text())["result"]
            for k in ("partition_nmi_M", "partition_ari_M", "partition_agreement_M"):
                if got[k] != result[k]:
                    print(f"  CHANGED {k} in {target.name}")
                    n_bad += 1
            if start is not None:
                for k in UNTRUSTED:
                    if any(v is not None for v in got[k][start:]):
                        print(f"  NOT NULLED {k}[{start}:] in {target.name}")
                        n_bad += 1
            continue

        if start is not None:
            for k in UNTRUSTED:
                result[k] = result[k][:start] + [None] * (len(result[k]) - start)
            n_nulled += 1
        payload["config"] = new_config
        payload["migrated_from"] = d.name

        n_p = len(result["p_values"])
        note = (f"fill from index {start}/{n_p} (p={result['p_values'][start]:.4g}), "
                f"{n_p - start} point(s) nulled" if start is not None
                else "stop never fired, nothing to null")
        print(f"  {d.name[:52]:<52} -> ...{key}  {note}")

        if args.dry_run:
            continue
        target.mkdir(parents=True, exist_ok=True)
        sentinel = target / SENTINEL
        if sentinel.exists():
            sentinel.unlink()
        tmp = target / "result.json.tmp"
        tmp.write_text(json.dumps(payload, indent=2, default=str))
        tmp.replace(target / "result.json")
        sentinel.touch()
        n_written += 1
        if args.drop_v3 and target != d:
            shutil.rmtree(d)

    if args.verify:
        print(f"\nverify: {n_v4} v4 entr{'y' if n_v4 == 1 else 'ies'}, {n_bad} problem(s)")
        return 1 if n_bad else 0
    print(f"\nv3 {n_v3}, already v4 {n_v4}, unrecognised {n_skipped}")
    print(f"{'would null' if args.dry_run else 'nulled'} a fabricated tail in {n_nulled}"
          f"{'' if args.dry_run else f', wrote {n_written}'}")
    if args.dry_run:
        print("dry run -- nothing written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
