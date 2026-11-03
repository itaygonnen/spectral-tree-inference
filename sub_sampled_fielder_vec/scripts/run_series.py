#!/usr/bin/env python3
"""Run a whole size series from one JSON: simulate -> screen -> sweep, per size.

    python scripts/run_series.py --config cluster/series.json
    python scripts/run_series.py --config cluster/series.json --only 8000
    python scripts/run_series.py --config cluster/series.json --dry-run

One file describes the job because the job has two halves that used to be configured
differently: the *simulation* (sizes, buckets, alignment length) took CLI flags while the
*sweep* took a ``config.json``. Splitting one run across two formats is how a series gets
launched with the caps you meant and the p-grid you did not.

Every stage is resumable, so this is safe to re-run: the fill continues from its
manifest's ``last_seed``, the screen reuses its cached rows, the sweep skips cached trees.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO = ROOT.parent
for _p in (str(ROOT), str(REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from src.models.sailfish_sim import SailfishSpec, preflight            # noqa: E402
from src.runners.experiment_run import RunSpec, execute, plan          # noqa: E402
from src.runners.sailfish_pool import fill_buckets                     # noqa: E402

SIM_KEYS = ("sizes", "seq_len", "eta_targets", "per_bucket", "eta_tol",
            "max_attempts", "iqtree_bin")


def load(path: Path) -> dict:
    cfg = json.loads(Path(path).read_text())
    missing = [k for k in ("sizes", "per_bucket", "eta_targets") if k not in cfg]
    if missing:
        raise SystemExit(f"{path}: missing {', '.join(missing)}")
    if len(cfg["per_bucket"]) == 1:
        cfg["per_bucket"] = cfg["per_bucket"] * len(cfg["eta_targets"])
    if len(cfg["per_bucket"]) != len(cfg["eta_targets"]):
        raise SystemExit(f"{path}: per_bucket and eta_targets are parallel lists")
    return cfg


def run_spec_from(cfg: dict) -> RunSpec:
    """The sweep half of the file, as a RunSpec. Unknown keys are the simulation's."""
    fields = {k: v for k, v in cfg.items()
              if not k.startswith("_") and k not in SIM_KEYS}
    fields.setdefault("eta_targets", cfg["eta_targets"])
    fields.setdefault("per_bucket", cfg["per_bucket"])
    fields.setdefault("eta_tol", cfg.get("eta_tol", 1.0))
    return RunSpec(**fields)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--only", type=int, action="append",
                    help="run just this size (repeatable); default is every size")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = load(args.config)
    sizes = args.only or cfg["sizes"]
    spec = run_spec_from(cfg)

    print(f"series: {sizes} taxa, L={cfg.get('seq_len', 1000)}, "
          f"{sum(cfg['per_bucket'])} trees each "
          f"({', '.join(f'{c}@eta{t}' for t, c in zip(cfg['eta_targets'], cfg['per_bucket']))})")
    print(f"sweep:  {spec.p_points} p in {spec.p_min:g}..1, {spec.reps} reps, "
          f"{spec.workers} workers, stage={spec.stage}")
    if args.dry_run:
        print("\n(dry run -- nothing written)")
        return 0

    missing = preflight(SailfishSpec(n_taxa=4, iqtree_bin=cfg.get("iqtree_bin")))
    if missing:
        print("\ncannot run -- missing dependencies:", file=sys.stderr)
        for line in missing:
            print(f"  - {line}", file=sys.stderr)
        return 2

    from analysis.utils.real_datasets import DATASETS_ROOT, get_dataset

    worst = 0
    for n in sizes:
        t0 = time.time()
        print(f"\n{'=' * 60}\n{n} taxa\n{'=' * 60}", flush=True)
        sim = SailfishSpec(n_taxa=int(n), root_seq_size=int(cfg.get("seq_len", 1000)),
                           pooled=True, iqtree_bin=cfg.get("iqtree_bin"))
        stats = fill_buckets(sim, DATASETS_ROOT, targets=cfg["eta_targets"],
                             per_bucket=cfg["per_bucket"],
                             tol=float(cfg.get("eta_tol", 1.0)),
                             workers=int(spec.workers),
                             max_attempts=int(cfg.get("max_attempts", 20000)))
        if not stats["complete"]:
            short = [f"eta{t:02d} {stats['counts'][int(t)]}/{c}"
                     for t, c in zip(cfg["eta_targets"], cfg["per_bucket"])
                     if stats["counts"][int(t)] < c]
            print(f"  WARNING: buckets short ({', '.join(short)}) -- sweeping what "
                  f"filled; re-run to continue from seed {stats['last_seed'] + 1}")
            worst = max(worst, 1)

        source = get_dataset(sim.name).source()
        rc = execute(spec, [source], plan(spec, [source]))
        worst = max(worst, rc)
        print(f"{n} taxa finished in {(time.time() - t0) / 3600:.2f} h", flush=True)
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
