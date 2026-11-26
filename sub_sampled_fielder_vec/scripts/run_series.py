#!/usr/bin/env python3
"""Run a whole size series from one JSON: simulate -> screen -> sweep, per size.

    python scripts/run_series.py --config cluster/series.json
    python scripts/run_series.py --config cluster/series.json --only 8000
    python scripts/run_series.py --config cluster/series.json --dry-run
    python scripts/run_series.py --config cluster/series.json --gather

``"source": "generated"`` runs the in-process builders instead (Kingman and the rest of
``src/models/tree_models.py``): no IQ-TREE, no FASTA, trees from the eta pool, which is
reused where it is already filled and topped up where it is short. Everything after that
-- screen, gate, buckets, sweep, export -- is the same code on the same terms.

ONE run directory per invocation, holding every size it ran: the sizes are filled first
and then swept in a single ``execute``, so ``runs/<timestamp>/`` carries one config.json,
one screening.csv, one per_tree.csv, one curves.csv (its first column is the dataset) and
one recovery_curve.png across all of them. It used to call ``execute`` per size, which
scattered one experiment over five timestamps for no reason -- the exporter has always
taken a list of sources.

``--gather`` is the same export with nothing computed: it skips the fill and the sweep
finds every tree cached, so it collects sizes that were run SEPARATELY -- which is what
the SLURM array does, one task per size, each in its own process and so its own run
directory -- into one directory in seconds.

One file describes the job because the job has two halves that used to be configured
differently: the *simulation* (sizes, buckets, alignment length) took CLI flags while the
*sweep* took a ``config.json``. Splitting one run across two formats is how a series gets
launched with the caps you meant and the p-grid you did not.

Every stage is resumable, so this is safe to re-run: the fill continues from its
manifest's ``last_seed``, the screen reuses its cached rows, the sweep skips cached trees.

``operators`` and ``gate`` drive BOTH halves. The fill screens every operator listed and
accepts on the gate's arm, so a series with ``["S", "Lsym"]`` writes both verdicts at fill
time while still selecting the eta_S buckets the figures are stated in; the sweep then
reads both arms, reusing whatever is already cached per arm. They used to reach the sweep
alone, which is how a run meant to carry L_sym filled a pool that had never heard of it.
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
            "max_attempts", "iqtree_bin",
            # the generated (in-process) source's own keys
            "source", "models", "mutation_rate", "params")


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


def generated_sources(cfg: dict, build: bool = True):
    """Sources for an in-process series (Kingman and the other `tree_models` builders).

    No IQ-TREE, no FASTA: the trees come from the eta pool, which is rejection-sampled
    to the targets. ``resolve_pool_ids`` REUSES every sample already on disk and builds
    only what the cells are short of, so re-running this is how a pool grows. One source
    per (model, size) cell, which is what puts every size in one run directory.
    """
    from src.runners.generated_source import sources_from_plan
    from src.utils.generated_prompts import (GeneratedPlan, describe_plan,
                                             resolve_pool_ids)

    per = sorted({int(c) for c in cfg["per_bucket"]})
    if len(per) > 1:
        raise SystemExit("a generated series has ONE per-cell count: the eta pool is "
                         f"built per cell, so per_bucket must be uniform, not {per}")
    plan = GeneratedPlan(
        models=[str(m) for m in cfg.get("models", ["kingman"])],
        n_values=[int(n) for n in cfg["sizes"]],
        seq_len=int(cfg.get("seq_len", 10000)),
        mutation_rate=float(cfg.get("mutation_rate", 0.1)),
        params={k: float(v) for k, v in dict(cfg.get("params", {})).items()},
        etas=[int(t) for t in cfg["eta_targets"]],
        per_cell=per[0],
        max_attempts=int(cfg.get("max_attempts", 2000)),
    )
    print(f"  {describe_plan(plan)}", flush=True)
    plan.tree_ids = resolve_pool_ids(plan, build=build)
    if not plan.tree_ids:
        return []
    return sources_from_plan(plan)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--only", type=int, action="append",
                    help="run just this size (repeatable); default is every size")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--gather", action="store_true",
                    help="skip the fill and sweep only what is cached: collects sizes "
                         "run separately (the SLURM array) into one run directory")
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

    if cfg.get("source", "sailfish") == "generated":
        # no preflight: nothing external is called, the trees are built in process
        print(f"\n{'=' * 60}\ngenerated trees: reusing the eta pool, building the rest"
              f"\n{'=' * 60}", flush=True)
        sources = generated_sources(cfg, build=not args.gather)
        if not sources:
            print("\nthe eta pool is empty for these parameters", file=sys.stderr)
            return 2
        print(f"\n{'=' * 60}\nsweeping {len(sources)} size(s) into one run directory\n"
              f"{'=' * 60}", flush=True)
        t0 = time.time()
        rc = execute(spec, sources, plan(spec, sources))
        print(f"series finished in {(time.time() - t0) / 3600:.2f} h", flush=True)
        return rc

    missing = preflight(SailfishSpec(n_taxa=4, iqtree_bin=cfg.get("iqtree_bin")))
    if missing:
        print("\ncannot run -- missing dependencies:", file=sys.stderr)
        for line in missing:
            print(f"  - {line}", file=sys.stderr)
        return 2

    from analysis.utils.real_datasets import DATASETS_ROOT, get_dataset

    worst = 0
    sources = []
    for n in sizes:
        t0 = time.time()
        print(f"\n{'=' * 60}\n{n} taxa\n{'=' * 60}", flush=True)
        sim = SailfishSpec(n_taxa=int(n), root_seq_size=int(cfg.get("seq_len", 1000)),
                           pooled=True, iqtree_bin=cfg.get("iqtree_bin"))
        if args.gather:
            print("  --gather: reading what is cached, filling nothing", flush=True)
        else:
            stats = fill_buckets(sim, DATASETS_ROOT, targets=cfg["eta_targets"],
                                 per_bucket=cfg["per_bucket"],
                                 tol=float(cfg.get("eta_tol", 1.0)),
                                 workers=int(spec.workers),
                                 max_attempts=int(cfg.get("max_attempts", 20000)),
                                 operators=spec.operators, gate=spec.gate)
            if not stats["complete"]:
                short = [f"eta{t:02d} {stats['counts'][int(t)]}/{c}"
                         for t, c in zip(cfg["eta_targets"], cfg["per_bucket"])
                         if stats["counts"][int(t)] < c]
                print(f"  WARNING: buckets short ({', '.join(short)}) -- sweeping what "
                      f"filled; re-run to continue from seed {stats['last_seed'] + 1}")
                worst = max(worst, 1)
            print(f"{n} taxa filled in {(time.time() - t0) / 3600:.2f} h", flush=True)
        try:
            sources.append(get_dataset(sim.name).source())
        except Exception as exc:              # a size that was never generated
            print(f"  WARNING: {sim.name} not on disk ({exc}) -- skipping it", flush=True)
            worst = max(worst, 1)

    if not sources:
        print("\nno dataset to run", file=sys.stderr)
        return 2
    # ONE run directory for every size: `execute` has always taken a list of sources and
    # written them into one, with the dataset as a column in every CSV.
    print(f"\n{'=' * 60}\nsweeping {len(sources)} size(s) into one run directory\n"
          f"{'=' * 60}", flush=True)
    t0 = time.time()
    worst = max(worst, execute(spec, sources, plan(spec, sources)))
    print(f"series finished in {(time.time() - t0) / 3600:.2f} h", flush=True)
    return worst


if __name__ == "__main__":
    raise SystemExit(main())
