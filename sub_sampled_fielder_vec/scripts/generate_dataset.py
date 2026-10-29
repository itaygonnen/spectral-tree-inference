#!/usr/bin/env python3
"""Simulate a dataset with IQ-TREE + Sailfish, in the layout the runners already read.

    python scripts/generate_dataset.py --n-taxa 1000 --n-trees 10 --dry-run
    python scripts/generate_dataset.py --n-taxa 1000 --n-trees 10
    python scripts/generate_dataset.py --check --name sailfish_yule_n1000_L1000_g1x4

Writes ``data/tree_sets/<name>/{fasta,newick}`` plus a ``manifest.json``, which is all it
takes to become a dataset: ``analysis.utils.real_datasets.list_datasets`` discovers any
directory with those two subdirectories, and from there the screen, the gate, the eta
buckets and the p-sweep are the code the downloaded datasets use. Nothing in ``src/``
needed changing to add this source.

**Generation is free; reading the result is not.** One tree costs about 0.2 s at 10,000
taxa, so ``--dry-run`` reports disk and the *downstream* screening cost instead -- that is
the number that decides whether a size is worth generating. ``src.runners.experiment_run``
owns that estimate and this script borrows it rather than keeping a second copy.

``--check`` is the contract between what this writes and what ``RealLoader`` can read:
rectangular alignments, {A,C,G,T} only, k == 4 through ``charmatrix2array``, stems paired
1:1, and FASTA headers identical to the Newick leaf labels. The last is the one that bites
-- ``RealLoader`` parses the tree into the alignment's taxon namespace, so a label that
does not match exactly is a hard failure rather than a warning.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

ROOT = Path(__file__).resolve().parent.parent          # sub_sampled_fielder_vec
REPO = ROOT.parent                                     # spectral-tree-inference
for _p in (str(ROOT), str(REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.utils.real_datasets import DATASETS_ROOT          # noqa: E402
from src.models.sailfish_sim import (SailfishSpec, generate,     # noqa: E402
                                     preflight)
from src.runners.eta_buckets import ETA_TARGETS, ETA_TOL, PER_BUCKET        # noqa: E402
from src.runners.experiment_run import screen_secs              # noqa: E402
from src.utils.dataset_check import check_dataset               # noqa: E402


# ---------------------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--n-taxa", help="taxa per tree; comma-separate for several "
                   "sizes, which writes one dataset each")
    p.add_argument("--n-trees", type=int, default=10, help="trees in the dataset")
    p.add_argument("--root-seq-size", type=int, default=1000,
                   help="root sequence length; with indels off this is the alignment "
                        "length L of every tree (default: 1000)")
    p.add_argument("--gamma-alpha", type=float, default=1.0)
    p.add_argument("--gamma-cats", type=int, default=4)
    p.add_argument("--rlen", type=float, nargs=3, default=[0.001, 0.1, 0.5],
                   metavar=("MIN", "MEAN", "MAX"), help="iqtree branch lengths")
    p.add_argument("--indels", type=float, nargs=2, default=[0.0, 0.0],
                   metavar=("INS", "DEL"),
                   help="indel rates; non-zero needs --i-know-about-k5 (default: 0 0)")
    p.add_argument("--i-know-about-k5", action="store_true",
                   help="acknowledge that gaps break the JC correction (see --indels)")
    p.add_argument("--zipf-alpha", type=float, default=1.7)
    p.add_argument("--zipf-maxlen", type=int, default=50)
    p.add_argument("--tree-seed-base", type=int, default=1000)
    p.add_argument("--msa-seed-base", type=int, default=90000)
    p.add_argument("--iqtree-bin", default=None,
                   help="override the IQ-TREE binary (default: first of iqtree3, "
                        "iqtree2, iqtree on PATH)")

    p.add_argument("--name", default=None,
                   help="dataset directory name (default: derived from the parameters)")
    p.add_argument("--out-root", type=Path, default=DATASETS_ROOT,
                   help=f"where datasets live (default: {DATASETS_ROOT})")
    p.add_argument("--resume", action="store_true",
                   help="add to an existing dataset instead of refusing it")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and the downstream cost, write nothing")
    p.add_argument("--check", action="store_true",
                   help="validate an existing dataset; non-zero exit on any violation")

    g = p.add_argument_group(
        "eta buckets",
        "Rejection-sample to the paper's eta targets instead of drawing blind: keep a "
        "tree only if its eta lands in a bucket with room AND its Fiedler partition is a "
        "real tree edge. Acceptance runs the SCREEN's eta (operator_screen.screen_one + "
        "eta_buckets.bucket_of), so an accepted tree is selected from the bucket it was "
        "accepted into -- and its screen row is cached, so it is never screened twice.")
    g.add_argument("--fill-buckets", action="store_true",
                   help="generate until every bucket is at its cap (ignores --n-trees)")
    g.add_argument("--eta-targets", default=",".join(str(t) for t in ETA_TARGETS),
                   metavar="A,B,...", help=f"default: {','.join(map(str, ETA_TARGETS))}")
    g.add_argument("--per-bucket", default=",".join(str(c) for c in PER_BUCKET),
                   metavar="A,B,...", help=f"one number per eta target, or a single number for all "
                        f"(default: {chr(44).join(map(str, PER_BUCKET))})")
    g.add_argument("--eta-tol", type=float, default=ETA_TOL,
                   help=f"bucket half-width (default: {ETA_TOL:g})")
    g.add_argument("--max-attempts", type=int, default=2000,
                   help="give up after this many trees (default: 2000)")
    g.add_argument("--workers", type=int, default=1,
                   help="parallel attempts; >1 pins BLAS to one thread per worker")
    g.add_argument("--calibrate", type=int, default=0, metavar="N",
                   help="run N attempts, report the eta distribution and hit rate, "
                        "write nothing -- the acceptance rate sets the whole cost")
    return p.parse_args(argv)


def _ints(text: str) -> List[int]:
    return [int(x) for x in str(text).replace(" ", "").split(",") if x]


def spec_from(args: argparse.Namespace) -> SailfishSpec:
    ins, dele = args.indels
    return SailfishSpec(
        pooled=bool(args.fill_buckets or args.calibrate),
        n_taxa=args.n_taxa, n_trees=args.n_trees,
        root_seq_size=args.root_seq_size,
        rlen=tuple(args.rlen),
        gamma_alpha=args.gamma_alpha, gamma_categories=args.gamma_cats,
        insertion_rate=ins, deletion_rate=dele,
        zipf_alpha=args.zipf_alpha, zipf_maxlen=args.zipf_maxlen,
        tree_seed_base=args.tree_seed_base, msa_seed_base=args.msa_seed_base,
        iqtree_bin=args.iqtree_bin,
    )


# ---------------------------------------------------------------------------
# --dry-run
# ---------------------------------------------------------------------------

def _duration(seconds: float) -> str:
    """Seconds, minutes or hours -- whichever keeps a significant figure.

    ``screen_secs`` is calibrated at m=6000 and reads low at small m (see its docstring),
    so printing a small estimate as '0.0 h' would hide an already-approximate number
    behind a rounding artefact.
    """
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    return f"{seconds / 3600:.1f} h"


def describe_buckets(spec: SailfishSpec, targets: List[int], caps: List[int],
                     tol: float) -> None:
    """The caps, in one line. The run's own progress reports the rest as it goes."""
    print(f"  buckets       "
          + ", ".join(f"{c} at eta{t}" for t, c in zip(targets, caps))
          + f"  (within +-{tol:g})")


def report_calibration(rows: List[dict], targets: List[int], caps: List[int],
                       tol: float, cap_ops, elapsed: float) -> None:
    """What N attempts say about whether the full fill is worth starting."""
    from src.runners.eta_buckets import bucket_of

    arm = cap_ops[0] if cap_ops else "S"
    ok = [r for r in rows if "error" not in r]
    etas = sorted(float(r[f"eta_{arm}"]) for r in ok if f"eta_{arm}" in r)
    valid = [r for r in ok if r.get(f"valid_{arm}")]
    # Counted twice on purpose: an empty bucket because no tree had that eta needs a
    # wider --eta-tol, while an empty bucket whose window was reached but never by a
    # valid split needs more attempts. The two call for opposite actions.
    hits = {int(t): 0 for t in targets}
    in_window = {int(t): 0 for t in targets}
    for r in ok:
        t = bucket_of(r, cap_ops, targets, tol)
        if t is None:
            continue
        in_window[int(t)] += 1
        if r.get(f"valid_{arm}"):
            hits[int(t)] += 1

    n = len(rows)
    print(f"\n{n} attempts in {_duration(elapsed)} "
          f"({_duration(elapsed / max(n, 1))} each)")
    print(f"  screened ok   {len(ok)}/{n}")
    print(f"  valid split   {len(valid)}/{len(ok) or 1}"
          f"  ({100 * len(valid) / max(len(ok), 1):.0f}%)")
    if etas:
        q = [etas[min(int(p * len(etas)), len(etas) - 1)] for p in (0.1, 0.5, 0.9)]
        print(f"  eta_{arm}         min={etas[0]:.1f}  p10={q[0]:.1f}  "
              f"median={q[1]:.1f}  p90={q[2]:.1f}  max={etas[-1]:.1f}")
    print("\n  bucket   in window   valid    rate    attempts for the cap")
    for t, cap in zip(targets, caps):
        h, w = hits[int(t)], in_window[int(t)]
        rate = h / n if n else 0.0
        need = f"{cap / rate:>6.0f}" if rate else "     -"
        if rate:
            flag = ""
        elif w:
            flag = "   <- reached, never valid: needs more attempts"
        else:
            flag = "   <- window never reached: widen --eta-tol"
        print(f"  eta{int(t):02d}   {w:>9}   {h:>5}   {100 * rate:>5.1f}%   {need}{flag}")
    total_rate = sum(hits.values()) / n if n else 0.0
    print(f"  overall accepted {100 * total_rate:.1f}%")


def describe(spec: SailfishSpec, dataset: Path, n_trees: Optional[int] = None) -> None:
    # In bucket mode the caps decide the count, so the caller passes it; printing
    # ``--n-trees`` there would state a number nothing uses.
    n = spec.n_trees if n_trees is None else n_trees
    fasta_mb = spec.n_taxa * (spec.root_seq_size + 12) * n / 1e6
    screen_s = screen_secs(spec.n_taxa, spec.root_seq_size) * n
    print(f"  name          {dataset.name}")
    print(f"  into          {dataset}")
    print(f"  trees         {n} x {spec.n_taxa} taxa")
    print(f"  alignment     L={spec.root_seq_size}, JC + Gamma("
          f"{spec.gamma_categories} cats, alpha={spec.gamma_alpha:g})")
    print(f"  indels        {'ins=%g del=%g  -- GAPPED, k>4' % (spec.insertion_rate, spec.deletion_rate) if spec.has_indels else 'off (gapless, k=4)'}")
    print(f"  branch len    min={spec.rlen[0]:g} mean={spec.rlen[1]:g} max={spec.rlen[2]:g}")
    print(f"  disk          ~{fasta_mb:.0f} MB of fasta")
    print(f"  to screen it  ~{_duration(screen_s)}  (generation itself is seconds)")


# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    if args.check:
        if not args.name:
            print("--check needs --name (or --name plus --out-root)", file=sys.stderr)
            return 2
        dataset = Path(args.out_root) / args.name
        if not dataset.is_dir():
            print(f"no such dataset: {dataset}", file=sys.stderr)
            return 2
        print(f"checking {dataset}")
        gapless = not _manifest_has_indels(dataset)
        fails = check_dataset(dataset, expect_gapless=gapless)
        for f in fails:
            print(f"  FAIL {f}")
        print("OK" if not fails else f"FAILED: {len(fails)} problem(s)")
        return 1 if fails else 0

    if not args.n_taxa:
        print("--n-taxa is required (or pass --check --name ...)", file=sys.stderr)
        return 2

    try:
        sizes = _ints(args.n_taxa)
    except ValueError:
        print(f"--n-taxa is not a list of whole numbers: {args.n_taxa!r}",
              file=sys.stderr)
        return 2
    if not sizes:
        print("--n-taxa is empty", file=sys.stderr)
        return 2
    if len(sizes) > 1 and args.name:
        print("--name fixes one directory, so it cannot be combined with several "
              "sizes; drop it and each size is named after its own parameters",
              file=sys.stderr)
        return 2

    worst = 0
    for n in sizes:
        if len(sizes) > 1:
            print(f"\n=== {n} taxa " + "=" * 46)
        worst = max(worst, _one_size(args, n))
    return worst


def _one_size(args: argparse.Namespace, n_taxa: int) -> int:
    """Everything for one taxon count. Several sizes = several datasets, run in turn."""
    args = argparse.Namespace(**{**vars(args), "n_taxa": n_taxa})
    spec = spec_from(args)
    if spec.has_indels and not args.i_know_about_k5:
        print(
            "refusing: non-zero indel rates produce gapped alignments, and\n"
            "  spectraltree/similarities.py:39 takes k from the data -- a gap becomes a\n"
            "  5th state, so the JC correction is wrong and gap-vs-base counts as a\n"
            "  substitution. Measured: data/tree_sets/'1000 taxa' is 88.3% gap at k=5.\n"
            "  Pass --i-know-about-k5 to generate one anyway.", file=sys.stderr)
        return 2

    dataset = Path(args.out_root) / (args.name or spec.name)
    targets, caps = _ints(args.eta_targets), _ints(args.per_bucket)
    if len(caps) == 1:
        caps = caps * len(targets)      # one number plainly means that many everywhere
    bucket_mode = bool(args.fill_buckets or args.calibrate)
    if bucket_mode and len(targets) != len(caps):
        named = ", ".join(f"eta{t}" for t in targets)
        print(f"--per-bucket needs one number per target ({named}), or a single "
              f"number for all; got {len(caps)} for {len(targets)} targets",
              file=sys.stderr)
        return 2

    print("plan")
    describe(spec, dataset, n_trees=sum(caps) if bucket_mode else None)
    if bucket_mode:
        describe_buckets(spec, targets, caps, args.eta_tol)
        if args.n_trees != 10:
            print("  note          --n-trees is ignored here; the caps decide the count")

    if args.dry_run:
        print("\n(dry run -- nothing written)")
        return 0

    missing = preflight(spec)
    if missing:
        print("\ncannot run -- missing dependencies:", file=sys.stderr)
        for line in missing:
            print(f"  - {line}", file=sys.stderr)
        return 2

    if args.calibrate:
        return run_calibration(spec, args, targets, caps)

    if args.fill_buckets:
        return run_fill(spec, args, dataset, targets, caps)

    if dataset.exists() and any(dataset.iterdir()) and not args.resume:
        print(f"\nrefusing: {dataset} already exists and is not empty. Pass --resume to "
              f"add to it, or --name to write a separate dataset.", file=sys.stderr)
        return 2

    print()
    generate(spec, Path(args.out_root), name=args.name)

    print("\nchecking what was written")
    fails = check_dataset(dataset, expect_gapless=not spec.has_indels)
    for f in fails:
        print(f"  FAIL {f}")
    if fails:
        print(f"FAILED: {len(fails)} problem(s)")
        return 1
    print(f"OK -- run it with: python scripts/interactive_run.py  "
          f"(appears as '{dataset.name}')")
    return 0


def run_calibration(spec: SailfishSpec, args, targets, caps) -> int:
    """``--calibrate N``: measure the hit rate, write nothing, keep no files."""
    import shutil
    import tempfile
    import time

    from src.runners.experiment_run import gate_operators
    from src.runners.sailfish_pool import _attempts
    from src.models.sailfish_sim import find_iqtree

    cap_ops = gate_operators("valid_L")
    staging = Path(tempfile.mkdtemp(prefix="sailfish_calib_"))
    binary = find_iqtree(spec)
    print(f"\ncalibrating on {args.calibrate} attempts "
          f"({args.workers} worker{'s' if args.workers > 1 else ''})")
    t0 = time.time()
    rows = []
    try:
        seeds = iter(range(spec.tree_seed_base, spec.tree_seed_base + args.calibrate))
        for row in _attempts(seeds, spec, staging, ("S",), binary, args.workers,
                             spec_policy="best"):
            rows.append(row)
            shutil.rmtree(row.get("_staging", ""), ignore_errors=True)
            print(".", end="", flush=True)
    except KeyboardInterrupt:
        print("\ninterrupted -- reporting what was measured")
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    if not rows:
        print("\nno attempts completed", file=sys.stderr)
        return 1
    report_calibration(rows, targets, caps, args.eta_tol, cap_ops, time.time() - t0)
    print("\n(calibration -- nothing written)")
    return 0


def run_fill(spec: SailfishSpec, args, dataset: Path, targets, caps) -> int:
    """``--fill-buckets``: rejection-sample until every bucket is at its cap."""
    from src.runners.sailfish_pool import fill_buckets

    print()
    stats = fill_buckets(
        spec, Path(args.out_root), name=args.name,
        targets=targets, per_bucket=caps, tol=args.eta_tol,
        max_attempts=args.max_attempts, workers=args.workers,
        seed_start=None)

    got = [stats["counts"][int(t)] for t in targets]
    print(f"\n{stats['accepted']} accepted of {stats['attempts']} attempts "
          f"({100 * stats['acceptance']:.1f}%) in {_duration(stats['elapsed_sec'])}")
    print(f"  buckets      {got} of {caps}")
    print(f"  rejected     {stats['no_bucket']} no bucket, {stats['invalid']} invalid "
          f"split, {stats['full']} bucket full, {stats['errors']} errors")

    # the manifest is written by fill_buckets on every accept, so a killed run still
    # knows which seeds it tried -- nothing to write here

    print("\nchecking what was written")
    fails = check_dataset(dataset, expect_gapless=not spec.has_indels)
    for f in fails:
        print(f"  FAIL {f}")
    if fails:
        return 1
    if not stats["complete"]:
        why = {"signal": "stopped by a signal (time limit?)",
               "exhausted": f"ran out of --max-attempts ({args.max_attempts})"}.get(
                   stats.get("stopped"), "stopped early")
        short = [f"eta{t:02d} {stats['counts'][int(t)]}/{c}"
                 for t, c in zip(targets, caps) if stats["counts"][int(t)] < c]
        print(f"\nINCOMPLETE after {stats['attempts']} attempts -- {why}. "
              f"Short: {', '.join(short)}.\n"
              f"  Re-run the same command to continue from seed "
              f"{stats['last_seed'] + 1} (nothing is re-drawn).\n"
              f"  A bucket that never fills wants a wider --eta-tol, not more attempts "
              f"-- check with --calibrate.")
        return 1
    print(f"OK -- every bucket full. Sweep it without screening:\n"
          f"  python scripts/run_sweep.py --dataset {dataset.name} --stage sweep")
    return 0


def _manifest_has_indels(dataset: Path) -> bool:
    """Whether a written dataset was generated with indels on, per its manifest."""
    import json
    path = dataset / "manifest.json"
    if not path.exists():
        return False
    try:
        spec = json.loads(path.read_text()).get("spec", {})
    except (OSError, ValueError):
        return False
    return bool(spec.get("insertion_rate") or spec.get("deletion_rate"))


if __name__ == "__main__":
    raise SystemExit(main())
