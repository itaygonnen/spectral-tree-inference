"""Build a pool of (M, v_pop, partition) triples binned by realized eta.

Generate Kingman trees, simulate JC69 sequences, build the similarity matrix
M and its reference Fiedler vector v_pop, estimate eta = max(n1, n2) /
min(n1, n2), and persist only those samples whose eta lands within
``--eta-tol`` of one of the ``--eta-targets`` AND whose sign(v_pop) partition
is a single-edge bipartition of the underlying tree. Each qualifying sample
lands in ``cache_root/pool_sample/<param_key>/eta<TT>/sample_NNNN/`` and
contains ``M.npz``, ``fiedler_ref.npz``, ``partition.npz`` (bool, len n),
``metadata.json``, and a ``.complete`` sentinel.

The pool is resumable: re-running the script with the same parameters reads
``manifest.json`` and continues seeding from ``last_seed + 1``, skipping
bins that are already full.

First-pass defaults match the user's request: n=500, mu=0.1, pop_size=1.0,
kingman, eta targets {1, 5, 10, 15} +/- 1.

Example::

    cd sub_sampled_fielder_vec
    python scripts/build_eta_pool.py
    python scripts/build_eta_pool.py --samples-per-bin 50 --max-attempts 5000
    python scripts/build_eta_pool.py --n 1000 --max-attempts 4000
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from src.runners.eta_pool import attempt_one
from src.utils.eta_pool_cache import (
    all_bins_full,
    bin_counts_from_disk,
    closest_target,
    init_manifest,
    load_manifest,
    next_sample_idx,
    param_key,
    remaining_capacity,
    save_manifest,
    save_pool_entry,
)
from src.utils.partition_validity import check_partition_valid_in_tree
from analysis.utils.tree_features import (
    estimate_features_from_M,
)


DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent.parent / "src" / "cache"


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--n", type=int, default=500, help="number of taxa")
    p.add_argument("--mu", type=float, default=0.1, help="JC69 mutation rate")
    p.add_argument("--pop-size", type=float, default=1.0, help="Kingman pop_size")
    p.add_argument("--seq-len", type=int, default=10_000, help="simulated sequence length")
    p.add_argument("--tree-model", default="kingman",
                   choices=["kingman", "kingman_mean", "lopsided", "birth_death", "balanced_binary"])
    p.add_argument("--seq-model", default="JC69")
    p.add_argument(
        "--eta-targets",
        type=int,
        nargs="+",
        default=[1, 5, 10, 15],
        help="integer target etas",
    )
    p.add_argument("--eta-tol", type=float, default=1.0, help="bin half-width on |eta - target|")
    p.add_argument(
        "--samples-per-bin",
        type=int,
        default=20,
        help="how many qualifying (M, v_pop) pairs per target eta",
    )
    p.add_argument(
        "--max-attempts",
        type=int,
        default=2000,
        help="cap on tree-generation attempts in this invocation",
    )
    p.add_argument(
        "--seed-start",
        type=int,
        default=None,
        help="override resume seed (default: manifest.last_seed + 1)",
    )
    p.add_argument(
        "--cache-root",
        type=Path,
        default=DEFAULT_CACHE_ROOT,
        help=f"cache root (default: {DEFAULT_CACHE_ROOT})",
    )
    p.add_argument(
        "--report-every",
        type=int,
        default=10,
        help="print progress every N attempts",
    )
    return p.parse_args()


def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


def _format_counts(counts_by_target: dict, target_total: int) -> str:
    parts = [f"eta{int(t):02d}={int(c)}/{target_total}" for t, c in sorted(counts_by_target.items(), key=lambda kv: int(kv[0]))]
    return "  ".join(parts)


def main() -> None:
    args = _parse_args()
    cache_root = Path(args.cache_root)
    eta_targets: List[int] = sorted(int(t) for t in args.eta_targets)

    key = param_key(
        n=args.n,
        seq_len=args.seq_len,
        mu=args.mu,
        tree_model=args.tree_model,
        pop_size=args.pop_size,
        seq_model=args.seq_model,
    )
    params = dict(
        n=args.n,
        seq_len=args.seq_len,
        mu=args.mu,
        pop_size=args.pop_size,
        tree_model=args.tree_model,
        seq_model=args.seq_model,
    )

    manifest = load_manifest(cache_root, key)
    if manifest is None:
        manifest = init_manifest(
            key=key,
            params=params,
            eta_targets=eta_targets,
            eta_tol=args.eta_tol,
            samples_per_bin=args.samples_per_bin,
        )
        _log(f"new pool: {key}")
    else:
        _log(f"resuming pool: {key}")
        manifest["samples_per_bin"] = int(args.samples_per_bin)
        manifest["eta_tol"] = float(args.eta_tol)
        for t in eta_targets:
            manifest["bin_counts"].setdefault(str(t), 0)
            manifest["seeds_per_bin"].setdefault(str(t), [])
        manifest["eta_targets"] = sorted({*eta_targets, *manifest.get("eta_targets", [])})

    disk_counts = bin_counts_from_disk(cache_root, key, manifest["eta_targets"])
    for t, c in disk_counts.items():
        manifest["bin_counts"][str(t)] = int(c)

    if all_bins_full(manifest):
        _log(f"all bins already full: {_format_counts(manifest['bin_counts'], manifest['samples_per_bin'])}")
        save_manifest(cache_root, key, manifest)
        return

    if args.seed_start is not None:
        seed = int(args.seed_start)
    else:
        seed = int(manifest.get("last_seed", -1)) + 1
    _log(
        f"params: n={args.n} mu={args.mu} pop_size={args.pop_size} L={args.seq_len} "
        f"model={args.tree_model}; targets={eta_targets} tol=±{args.eta_tol}; "
        f"samples/bin={args.samples_per_bin}; starting seed={seed}; "
        f"max_attempts={args.max_attempts}"
    )
    _log(f"current bin counts: {_format_counts(manifest['bin_counts'], manifest['samples_per_bin'])}")

    t0 = time.time()
    attempts_this_run = 0
    saves_this_run = 0
    invalid_skips = 0
    try:
        for _ in range(args.max_attempts):
            if all_bins_full(manifest):
                break
            attempts_this_run += 1
            manifest["attempts_used"] = int(manifest.get("attempts_used", 0)) + 1
            manifest["last_seed"] = seed

            try:
                M, v_pop, tree = attempt_one(
                    seed=seed,
                    n=args.n,
                    mu=args.mu,
                    pop_size=args.pop_size,
                    seq_len=args.seq_len,
                    tree_model=args.tree_model,
                    seq_model=args.seq_model,
                )
            except Exception as exc:
                _log(f"attempt seed={seed} failed: {exc!r}; skipping")
                seed += 1
                continue

            feats = estimate_features_from_M(M, v_pop)
            target = closest_target(feats["eta"], eta_targets, args.eta_tol)

            if target is None or remaining_capacity(manifest, target) == 0:
                if attempts_this_run % args.report_every == 0:
                    elapsed = time.time() - t0
                    _log(
                        f"seed={seed} eta={feats['eta']:.3f} -> "
                        f"{'no match' if target is None else f'eta{target:02d} full'};"
                        f" attempts={attempts_this_run} saves={saves_this_run} "
                        f"invalid={invalid_skips} elapsed={elapsed:.1f}s"
                    )
                seed += 1
                continue

            partition = (v_pop > 0).astype(bool)
            if not check_partition_valid_in_tree(tree, partition):
                invalid_skips += 1
                if attempts_this_run % args.report_every == 0:
                    elapsed = time.time() - t0
                    _log(
                        f"seed={seed} eta={feats['eta']:.3f} -> "
                        f"eta{target:02d} invalid bipartition;"
                        f" attempts={attempts_this_run} saves={saves_this_run} "
                        f"invalid={invalid_skips} elapsed={elapsed:.1f}s"
                    )
                seed += 1
                continue

            idx = next_sample_idx(cache_root, key, target)
            metadata = {
                "seed": int(seed),
                "eta_target": int(target),
                "eta": float(feats["eta"]),
                "n1": int(feats["n1"]),
                "n2": int(feats["n2"]),
                "S_in_max": float(feats["S_in_max"]),
                "S_out_max": float(feats["S_out_max"]),
                "S_out_min": float(feats["S_out_min"]),
                "rho": float(feats["rho"]),
                "margin": float(feats["margin"]),
                **params,
            }
            tree_newick = tree.as_string(
                schema="newick", suppress_internal_node_labels=True,
            ).strip()
            save_pool_entry(
                cache_root, key, target, idx,
                M, v_pop, metadata,
                partition=partition, tree_newick=tree_newick,
            )
            saves_this_run += 1

            manifest["bin_counts"][str(int(target))] = (
                int(manifest["bin_counts"].get(str(int(target)), 0)) + 1
            )
            manifest["seeds_per_bin"].setdefault(str(int(target)), []).append(int(seed))
            save_manifest(cache_root, key, manifest)

            elapsed = time.time() - t0
            _log(
                f"saved seed={seed} eta={feats['eta']:.3f} -> eta{target:02d} "
                f"(idx={idx}, n1={feats['n1']}, n2={feats['n2']}, margin={feats['margin']:.4f}); "
                f"counts: {_format_counts(manifest['bin_counts'], manifest['samples_per_bin'])} "
                f"[attempts={attempts_this_run}, elapsed={elapsed:.1f}s]"
            )
            seed += 1
    except KeyboardInterrupt:
        _log("interrupted by user; persisting manifest")
    finally:
        save_manifest(cache_root, key, manifest)

    elapsed = time.time() - t0
    _log(
        f"done. attempts_this_run={attempts_this_run} saves_this_run={saves_this_run} "
        f"invalid_bipartitions={invalid_skips} "
        f"total_attempts={manifest['attempts_used']} elapsed={elapsed:.1f}s"
    )
    _log(f"final counts: {_format_counts(manifest['bin_counts'], manifest['samples_per_bin'])}")
    if all_bins_full(manifest):
        _log("all bins are full.")
    else:
        missing = {
            t: max(0, manifest["samples_per_bin"] - int(manifest["bin_counts"].get(str(t), 0)))
            for t in manifest["eta_targets"]
        }
        unfilled = {t: m for t, m in missing.items() if m > 0}
        _log(f"unfilled bins: {unfilled}; rerun to continue.")


if __name__ == "__main__":
    main()
