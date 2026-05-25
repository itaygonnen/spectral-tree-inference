"""Parallel multi-process pool builder for the eta-binned matrix pool.

Mirrors ``build_eta_pool.py`` but runs ``attempt_one`` across N workers via
``multiprocessing.Pool``. Workers stage qualifying (M, v_pop) under
``<param_dir>/_staging/pid<P>_seed<S>/``; the main process is the sole
writer of canonical ``sample_NNNN/`` dirs so there is no idx race.

Saves a single shot of work for hard targets like (n=4000, eta=15) where
the single-process variant tops out at ~5min/attempt × ~1% hit rate.

Limits BLAS threads per worker to 1 (set BEFORE importing numpy) so 8
workers on 10 cores don't oversubscribe.

Example::

    cd sub_sampled_fielder_vec
    python scripts/build_eta_pool_parallel.py --n 4000 --samples-per-bin 10 \
        --num-workers 8 --max-attempts 4000 --seed-start 10000
"""
from __future__ import annotations

import os
# Limit BLAS threads per worker BEFORE numpy is imported anywhere.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import multiprocessing as mp

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from src.runners.eta_pool import attempt_one  # noqa: E402
from src.utils.eta_pool_cache import (  # noqa: E402
    all_bins_full, bin_counts_from_disk, bin_name, closest_target,
    init_manifest, load_manifest, next_sample_idx, param_key,
    pool_root, sample_dir, save_manifest,
)
from src.utils.partition_validity import check_partition_valid_in_tree  # noqa: E402
from analysis.theoretical_interpretation.utils.tree_features import (  # noqa: E402
    estimate_features_from_M,
)
import numpy as np  # noqa: E402


DEFAULT_CACHE_ROOT = Path(__file__).resolve().parent.parent / "src" / "cache"


def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


# Globals set in each worker via initializer; keeps imap payload small.
_PARAMS: Dict[str, Any] = {}
_STAGING_ROOT: Optional[Path] = None
_ETA_TARGETS: list = []
_ETA_TOL: float = 1.0


def _worker_init(params: Dict[str, Any], staging_root: str,
                 eta_targets: list, eta_tol: float) -> None:
    global _PARAMS, _STAGING_ROOT, _ETA_TARGETS, _ETA_TOL
    _PARAMS = dict(params)
    _STAGING_ROOT = Path(staging_root)
    _ETA_TARGETS = list(eta_targets)
    _ETA_TOL = float(eta_tol)
    _STAGING_ROOT.mkdir(parents=True, exist_ok=True)


def _worker(seed: int) -> Dict[str, Any]:
    """Generate one (M, v_pop, tree); stage to disk only if eta qualifies AND
    the Fiedler-sign partition is a single-edge bipartition of the tree."""
    try:
        M, v_pop, tree = attempt_one(seed=int(seed), **_PARAMS)
    except Exception as exc:  # noqa: BLE001
        return {"seed": int(seed), "error": f"{type(exc).__name__}: {exc}"}
    feats = estimate_features_from_M(M, v_pop)
    target = closest_target(float(feats["eta"]), _ETA_TARGETS, _ETA_TOL)
    if target is None:
        return {"seed": int(seed), "eta": float(feats["eta"]), "target": None}

    partition = (v_pop > 0).astype(bool)
    if not check_partition_valid_in_tree(tree, partition):
        return {
            "seed": int(seed),
            "eta": float(feats["eta"]),
            "target": int(target),
            "invalid_bipartition": True,
        }

    stage = _STAGING_ROOT / f"pid{os.getpid()}_seed{int(seed)}"
    stage.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(stage / "M.npz", M=M)
    np.savez_compressed(stage / "fiedler_ref.npz", fiedler_ref=v_pop)
    np.savez_compressed(stage / "partition.npz", partition=partition)
    tree_newick = tree.as_string(
        schema="newick", suppress_internal_node_labels=True,
    ).strip()
    (stage / "tree.txt").write_text(tree_newick)
    return {
        "seed": int(seed),
        "eta": float(feats["eta"]),
        "target": int(target),
        "feats": {k: (float(v) if hasattr(v, "__float__") else v)
                  for k, v in feats.items()},
        "staging": str(stage),
    }


def _consume(res: Dict[str, Any], cache_root: Path, key: str,
             bin_counts: Dict[int, int], samples_per_bin: int,
             params: Dict[str, Any], manifest: Dict[str, Any]) -> Tuple[bool, str]:
    """Move staging → canonical if bin still has capacity; update manifest.

    Returns (saved, label) where label is a short reason for the log line.
    """
    if "error" in res:
        return False, f"err seed={res['seed']} {res['error']}"
    target = res.get("target")
    if target is None:
        return False, f"seed={res['seed']} eta={res['eta']:.3f} no match"
    target = int(target)
    if res.get("invalid_bipartition"):
        return False, (f"seed={res['seed']} eta={res['eta']:.3f} "
                       f"-> {bin_name(target)} invalid bipartition")
    if bin_counts.get(target, 0) >= samples_per_bin:
        shutil.rmtree(res["staging"], ignore_errors=True)
        return False, (f"seed={res['seed']} eta={res['eta']:.3f} "
                       f"-> {bin_name(target)} already full")

    idx = next_sample_idx(cache_root, key, target)
    dst = sample_dir(cache_root, key, target, idx)
    dst.mkdir(parents=True, exist_ok=True)
    src = Path(res["staging"])
    # rename within the same filesystem; falls back to copy if cross-device.
    for fname in ("M.npz", "fiedler_ref.npz", "partition.npz", "tree.txt"):
        os.replace(src / fname, dst / fname)
    metadata = {
        "seed": int(res["seed"]),
        "eta_target": int(target),
        "eta": float(res["feats"]["eta"]),
        "n1": int(res["feats"]["n1"]),
        "n2": int(res["feats"]["n2"]),
        "S_in_max": float(res["feats"]["S_in_max"]),
        "S_out_max": float(res["feats"]["S_out_max"]),
        "S_out_min": float(res["feats"]["S_out_min"]),
        "rho": float(res["feats"]["rho"]),
        "margin": float(res["feats"]["margin"]),
        **params,
    }
    with open(dst / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    (dst / ".complete").touch()
    try:
        src.rmdir()
    except OSError:
        pass

    bin_counts[target] = bin_counts.get(target, 0) + 1
    manifest["bin_counts"][str(target)] = int(bin_counts[target])
    manifest["seeds_per_bin"].setdefault(str(target), []).append(int(res["seed"]))
    save_manifest(cache_root, key, manifest)
    return True, (f"saved seed={res['seed']} eta={res['eta']:.3f} "
                  f"-> {bin_name(target)} idx={idx} "
                  f"(n1={res['feats']['n1']}, n2={res['feats']['n2']}, "
                  f"margin={res['feats']['margin']:.4f})")


def _format_counts(bc: Dict[int, int], target: int) -> str:
    return "  ".join(f"eta{int(t):02d}={int(c)}/{target}"
                     for t, c in sorted(bc.items()))


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--n", type=int, default=4000)
    p.add_argument("--mu", type=float, default=0.1)
    p.add_argument("--pop-size", type=float, default=1.0)
    p.add_argument("--seq-len", type=int, default=10_000)
    p.add_argument("--tree-model", default="kingman")
    p.add_argument("--seq-model", default="JC69")
    p.add_argument("--eta-targets", type=int, nargs="+", default=[1, 5, 10, 15])
    p.add_argument("--eta-tol", type=float, default=1.0)
    p.add_argument("--samples-per-bin", type=int, default=10)
    p.add_argument("--max-attempts", type=int, default=4000)
    p.add_argument("--seed-start", type=int, default=None,
                   help="default: manifest.last_seed + 1")
    p.add_argument("--num-workers", type=int, default=8)
    p.add_argument("--chunksize", type=int, default=4)
    p.add_argument("--report-every", type=int, default=20)
    p.add_argument("--cache-root", type=Path, default=DEFAULT_CACHE_ROOT)
    return p.parse_args()


def main() -> None:
    args = _parse_args()
    cache_root = Path(args.cache_root)
    eta_targets = sorted(int(t) for t in args.eta_targets)
    params = dict(
        n=args.n, seq_len=args.seq_len, mu=args.mu,
        pop_size=args.pop_size, tree_model=args.tree_model,
        seq_model=args.seq_model,
    )
    worker_params = dict(
        n=args.n, mu=args.mu, pop_size=args.pop_size, seq_len=args.seq_len,
        tree_model=args.tree_model, seq_model=args.seq_model,
    )
    key = param_key(**params)
    staging_root = pool_root(cache_root) / key / "_staging"

    manifest = load_manifest(cache_root, key)
    if manifest is None:
        manifest = init_manifest(key=key, params=params,
                                 eta_targets=eta_targets,
                                 eta_tol=args.eta_tol,
                                 samples_per_bin=args.samples_per_bin)
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
    bin_counts: Dict[int, int] = {int(t): int(c) for t, c in disk_counts.items()}
    for t, c in bin_counts.items():
        manifest["bin_counts"][str(t)] = c
    save_manifest(cache_root, key, manifest)

    if all(c >= args.samples_per_bin for c in bin_counts.values()):
        _log(f"all bins already full: {_format_counts(bin_counts, args.samples_per_bin)}")
        return

    seed = (int(args.seed_start) if args.seed_start is not None
            else int(manifest.get("last_seed", -1)) + 1)

    _log(f"params: n={args.n} mu={args.mu} pop_size={args.pop_size} "
         f"L={args.seq_len} model={args.tree_model}; targets={eta_targets} "
         f"tol=±{args.eta_tol}; samples/bin={args.samples_per_bin}; "
         f"seed_start={seed}; max_attempts={args.max_attempts}; "
         f"num_workers={args.num_workers}")
    _log(f"current bin counts: {_format_counts(bin_counts, args.samples_per_bin)}")

    seeds = list(range(seed, seed + args.max_attempts))

    t0 = time.time()
    attempts = 0
    saves = 0
    invalid_skips = 0
    last_seen_seed = seed - 1

    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=args.num_workers,
                  initializer=_worker_init,
                  initargs=(worker_params, str(staging_root),
                            eta_targets, args.eta_tol)) as pool:
        try:
            for res in pool.imap_unordered(_worker, seeds, chunksize=args.chunksize):
                attempts += 1
                last_seen_seed = max(last_seen_seed, int(res["seed"]))
                manifest["attempts_used"] = int(manifest.get("attempts_used", 0)) + 1
                manifest["last_seed"] = last_seen_seed

                if res.get("invalid_bipartition"):
                    invalid_skips += 1

                saved, label = _consume(res, cache_root, key, bin_counts,
                                        args.samples_per_bin, params, manifest)
                if saved:
                    saves += 1
                    elapsed = time.time() - t0
                    _log(f"{label}; counts: {_format_counts(bin_counts, args.samples_per_bin)} "
                         f"[attempts={attempts}, saves={saves}, "
                         f"invalid={invalid_skips}, elapsed={elapsed:.1f}s]")
                elif attempts % args.report_every == 0:
                    elapsed = time.time() - t0
                    _log(f"{label}; attempts={attempts} saves={saves} "
                         f"invalid={invalid_skips} elapsed={elapsed:.1f}s")

                if all(c >= args.samples_per_bin for c in bin_counts.values()):
                    _log("all bins full; terminating workers")
                    pool.terminate()
                    break
        except KeyboardInterrupt:
            _log("interrupted; terminating workers")
            pool.terminate()
        except Exception:  # noqa: BLE001
            traceback.print_exc()
            pool.terminate()
        finally:
            save_manifest(cache_root, key, manifest)

    # Best-effort cleanup of any leftover staging dirs (workers that
    # returned a qualifying result but main exited before consuming).
    if staging_root.exists():
        for child in staging_root.iterdir():
            shutil.rmtree(child, ignore_errors=True)

    elapsed = time.time() - t0
    _log(f"done. attempts={attempts} saves={saves} "
         f"invalid_bipartitions={invalid_skips} elapsed={elapsed:.1f}s")
    _log(f"final counts: {_format_counts(bin_counts, args.samples_per_bin)}")


if __name__ == "__main__":
    main()
