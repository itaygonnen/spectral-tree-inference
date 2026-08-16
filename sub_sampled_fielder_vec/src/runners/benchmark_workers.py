"""Per-tree units of work for the benchmark, and the process pool that runs them.

Everything here is module level and communicates through globals set by the pool
initializer: macOS spawns workers, so the mapped function must be importable by
name and the per-item payload must stay small (a tree id, not a matrix).

Workers only *compute* — the parent process remains the sole writer of every
cache file, so a worker dying mid-task can never leave a half-written row.
"""
from __future__ import annotations

import os

# Pin BLAS to one thread per process before numpy is imported — N workers each
# running a dense eigh on an n×n matrix would otherwise oversubscribe the box.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import multiprocessing as mp  # noqa: E402
import traceback  # noqa: E402
from typing import Callable, Dict, Optional, Sequence  # noqa: E402

import numpy as np  # noqa: E402

from ..core.utils import (  # noqa: E402
    compute_fiedler_from_laplacian, compute_laplacian,
    compute_normalized_laplacian,
)
from ..utils.benchmark_cache import (  # noqa: E402
    SWEEP_METRIC_SOURCES, BenchmarkConfig, curve_key,
)
from ..utils.bpart_sweep_cache import bpart_sweep_raw  # noqa: E402
from ..utils.griffing import griffing_leading_eigvec  # noqa: E402
from ..utils.metrics import compute_reference_partition_and_quality  # noqa: E402
from ..utils.partition_validity import check_partition_valid_in_tree  # noqa: E402
from ..utils.screening import _eta, _tree_leaf_index  # noqa: E402
from .p_sweep_inner import bootstrap_p_sweep_simple  # noqa: E402

_LOADER: Optional[Callable] = None
_CFG: Optional[BenchmarkConfig] = None


def init_worker(loader, cfg: BenchmarkConfig) -> None:
    global _LOADER, _CFG
    _LOADER, _CFG = loader, cfg


def canonical_partitions(S: np.ndarray, D: np.ndarray,
                         cfg: BenchmarkConfig) -> Dict[str, np.ndarray]:
    """The selected (operator, threshold) bipartitions of one tree.

    Fiedler-on-S and Fiedler-on-L_sym are cut at the σ₂ gap; Griffing-on-D by the
    sign of the leading eigenvector of ``B = HDH``. One failing operator must not
    cost the whole tree, so each is guarded separately. Operators the config did
    not select are not computed at all.
    """
    wanted = set(cfg.screen_ops())
    out: Dict[str, np.ndarray] = {}
    for key, vec_fn in (
        ("S", lambda: compute_fiedler_from_laplacian(compute_laplacian(S))),
        ("Lsym", lambda: compute_fiedler_from_laplacian(
            compute_normalized_laplacian(S))),
    ):
        if key not in wanted:
            continue
        try:
            part, _, _ = compute_reference_partition_and_quality(
                vec_fn(), S, num_gaps=cfg.num_gaps, min_split=cfg.min_split)
            out[key] = np.asarray(part).astype(bool)
        except Exception:  # noqa: BLE001
            pass
    if "D" in wanted:
        try:
            out["D"] = griffing_leading_eigvec(D) >= 0
        except Exception:  # noqa: BLE001
            pass
    return out


def screen_one(tree_id: str) -> Optional[dict]:
    """One row of ``screen_table.csv``: per-operator validity + imbalance eta."""
    try:
        loaded = _LOADER(tree_id)
        if loaded is None:
            return None
        S, labels, tree, D = loaded
        row = {"tree": tree_id, "n_taxa": len(labels)}
        sidx = _tree_leaf_index(tree, labels) if tree is not None else None
        for key, part in canonical_partitions(S, D, _CFG).items():
            row[f"eta_{key}"] = float(_eta(part))
            if sidx is not None:
                try:
                    row[f"valid_{key}"] = bool(
                        check_partition_valid_in_tree(tree, part[sidx]))
                except Exception:  # noqa: BLE001
                    row[f"valid_{key}"] = False
        return row
    except Exception as exc:  # noqa: BLE001
        return {"tree": tree_id, "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def _curves_from_sweep(method: str, out: dict) -> dict:
    """Pull every per-p metric the Fiedler sweep returned, under canonical keys."""
    return {curve_key(method, metric): np.asarray(out[src], float)
            for metric, src in SWEEP_METRIC_SOURCES.items() if src in out}


def sweep_one(item) -> dict:
    """The requested sweeps for one tree, plus r(T).

    ``item`` = ``(tree_id, ti)`` or ``(tree_id, ti, methods)``. ``ti`` is the
    tree's index in the frozen cohort, so the per-tree seed ``1000 * ti`` is
    stable across resumes and matches the notebook. ``methods`` defaults to
    everything the config selects; passing a subset matters because a tree only
    needs the operators it was found *valid* for — the others would contribute a
    curve to no figure.
    """
    tree_id, ti, *rest = item
    want = set(rest[0]) if rest else set(_CFG.methods())
    try:
        loaded = _LOADER(tree_id)
        if loaded is None:
            return {"tree": tree_id, "error": "loader returned None"}
        S, _labels, _tree, D = loaded
        cfg, seed = _CFG, 1000 * int(ti)
        res: dict = {"tree": tree_id, "rT": float(np.max(D))}

        if "G" in want:
            raw = bpart_sweep_raw(D, cfg.p_values, reps=cfg.bootstrap_reps,
                                  seed_base=seed)
            # Mean over bootstrap reps, one value per p, for every metric the
            # sweep produced — not just NMI. results.json reports all of them.
            for metric, src in (("nmi", "nmi"), ("ari", "ari"),
                                ("agreement", "agreement"), ("dot", "dot")):
                if src in raw["per_p"][0]:
                    res[curve_key("G", metric)] = np.array(
                        [float(np.mean(pp[src])) for pp in raw["per_p"]], float)

        if "L" in want:
            out_L = bootstrap_p_sweep_simple(
                S, compute_fiedler_from_laplacian(compute_laplacian(S)),
                cfg.p_values, bootstrap_reps=cfg.bootstrap_reps, seed=seed,
                num_gaps=cfg.num_gaps, min_split=cfg.min_split,
                partition_method="kmeans", laplacian="unnormalized")
            res.update(_curves_from_sweep("L", out_L))

        if "Lsym" in want:
            out_LS = bootstrap_p_sweep_simple(
                S, compute_fiedler_from_laplacian(compute_normalized_laplacian(S)),
                cfg.p_values, bootstrap_reps=cfg.bootstrap_reps, seed=seed,
                num_gaps=cfg.num_gaps, min_split=cfg.min_split,
                partition_method="kmeans", laplacian="normalized")
            res.update(_curves_from_sweep("Lsym", out_LS))

        return res
    except Exception as exc:  # noqa: BLE001
        return {"tree": tree_id, "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def run_pool(fn, items: Sequence, n_workers: int, loader, cfg: BenchmarkConfig,
             on_result: Callable[[dict], None], label: str) -> None:
    """Map ``fn`` over ``items``, handing each result to ``on_result`` in the parent."""
    total = len(items)
    if not total:
        return
    every = max(1, min(10, total // 10 or 1))

    def _report(res, k: int) -> None:
        if res is not None:
            on_result(res)
        if k % every == 0 or k == total:
            print(f"  {label} {k}/{total}", flush=True)

    if n_workers <= 1:
        init_worker(loader, cfg)
        for k, item in enumerate(items, 1):
            _report(fn(item), k)
        return

    ctx = mp.get_context("spawn")
    with ctx.Pool(n_workers, initializer=init_worker,
                  initargs=(loader, cfg)) as pool:
        for k, res in enumerate(pool.imap_unordered(fn, items), 1):
            _report(res, k)
