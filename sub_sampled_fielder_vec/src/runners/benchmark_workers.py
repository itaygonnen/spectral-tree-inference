"""Per-tree units of work for the benchmark, and the process pool that runs them.

The *compute* is no longer here: screening a tree and sweeping a tree are
:mod:`src.runners.operator_screen` and :mod:`src.runners.operator_sweep`, which the
real-data pipeline also uses. This module is the adapter — it calls them with the
settings this benchmark is fixed to, and reshapes the result into the row and curve
schema ``benchmark_cache`` writes.

Those settings are not defaults and must not drift to them:

    screen   the sigma2 gap search  (``rule_policy="sigma2"``)
    sweep    k-means on the Fiedler arms, and ``per_rep`` aggregation with the dense
             solver on the distance arm — ``bpart_sweep_cache`` records that this
             accounting is like-for-like with a published figure.

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
from typing import Callable, Optional, Sequence  # noqa: E402

import numpy as np  # noqa: E402

from ..utils.benchmark_cache import BenchmarkConfig, curve_key  # noqa: E402
from ..utils.griffing import DEFAULT_SOLVER  # noqa: E402
from .operator_screen import screen_one as _screen_one  # noqa: E402
from .operator_sweep import sweep_one_tree as _sweep_one_tree  # noqa: E402
from .operators import ARM_OF  # noqa: E402

# This benchmark names the distance operator "D" and the similarity Fiedler arm's
# curve "L"; the shared table calls them "B" and "L". One mapping, here.
_OP_OF_SCREEN = {"S": "S", "Lsym": "Lsym", "D": "B"}
_SCREEN_OF_OP = {v: k for k, v in _OP_OF_SCREEN.items()}
_METHOD_OF_ARM = {"L": "L", "Lsym": "Lsym", "B": "G"}

_LOADER: Optional[Callable] = None
_CFG: Optional[BenchmarkConfig] = None


def init_worker(loader, cfg: BenchmarkConfig) -> None:
    global _LOADER, _CFG
    _LOADER, _CFG = loader, cfg


def screen_one(tree_id: str) -> Optional[dict]:
    """One row of ``screen_table.csv``: per-operator validity + imbalance eta.

    The shared screen keys its row by operator ("S", "Lsym", "B") and records both
    candidate rules; this benchmark's table wants one eta and one flag per operator
    under its own key ("S", "Lsym", "D"), from the sigma2 gap search.
    """
    try:
        ops = [_OP_OF_SCREEN[k] for k in _CFG.screen_ops()]
        rec = _screen_one(tree_id, _LOADER, operators=ops,
                          min_split=_CFG.min_split, rule_policy="sigma2",
                          num_gaps=_CFG.num_gaps)
        if rec is None:
            return None
        if "error" in rec:
            return rec
        row = {"tree": tree_id, "n_taxa": int(rec.get("m", 0))}
        for op in ops:
            key = _SCREEN_OF_OP[op]
            if f"eta_{op}" in rec:
                row[f"eta_{key}"] = float(rec[f"eta_{op}"])
            if f"valid_{op}" in rec:
                row[f"valid_{key}"] = bool(rec[f"valid_{op}"])
        return row
    except Exception as exc:  # noqa: BLE001
        return {"tree": tree_id, "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=3)}


def sweep_one(item) -> dict:
    """The requested sweeps for one tree, plus r(T).

    ``item`` = ``(tree_id, ti)`` or ``(tree_id, ti, methods)``. ``ti`` is the
    tree's index in the frozen cohort, so the per-tree seed ``1000 * ti`` is
    stable across resumes and matches the notebook. ``methods`` defaults to
    everything the config selects; passing a subset matters because a tree only
    needs the operators it was found *valid* for — the others would contribute a
    curve to no figure.

    The sweep itself is the shared one; this maps its arm names onto the curve keys
    ``benchmark_cache`` stores ("G" for the distance arm) and re-reads r(T), which is
    the scale plot's x-axis and not something a recovery sweep produces.
    """
    tree_id, ti, *rest = item
    want = set(rest[0]) if rest else set(_CFG.methods())
    ops = [k for k, arm in ARM_OF.items() if _METHOD_OF_ARM[arm] in want]
    try:
        loaded = _LOADER(tree_id)
        if loaded is None:
            return {"tree": tree_id, "error": "loader returned None"}
        _S, _labels, _tree, D = loaded
        cfg = _CFG
        curves = _sweep_one_tree(
            tree_id, _LOADER, seed=1000 * int(ti), p_values=cfg.p_values,
            reps=cfg.bootstrap_reps, num_gaps=cfg.num_gaps,
            min_split=cfg.min_split, operators=ops,
            # fixed for this benchmark; see the module docstring
            rule_policy="kmeans", eigsolver=DEFAULT_SOLVER, aggregation="per_rep")

        res: dict = {"tree": tree_id, "rT": float(np.max(D))}
        for op in ops:
            method = _METHOD_OF_ARM[ARM_OF[op]]
            for metric in ("nmi", "ari", "agreement", "dot"):
                col = f"{metric}_{ARM_OF[op]}_vs_fullmatrix"
                if col in curves:
                    res[curve_key(method, metric)] = np.asarray(curves[col], float)
            sign_col = f"signagreement_{ARM_OF[op]}_vs_fullmatrix"
            if sign_col in curves:
                res[curve_key(method, "sign")] = np.asarray(curves[sign_col], float)
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
