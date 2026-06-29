"""Shared "discover samples -> run bootstrap p-sweep -> collect results" loop.

Single source of truth for the loop that was copy-pasted across the
``02_real_data_sweeps`` notebooks and ``scripts/build_sweeps.py``. Every value
here feeds ``compute_sweep_key`` (or the pool ``param_key``), so changing any
of them moves the cache keys -- keep them byte-identical to the notebook.

Canonical cache keys at the defaults below:
``sign=acb9d64195e6``, ``sigma2=6f5687277cea``, ``kmeans=234a0bd360f7``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from .eta_pool_cache import (
    param_key, list_completed_samples, load_pool_entry, pool_root,
)
from .sweep_cache import compute_or_load_sweep, compute_sweep_key, load_sweep_result
from ..core.utils import compute_normalized_laplacian, compute_fiedler_from_laplacian

# Canonical cache location. The eta_pool / sweep caches ignore the cache_root
# arg (they use src.cache_io.CACHE_ROOT), but we pass this for compatibility.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_CACHE_ROOT = _PROJECT_ROOT / "src" / "cache"

# --- method specs: name -> (laplacian, partition_method, min_split) ----------
METHOD_SPECS: Dict[str, Tuple[str, str, int]] = {
    "sign":   ("unnormalized", "sign",   1),
    "sigma2": ("unnormalized", "sigma2", 5),
    "kmeans": ("normalized",   "kmeans", 5),
}

# --- eta bins + per-eta sample caps (same caps the notebook reads with) ------
ETA_TARGETS: List[int] = [1, 5, 10, 15]
ETA_SAMPLE_CAP: Dict[int, int] = {1: 10, 5: 10, 10: 3, 15: 3}

# --- default sweep params (MUST match the notebook for cache keys) -----------
P_VALUES: List[float] = np.geomspace(1e-3, 1.0, 25).tolist()
BOOTSTRAP_REPS = 10
SEED = 0
NUM_GAPS = 10
EARLY_STOP_CONSECUTIVE_100 = 2

# --- default pool-key params (MUST match build_eta_pool defaults) ------------
MU, POP_SIZE, SEQ_LEN, TREE_MODEL, SEQ_MODEL = 0.1, 1.0, 10_000, "kingman", "JC69"

_METRIC_KEYS = {
    "nmi": "partition_nmi_M",
    "ari": "partition_ari_M",
    "agreement": "partition_agreement_M",
}


def _default_pool_params() -> Dict[str, object]:
    return dict(seq_len=SEQ_LEN, mu=MU, tree_model=TREE_MODEL,
               pop_size=POP_SIZE, seq_model=SEQ_MODEL)


def _default_sweep_params() -> Dict[str, object]:
    return dict(p_values=P_VALUES, bootstrap_reps=BOOTSTRAP_REPS, seed=SEED,
               num_gaps=NUM_GAPS,
               early_stop_consecutive_100=EARLY_STOP_CONSECUTIVE_100)


def _norm_ref(M: np.ndarray) -> np.ndarray:
    return compute_fiedler_from_laplacian(compute_normalized_laplacian(M))


def discover_ns_and_samples(
    cache_root: Path,
    pool_params: Optional[Dict[str, object]] = None,
    eta_targets: List[int] = ETA_TARGETS,
    max_per_eta: Optional[int] = None,
) -> Tuple[List[int], Dict[Tuple[int, int], List[int]]]:
    """Scan the pool for completed samples.

    Returns ``(ns_sorted, samples_by_n_eta)`` where
    ``samples_by_n_eta[(n, eta)]`` is the list of completed sample idxs,
    capped by ``max_per_eta`` (override) or :data:`ETA_SAMPLE_CAP`.
    """
    pp = pool_params if pool_params is not None else _default_pool_params()
    key_suffix = param_key(n=0, **pp).split("_", 1)[1]

    root = pool_root(cache_root)
    ns: List[int] = []
    samples_by_n_eta: Dict[Tuple[int, int], List[int]] = {}
    for d in (sorted(root.iterdir()) if root.exists() else []):
        if not d.is_dir() or not d.name.endswith(key_suffix):
            continue
        head = d.name.split("_", 1)[0]  # e.g. 'n0500'
        try:
            n = int(head[1:])
        except ValueError:
            continue
        have_any = False
        for eta in eta_targets:
            cap = max_per_eta if max_per_eta is not None else ETA_SAMPLE_CAP.get(eta, 10)
            idxs = list_completed_samples(cache_root, d.name, eta)[:cap]
            samples_by_n_eta[(n, eta)] = idxs
            have_any = have_any or bool(idxs)
        if have_any:
            ns.append(n)
    return sorted(set(ns)), samples_by_n_eta


def run_and_collect_sweeps(
    ns: List[int],
    eta_targets: List[int] = ETA_TARGETS,
    methods: List[str] = list(METHOD_SPECS),
    pool_params: Optional[Dict[str, object]] = None,
    sweep_params: Optional[Dict[str, object]] = None,
    metric: str = "nmi",
    max_per_eta: Optional[int] = None,
    use_cache: bool = True,
    cached_only_methods: Optional[set] = None,
    verbose: bool = True,
    on_sweep: Optional[Callable[[str, int, int, int, bool], None]] = None,
    on_group: Optional[Callable[[int, int, int], None]] = None,
) -> Dict[str, Dict[Tuple[int, int], np.ndarray]]:
    """Run (or load) the bootstrap p-sweep for every (method, n, eta, idx).

    Returns ``results[method][(n, eta)]`` = ndarray of shape
    ``(num_samples x len(p_values))`` -- one row per pool sample idx, stacking
    each sample's per-p metric curve (matches the notebook exactly).

    ``cached_only_methods`` (e.g. ``{"sigma2"}``) are plotted only where their
    sweep is already on disk and NEVER recomputed -- so adding new ``n`` cannot
    trigger an expensive fresh run for those operators.

    ``on_sweep(name, n, eta, idx, was_cached)`` is called per computed/loaded
    sweep (``was_cached`` is ``None`` for a missing pool entry). ``on_group(n,
    eta, n_idxs)`` fires once per (n, eta). Both used by build_sweeps.py.
    """
    if metric not in _METRIC_KEYS:
        raise ValueError(f"metric must be one of {list(_METRIC_KEYS)}, got {metric!r}")
    metric_key = _METRIC_KEYS[metric]
    cached_only = cached_only_methods or set()
    pp = pool_params if pool_params is not None else _default_pool_params()
    sp = sweep_params if sweep_params is not None else _default_sweep_params()
    p_values = sp["p_values"]

    results: Dict[str, Dict[Tuple[int, int], np.ndarray]] = {m: {} for m in methods}
    for n in ns:
        key = param_key(n=n, **pp)
        for eta in eta_targets:
            cap = max_per_eta if max_per_eta is not None else ETA_SAMPLE_CAP.get(eta, 10)
            idxs = list_completed_samples(_CACHE_ROOT, key, eta)[:cap]
            if on_group is not None:
                on_group(n, eta, len(idxs))
            for name in methods:
                lap, pm, min_split = METHOD_SPECS[name]
                is_cached_only = name in cached_only
                sk = compute_sweep_key(
                    laplacian=lap, partition_method=pm, min_split=min_split, **sp,
                )[0] if is_cached_only else None
                arrs: List[List[float]] = []
                for idx in idxs:
                    if is_cached_only and load_sweep_result(
                        _CACHE_ROOT, key, eta, idx, sk) is None:
                        continue  # never compute a fresh (expensive) cached-only sweep
                    def _loader(_k=key, _e=eta, _i=idx, _lap=lap):
                        entry = load_pool_entry(_CACHE_ROOT, _k, _e, _i)
                        if entry is None:
                            return None
                        M, v_pop, *_ = entry
                        return (M, _norm_ref(M)) if _lap == "normalized" else (M, v_pop)

                    out = compute_or_load_sweep(
                        cache_root=_CACHE_ROOT, param_key=key, eta_target=eta, idx=idx,
                        M_loader=_loader, min_split=min_split,
                        partition_method=pm, laplacian=lap, use_cache=use_cache, **sp,
                    )
                    if out is None:
                        if on_sweep is not None:
                            on_sweep(name, n, eta, idx, None)
                        continue
                    res, was_cached = out
                    if on_sweep is not None:
                        on_sweep(name, n, eta, idx, was_cached)
                    c = res.get(metric_key)
                    if c is not None:
                        arrs.append(c)
                results[name][(n, eta)] = (
                    np.asarray(arrs) if arrs else np.empty((0, len(p_values)))
                )
            if verbose:
                print(f"  n={n:>4} eta={eta:>2} done", flush=True)
    return results
