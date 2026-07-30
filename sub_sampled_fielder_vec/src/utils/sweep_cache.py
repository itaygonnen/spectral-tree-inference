"""Disk cache for ``bootstrap_p_sweep_simple`` outputs.

Thin shim over :mod:`src.cache_io`. Backing scope is ``bootstrap_sweep``;
``cache_root`` is accepted for backward compat but ignored.

Each sweep entry is addressed by the flat key
``<param_key>__<eta_name>__<sample_name>__<sweep_key>`` to keep all
sweeps in one scope-level directory rather than nested inside pool
samples (where they lived in the legacy layout).
"""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from ..cache_io import bootstrap_sweep as _scope, SENTINEL
from .eta_pool_cache import bin_name


SCHEMA_VERSION = "v3"


def _flat_key(param_key: str, eta_target: int, idx: int, sweep_key: str) -> str:
    return f"{param_key}__{bin_name(eta_target)}__sample_{int(idx):04d}__{sweep_key}"


def sweep_cache_dir(
    cache_root: Path, param_key: str, eta_target: int, idx: int
) -> Path:
    """Legacy helper: previously the per-sample ``sweeps/`` parent.

    Now redundant — sweep entries are flat under :data:`bootstrap_sweep`.
    Returned path stays compatible (the parent of any sweep this trio
    would have produced), but it is no longer used by ``load``/``save``.
    """
    return _scope.root / f"{param_key}__{bin_name(eta_target)}__sample_{int(idx):04d}"


def _build_config(
    p_values: List[float],
    bootstrap_reps: int,
    seed: int,
    num_gaps: int,
    min_split: int,
    early_stop_consecutive_100: int,
    partition_method: str,
    matrix_kind: str = "similarity",
    distance_alpha: float = 1.0,
    laplacian: str = "unnormalized",
) -> Dict[str, Any]:
    config = {
        "p_values": [round(float(p), 9) for p in p_values],
        "bootstrap_reps": int(bootstrap_reps),
        "seed": int(seed),
        "num_gaps": int(num_gaps),
        "min_split": int(min_split),
        "early_stop_consecutive_100": int(early_stop_consecutive_100),
        "partition_method": str(partition_method),
        "matrix_kind": str(matrix_kind),
        "distance_alpha": round(float(distance_alpha), 9),
        "schema_version": SCHEMA_VERSION,
    }
    # Only tag non-default Laplacians so existing unnormalized keys stay byte-identical.
    if laplacian != "unnormalized":
        config["laplacian"] = str(laplacian)
    return config


def compute_sweep_key(
    p_values: List[float],
    bootstrap_reps: int,
    seed: int,
    num_gaps: int,
    min_split: int,
    early_stop_consecutive_100: int,
    partition_method: str,
    matrix_kind: str = "similarity",
    distance_alpha: float = 1.0,
    laplacian: str = "unnormalized",
) -> Tuple[str, Dict[str, Any]]:
    config = _build_config(
        p_values, bootstrap_reps, seed, num_gaps, min_split,
        early_stop_consecutive_100, partition_method,
        matrix_kind=matrix_kind, distance_alpha=distance_alpha,
        laplacian=laplacian,
    )
    blob = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha1(blob).hexdigest()[:12], config


def load_sweep_result(
    cache_root: Path, param_key: str, eta_target: int, idx: int, sweep_key: str
) -> Optional[Dict[str, Any]]:
    flat = _flat_key(param_key, eta_target, idx, sweep_key)
    d = _scope.path(flat)
    if not (d / SENTINEL).exists():
        return None
    with open(d / "result.json", "r") as f:
        payload = json.load(f)
    result = payload["result"]
    split = result.get("partition_split_ref")
    if isinstance(split, list) and len(split) == 2:
        result["partition_split_ref"] = (int(split[0]), int(split[1]))
    return result


def save_sweep_result(
    cache_root: Path,
    param_key: str,
    eta_target: int,
    idx: int,
    sweep_key: str,
    config: Dict[str, Any],
    result: Dict[str, Any],
) -> Path:
    flat = _flat_key(param_key, eta_target, idx, sweep_key)
    d = _scope.path(flat)
    d.mkdir(parents=True, exist_ok=True)
    sentinel = d / SENTINEL
    if sentinel.exists():
        sentinel.unlink()
    payload = {
        "config": config,
        "result": result,
        "saved_at": datetime.now(timezone.utc).isoformat(),
    }
    tmp = d / "result.json.tmp"
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    os.replace(tmp, d / "result.json")
    sentinel.touch()
    return d


def compute_or_load_sweep(
    cache_root: Path,
    param_key: str,
    eta_target: int,
    idx: int,
    M_loader: Callable[[], Optional[Tuple[Any, Any]]],
    *,
    p_values: List[float],
    bootstrap_reps: int = 50,
    seed: int = 0,
    num_gaps: int = 10,
    min_split: int = 1,
    early_stop_consecutive_100: int = 0,
    partition_method: str = "sigma2",
    matrix_kind: str = "similarity",
    distance_alpha: float = 1.0,
    laplacian: str = "unnormalized",
    use_cache: bool = True,
) -> Optional[Tuple[Dict[str, Any], bool]]:
    """Return ``(result, was_cached)`` for one sweep, or ``None`` if M_loader returns None."""
    from ..runners.p_sweep_inner import bootstrap_p_sweep_simple

    sweep_key, config = compute_sweep_key(
        p_values=p_values, bootstrap_reps=bootstrap_reps, seed=seed,
        num_gaps=num_gaps, min_split=min_split,
        early_stop_consecutive_100=early_stop_consecutive_100,
        partition_method=partition_method,
        matrix_kind=matrix_kind, distance_alpha=distance_alpha,
        laplacian=laplacian,
    )

    if use_cache:
        cached = load_sweep_result(cache_root, param_key, eta_target, idx, sweep_key)
        if cached is not None:
            return cached, True

    loaded = M_loader()
    if loaded is None:
        return None
    M, v_pop = loaded

    result = bootstrap_p_sweep_simple(
        M=M, fiedler_ref=v_pop, p_values=p_values,
        bootstrap_reps=bootstrap_reps, seed=seed, num_gaps=num_gaps,
        min_split=min_split,
        early_stop_consecutive_100=early_stop_consecutive_100,
        partition_method=partition_method,
        laplacian=laplacian,
    )
    save_sweep_result(
        cache_root, param_key, eta_target, idx, sweep_key, config, result,
    )
    return result, False
