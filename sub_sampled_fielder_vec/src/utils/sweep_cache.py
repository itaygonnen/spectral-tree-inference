"""Disk cache for ``bootstrap_p_sweep_simple`` outputs.

Keyed by ``(param_key, eta_target, sample_idx)`` (which locates the underlying
``(M, v_pop)`` pair via ``eta_pool_cache.sample_dir``) plus a deterministic
``sweep_key`` derived from the sweep config.

Layout, anchored under each pool sample dir::

    src/cache/eta_pool/<param_key>/eta<TT>/sample_NNNN/
        M.npz, fiedler_ref.npz, metadata.json, .complete   # eta_pool entries
        sweeps/
            <sweep_key>/
                result.json     # {"config": ..., "result": ..., "saved_at": ...}
                .complete       # atomic-write sentinel

A sweep entry is written atomically: payload first via temp-file + rename,
then a ``.complete`` sentinel touch. ``load_sweep_result`` returns ``None``
unless the sentinel exists.

Schema version is embedded in the config dict (``schema_version``) so a
behaviour change in ``bootstrap_p_sweep_simple`` can be reflected by bumping
it — old keys then re-derive into different sweep dirs and old caches stay
inert until cleaned up.
"""
from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from .eta_pool_cache import sample_dir
from ..runners.p_sweep_inner import bootstrap_p_sweep_simple


SCHEMA_VERSION = "v1"


def sweep_cache_dir(
    cache_root: Path, param_key: str, eta_target: int, idx: int
) -> Path:
    """``sweeps/`` subdir under the pool sample dir."""
    return sample_dir(cache_root, param_key, eta_target, idx) / "sweeps"


def _build_config(
    p_values: List[float],
    bootstrap_reps: int,
    seed: int,
    num_gaps: int,
    min_split: int,
    early_stop_consecutive_100: int,
    partition_method: str,
) -> Dict[str, Any]:
    return {
        "p_values": [round(float(p), 9) for p in p_values],
        "bootstrap_reps": int(bootstrap_reps),
        "seed": int(seed),
        "num_gaps": int(num_gaps),
        "min_split": int(min_split),
        "early_stop_consecutive_100": int(early_stop_consecutive_100),
        "partition_method": str(partition_method),
        "schema_version": SCHEMA_VERSION,
    }


def compute_sweep_key(
    p_values: List[float],
    bootstrap_reps: int,
    seed: int,
    num_gaps: int,
    min_split: int,
    early_stop_consecutive_100: int,
    partition_method: str,
) -> Tuple[str, Dict[str, Any]]:
    """Return ``(sweep_key, config)`` where sweep_key is sha1[:12] of the config."""
    config = _build_config(
        p_values, bootstrap_reps, seed, num_gaps, min_split,
        early_stop_consecutive_100, partition_method,
    )
    blob = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha1(blob).hexdigest()[:12], config


def load_sweep_result(
    cache_root: Path, param_key: str, eta_target: int, idx: int, sweep_key: str
) -> Optional[Dict[str, Any]]:
    """Return the deserialized ``result`` dict, or ``None`` on cache miss."""
    d = sweep_cache_dir(cache_root, param_key, eta_target, idx) / sweep_key
    if not (d / ".complete").exists():
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
    """Atomically persist (config, result) under the sweep dir."""
    d = sweep_cache_dir(cache_root, param_key, eta_target, idx) / sweep_key
    d.mkdir(parents=True, exist_ok=True)
    payload = {
        "config": config,
        "result": result,
        "saved_at": datetime.now(timezone.utc).isoformat(),
    }
    tmp = d / "result.json.tmp"
    with open(tmp, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    os.replace(tmp, d / "result.json")
    (d / ".complete").touch()
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
    use_cache: bool = True,
) -> Optional[Tuple[Dict[str, Any], bool]]:
    """Return ``(result, was_cached)`` for one sweep, or ``None`` if M_loader returns None.

    ``M_loader`` is a zero-arg callable returning ``(M, v_pop)`` or ``None``;
    on a cache hit it is never invoked, so the heavy ``M.npz`` load is skipped.
    """
    sweep_key, config = compute_sweep_key(
        p_values=p_values, bootstrap_reps=bootstrap_reps, seed=seed,
        num_gaps=num_gaps, min_split=min_split,
        early_stop_consecutive_100=early_stop_consecutive_100,
        partition_method=partition_method,
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
    )
    save_sweep_result(
        cache_root, param_key, eta_target, idx, sweep_key, config, result,
    )
    return result, False
