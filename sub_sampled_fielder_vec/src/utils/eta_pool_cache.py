"""Disk cache for the eta-binned matrix pool.

Thin shim over :mod:`src.cache_io`. Backing scope is ``pool_sample``;
the ``cache_root`` argument is accepted for backward compat but ignored —
the canonical location is :data:`src.cache_io.CACHE_ROOT`.

Layout under ``cache/pool_sample/<param_key>/``::

    eta01/sample_0001/{M.npz, fiedler_ref.npz, metadata.json, .complete}
    eta05/sample_0001/...
    manifest.json

The manifest sits next to the bins and is updated via temp-file + rename.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

from ..cache_io import pool_sample as _scope, SENTINEL


def _fmt_value(value: Any) -> str:
    """Filesystem-safe stringification: dots become 'p', floats get 4 decimals."""
    if isinstance(value, float):
        return f"{value:.4f}".replace(".", "p")
    if isinstance(value, int):
        return f"{value:04d}"
    return str(value).replace(".", "p")


def param_key(
    n: int,
    seq_len: int,
    mu: float,
    tree_model: str,
    pop_size: float,
    seq_model: str = "JC69",
    matrix_kind: str = "similarity",
    distance_alpha: float = 1.0,
) -> str:
    """Deterministic key identifying one pool configuration."""
    base = (
        f"n{_fmt_value(int(n))}"
        f"_L{_fmt_value(int(seq_len))}"
        f"_mu{_fmt_value(float(mu))}"
        f"_{tree_model}"
        f"_pop{_fmt_value(float(pop_size))}"
        f"_{seq_model}"
    )
    if matrix_kind == "distance":
        return base + f"_dist_a{_fmt_value(float(distance_alpha))}"
    return base


def bin_name(eta_target: int) -> str:
    return f"eta{int(eta_target):02d}"


# ----- paths --------------------------------------------------------------


def pool_root(cache_root: Path) -> Path:
    """Legacy alias. Returns the canonical scope root regardless of input."""
    return _scope.root


def param_dir(cache_root: Path, key: str) -> Path:
    return _scope.path(key)


def bin_dir(cache_root: Path, key: str, eta_target: int) -> Path:
    return _scope.path(key, bin_name(eta_target))


def sample_dir(cache_root: Path, key: str, eta_target: int, idx: int) -> Path:
    return _scope.path(key, bin_name(eta_target), f"sample_{int(idx):04d}")


# ----- per-sample I/O -----------------------------------------------------


def count_completed_samples(cache_root: Path, key: str, eta_target: int) -> int:
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return 0
    return sum(
        1 for child in d.iterdir()
        if child.is_dir() and (child / SENTINEL).exists()
    )


def next_sample_idx(cache_root: Path, key: str, eta_target: int) -> int:
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return 1
    used: set[int] = set()
    for child in d.iterdir():
        if child.is_dir() and child.name.startswith("sample_"):
            try:
                used.add(int(child.name[len("sample_"):]))
            except ValueError:
                pass
    idx = 1
    while idx in used:
        idx += 1
    return idx


def save_pool_entry(
    cache_root: Path,
    key: str,
    eta_target: int,
    idx: int,
    M: np.ndarray,
    v_pop: np.ndarray,
    metadata: Dict[str, Any],
    partition: Optional[np.ndarray] = None,
    tree_newick: Optional[str] = None,
) -> Path:
    """Atomically persist (M, v_pop, partition, tree, metadata) under sample_dir.

    ``partition`` is a boolean array of length n representing the Fiedler-sign
    partition (``v_pop > 0``). Stored so downstream readers don't have to
    recompute the sign convention.

    ``tree_newick`` is the actual tree used at sample time, serialized as
    Newick. Stored because ``np.random.seed`` does NOT make dendropy's Kingman
    tree generation deterministic (dendropy uses its own ``GLOBAL_RNG``), so
    the seed alone is insufficient to reproduce the tree downstream.

    Both args optional for backwards compat; new callers should pass them.
    """
    artifacts: Dict[str, Any] = {"M": M, "fiedler_ref": v_pop, "metadata": metadata}
    if partition is not None:
        artifacts["partition"] = np.asarray(partition, dtype=bool)
    if tree_newick is not None:
        artifacts["tree"] = str(tree_newick)
    return _scope.save(
        key, bin_name(eta_target), f"sample_{int(idx):04d}",
        **artifacts,
    )


def load_pool_entry(
    cache_root: Path,
    key: str,
    eta_target: int,
    idx: int,
) -> Optional[Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[str], Dict[str, Any]]]:
    """Return ``(M, v_pop, partition, tree_newick, metadata)`` or ``None`` on cache miss.

    ``partition`` and ``tree_newick`` are ``None`` for legacy samples that
    pre-date the partition + tree persistence changes. Callers that need them
    strictly should assert.
    """
    data = _scope.try_load(key, bin_name(eta_target), f"sample_{int(idx):04d}")
    if data is None:
        return None
    try:
        M = data["M"]
        v_pop = data["fiedler_ref"]
        metadata = data["metadata"]
    except KeyError:
        return None
    partition = data.get("partition")
    if partition is not None:
        partition = np.asarray(partition, dtype=bool)
    tree_newick = data.get("tree")
    if tree_newick is not None:
        tree_newick = str(tree_newick)
    return M, v_pop, partition, tree_newick, metadata


def list_completed_samples(
    cache_root: Path, key: str, eta_target: int
) -> List[int]:
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return []
    out: List[int] = []
    for child in sorted(d.iterdir()):
        if (child.is_dir() and child.name.startswith("sample_")
                and (child / SENTINEL).exists()):
            try:
                out.append(int(child.name[len("sample_"):]))
            except ValueError:
                pass
    return sorted(out)


# ----- manifest -----------------------------------------------------------


def manifest_path(cache_root: Path, key: str) -> Path:
    return _scope.path(key) / "manifest.json"


def load_manifest(cache_root: Path, key: str) -> Optional[Dict[str, Any]]:
    p = manifest_path(cache_root, key)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def save_manifest(cache_root: Path, key: str, manifest: Dict[str, Any]) -> None:
    """Write manifest atomically via temp file + rename."""
    p = manifest_path(cache_root, key)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".json.tmp")
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    os.replace(tmp, p)


def init_manifest(
    key: str,
    params: Dict[str, Any],
    eta_targets: Iterable[int],
    eta_tol: float,
    samples_per_bin: int,
) -> Dict[str, Any]:
    targets = [int(t) for t in eta_targets]
    return {
        "param_key": key,
        "params": dict(params),
        "eta_targets": targets,
        "eta_tol": float(eta_tol),
        "samples_per_bin": int(samples_per_bin),
        "bin_counts": {str(t): 0 for t in targets},
        "seeds_per_bin": {str(t): [] for t in targets},
        "attempts_used": 0,
        "last_seed": -1,
    }


def bin_counts_from_disk(
    cache_root: Path, key: str, eta_targets: Iterable[int]
) -> Dict[int, int]:
    return {
        int(t): count_completed_samples(cache_root, key, int(t))
        for t in eta_targets
    }


def all_bins_full(manifest: Dict[str, Any]) -> bool:
    target = int(manifest["samples_per_bin"])
    return all(int(c) >= target for c in manifest["bin_counts"].values())


def remaining_capacity(manifest: Dict[str, Any], eta_target: int) -> int:
    target = int(manifest["samples_per_bin"])
    have = int(manifest["bin_counts"].get(str(int(eta_target)), 0))
    return max(0, target - have)


def closest_target(
    eta: float, eta_targets: List[int], eta_tol: float
) -> Optional[int]:
    if not eta_targets:
        return None
    diffs = [(abs(eta - float(t)), int(t)) for t in eta_targets]
    diffs.sort()
    best_diff, best_t = diffs[0]
    if best_diff <= float(eta_tol):
        return best_t
    return None
