"""Disk cache for the eta-binned matrix pool.

Layout under ``cache_root / "eta_pool" / <param_key> /``:

    eta01/sample_0001/{M.npz, fiedler_ref.npz, metadata.json, .complete}
    eta05/sample_0001/...
    eta10/...
    eta15/...
    manifest.json

A pool entry is written atomically: payload first, then a ``.complete``
sentinel. Loaders treat a missing sentinel as a miss, so an interrupted run
won't leave a half-written sample lying around. The manifest is updated
after each successful save and persisted via temp-file + rename.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np


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
) -> str:
    """Deterministic key identifying one (n, L, mu, model, pop_size) configuration."""
    return (
        f"n{_fmt_value(int(n))}"
        f"_L{_fmt_value(int(seq_len))}"
        f"_mu{_fmt_value(float(mu))}"
        f"_{tree_model}"
        f"_pop{_fmt_value(float(pop_size))}"
        f"_{seq_model}"
    )


def bin_name(eta_target: int) -> str:
    """Two-digit bin label, e.g. eta_target=5 -> 'eta05'."""
    return f"eta{int(eta_target):02d}"


def pool_root(cache_root: Path) -> Path:
    return Path(cache_root) / "eta_pool"


def param_dir(cache_root: Path, key: str) -> Path:
    return pool_root(cache_root) / key


def bin_dir(cache_root: Path, key: str, eta_target: int) -> Path:
    return param_dir(cache_root, key) / bin_name(eta_target)


def sample_dir(cache_root: Path, key: str, eta_target: int, idx: int) -> Path:
    return bin_dir(cache_root, key, eta_target) / f"sample_{int(idx):04d}"


def count_completed_samples(cache_root: Path, key: str, eta_target: int) -> int:
    """Number of sample subdirs in this bin that have a .complete sentinel."""
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return 0
    return sum(1 for child in d.iterdir() if child.is_dir() and (child / ".complete").exists())


def next_sample_idx(cache_root: Path, key: str, eta_target: int) -> int:
    """Lowest unused 1-based sample index for this bin."""
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return 1
    used = set()
    for child in d.iterdir():
        if not child.is_dir():
            continue
        name = child.name
        if name.startswith("sample_"):
            try:
                used.add(int(name[len("sample_"):]))
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
) -> Path:
    """Atomically persist (M, v_pop, metadata) under sample_dir; return that dir."""
    d = sample_dir(cache_root, key, eta_target, idx)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "M.npz", M=M)
    np.savez_compressed(d / "fiedler_ref.npz", fiedler_ref=v_pop)
    with open(d / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    (d / ".complete").touch()
    return d


def load_pool_entry(
    cache_root: Path,
    key: str,
    eta_target: int,
    idx: int,
) -> Optional[Tuple[np.ndarray, np.ndarray, Dict[str, Any]]]:
    """Inverse of save_pool_entry. Returns None if .complete sentinel missing."""
    d = sample_dir(cache_root, key, eta_target, idx)
    if not (d / ".complete").exists():
        return None
    M = np.load(d / "M.npz")["M"]
    v_pop = np.load(d / "fiedler_ref.npz")["fiedler_ref"]
    with open(d / "metadata.json", "r") as f:
        metadata = json.load(f)
    return M, v_pop, metadata


def list_completed_samples(
    cache_root: Path, key: str, eta_target: int
) -> List[int]:
    """1-based indices of completed samples in this bin (sorted ascending)."""
    d = bin_dir(cache_root, key, eta_target)
    if not d.exists():
        return []
    out: List[int] = []
    for child in sorted(d.iterdir()):
        if child.is_dir() and child.name.startswith("sample_") and (child / ".complete").exists():
            try:
                out.append(int(child.name[len("sample_"):]))
            except ValueError:
                pass
    return sorted(out)


def manifest_path(cache_root: Path, key: str) -> Path:
    return param_dir(cache_root, key) / "manifest.json"


def load_manifest(cache_root: Path, key: str) -> Optional[Dict[str, Any]]:
    p = manifest_path(cache_root, key)
    if not p.exists():
        return None
    try:
        with open(p, "r") as f:
            return json.load(f)
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
    """Build a fresh manifest skeleton for a brand-new param_key."""
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
    """Recount completed samples per bin from disk (ground truth, not manifest)."""
    return {int(t): count_completed_samples(cache_root, key, int(t)) for t in eta_targets}


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
    """Return the target whose |eta - target| is smallest and within tol; else None."""
    if not eta_targets:
        return None
    diffs = [(abs(eta - float(t)), int(t)) for t in eta_targets]
    diffs.sort()
    best_diff, best_t = diffs[0]
    if best_diff <= float(eta_tol):
        return best_t
    return None
