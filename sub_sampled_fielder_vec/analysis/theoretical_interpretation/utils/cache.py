"""Two-level disk cache for theoretical-interpretation experiments.

Layout (under ``cache_dir``):

    full/<full_key>/
        S.npz, fiedler_full.npz, metadata.json, .complete
    subsampled/<full_key>/p<p>_seed<seed>/
        fiedler_hat.npz, agreement.json, .complete

Atomicity: every write touches the ``.complete`` sentinel last. Loaders return
``None`` when the sentinel is missing, so a half-written cache entry behaves
like a miss and is re-computed on the next call.

The cache layer is fully topology-agnostic: callers supply the ``full_key`` and
a builder closure. Only ``S_hat`` is intentionally NOT persisted — it is cheap
to regenerate from full ``S`` and the ``(p, seed)`` and would otherwise blow up
disk usage by orders of magnitude on the production sweep.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Callable, Dict, Literal, Optional, Tuple

import numpy as np


def _fmt_value(value: Any) -> str:
    """Filesystem-safe stringification: dots become 'p', floats get 4 decimals."""
    if isinstance(value, float):
        return f"{value:.4f}".replace(".", "p")
    if isinstance(value, int):
        return f"{value:04d}"
    return str(value).replace(".", "p")


def make_full_key(prefix: str, **kwargs: Any) -> str:
    """Deterministic, filesystem-safe key from a prefix and sorted kwargs.

    Example: ``make_full_key("balanced_binary", alpha=0.9, n=512)`` returns
    ``"balanced_binary_alpha0p9000_n0512"``.
    """
    parts = [prefix]
    for key in sorted(kwargs):
        parts.append(f"{key}{_fmt_value(kwargs[key])}")
    return "_".join(parts)


def make_subsample_key(p: float, seed: int) -> str:
    """Per-trial key, e.g. ``"p0p0100_seed0002"``."""
    return f"p{_fmt_value(float(p))}_seed{_fmt_value(int(seed))}"


def _full_dir(cache_dir: Path, full_key: str) -> Path:
    return Path(cache_dir) / "full" / full_key


def _subsample_dir(cache_dir: Path, full_key: str, p: float, seed: int) -> Path:
    return Path(cache_dir) / "subsampled" / full_key / make_subsample_key(p, seed)


def load_full(cache_dir: Path, full_key: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    d = _full_dir(cache_dir, full_key)
    if not (d / ".complete").exists():
        return None
    try:
        S = np.load(d / "S.npz")["S"]
        fiedler = np.load(d / "fiedler_full.npz")["fiedler_full"]
    except Exception:
        return None
    return S, fiedler


def save_full(
    cache_dir: Path,
    full_key: str,
    S: np.ndarray,
    fiedler_full: np.ndarray,
    metadata: Dict[str, Any],
) -> None:
    d = _full_dir(cache_dir, full_key)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "S.npz", S=S)
    np.savez_compressed(d / "fiedler_full.npz", fiedler_full=fiedler_full)
    with open(d / "metadata.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    (d / ".complete").touch()


def get_or_compute_full(
    cache_dir: Path,
    full_key: str,
    builder_fn: Callable[[], Tuple[np.ndarray, np.ndarray, Dict[str, Any]]],
    force: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(S, fiedler_full)`` from cache or compute via ``builder_fn``."""
    if not force:
        hit = load_full(cache_dir, full_key)
        if hit is not None:
            return hit
    S, fiedler_full, metadata = builder_fn()
    save_full(cache_dir, full_key, S, fiedler_full, metadata)
    return S, fiedler_full


def load_subsample(
    cache_dir: Path, full_key: str, p: float, seed: int
) -> Optional[Tuple[np.ndarray, float]]:
    d = _subsample_dir(cache_dir, full_key, p, seed)
    if not (d / ".complete").exists():
        return None
    try:
        fiedler_hat = np.load(d / "fiedler_hat.npz")["fiedler_hat"]
        with open(d / "agreement.json", "r") as f:
            agreement = float(json.load(f)["agreement"])
    except Exception:
        return None
    return fiedler_hat, agreement


def save_subsample(
    cache_dir: Path,
    full_key: str,
    p: float,
    seed: int,
    fiedler_hat: np.ndarray,
    agreement: float,
) -> None:
    d = _subsample_dir(cache_dir, full_key, p, seed)
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "fiedler_hat.npz", fiedler_hat=fiedler_hat)
    with open(d / "agreement.json", "w") as f:
        json.dump({"p": float(p), "seed": int(seed), "agreement": float(agreement)}, f)
    (d / ".complete").touch()


def get_or_compute_subsample(
    cache_dir: Path,
    full_key: str,
    p: float,
    seed: int,
    compute_fn: Callable[[], Tuple[np.ndarray, float]],
    force: bool = False,
) -> Tuple[np.ndarray, float]:
    """Return ``(fiedler_hat, agreement)`` from cache or compute via ``compute_fn``."""
    if not force:
        hit = load_subsample(cache_dir, full_key, p, seed)
        if hit is not None:
            return hit
    fiedler_hat, agreement = compute_fn()
    save_subsample(cache_dir, full_key, p, seed, fiedler_hat, float(agreement))
    return fiedler_hat, float(agreement)


def clear_theoretical_cache(
    cache_dir: Path, scope: Literal["all", "full", "subsampled"] = "all"
) -> None:
    """Remove cached entries. ``scope`` selects which level to wipe."""
    cache_dir = Path(cache_dir)
    targets = []
    if scope in ("all", "full"):
        targets.append(cache_dir / "full")
    if scope in ("all", "subsampled"):
        targets.append(cache_dir / "subsampled")
    for t in targets:
        if t.exists():
            shutil.rmtree(t)
