"""Two-level disk cache for theoretical-interpretation experiments.

Thin shim over :mod:`src.cache_io`. The unified layer owns the on-disk
layout, the ``.complete`` sentinel discipline, and the key format. Public
function names here are preserved for backward compat with existing
notebooks; the ``cache_dir`` argument is accepted but ignored — the scope
under :data:`src.cache_io.CACHE_ROOT` is canonical.

Old API → new scope:
    * ``full/<key>/`` → ``cache_io.full_matrix``
    * ``subsampled/<key>/<sub_key>/`` → ``cache_io.sweep_trial``

The ``agreement`` scalar that ``save_subsample``/``load_subsample`` exposes
is internally stored as ``metrics["sign_agreement"]`` so new metrics can
be added alongside without re-sweeping. Use
:func:`load_or_extend_metrics` in ``sweep.py`` for the metric-agnostic
path.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, Literal, Optional, Tuple

import numpy as np

# Local imports go through the project root via the same sys.path bootstrap
# the notebooks already do; cache_io is the single source of truth.
import sys
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.cache_io import (  # noqa: E402
    full_matrix as _full_matrix,
    sweep_trial as _sweep_trial,
    make_key as _make_key_unified,
)


def _fmt_value(value: Any) -> str:
    """Filesystem-safe stringification: dots become 'p', floats get 4 decimals."""
    if isinstance(value, float):
        return f"{value:.4f}".replace(".", "p")
    if isinstance(value, int):
        return f"{value:04d}"
    return str(value).replace(".", "p")


def make_full_key(prefix: str, **kwargs: Any) -> str:
    """Deterministic, filesystem-safe key from a prefix and sorted kwargs."""
    return _make_key_unified(prefix, **kwargs)


def make_subsample_key(p: float, seed: int) -> str:
    """Per-trial key, e.g. ``"p0p0100_seed0002"``."""
    return f"p{_fmt_value(float(p))}_seed{_fmt_value(int(seed))}"


def load_full(cache_dir: Path, full_key: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    data = _full_matrix.try_load(full_key)
    if data is None:
        return None
    try:
        return data["S"], data["fiedler"]
    except KeyError:
        return None


def save_full(
    cache_dir: Path,
    full_key: str,
    S: np.ndarray,
    fiedler_full: np.ndarray,
    metadata: Dict[str, Any],
) -> None:
    _full_matrix.save(full_key, S=S, fiedler=fiedler_full, metadata=metadata)


def get_or_compute_full(
    cache_dir: Path,
    full_key: str,
    builder_fn: Callable[[], Tuple[np.ndarray, np.ndarray, Dict[str, Any]]],
    force: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
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
    sub_key = make_subsample_key(p, seed)
    data = _sweep_trial.try_load(full_key, sub_key)
    if data is None:
        return None
    try:
        fiedler_hat = data["fiedler_hat"]
        metrics = data.get("metrics") or {}
        if "sign_agreement" not in metrics:
            return None
        return fiedler_hat, float(metrics["sign_agreement"])
    except KeyError:
        return None


def save_subsample(
    cache_dir: Path,
    full_key: str,
    p: float,
    seed: int,
    fiedler_hat: np.ndarray,
    agreement: float,
) -> None:
    sub_key = make_subsample_key(p, seed)
    _sweep_trial.save(
        full_key, sub_key,
        fiedler_hat=fiedler_hat,
        metrics={"sign_agreement": float(agreement)},
    )


def get_or_compute_subsample(
    cache_dir: Path,
    full_key: str,
    p: float,
    seed: int,
    compute_fn: Callable[[], Tuple[np.ndarray, float]],
    force: bool = False,
) -> Tuple[np.ndarray, float]:
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
    if scope in ("all", "full"):
        _full_matrix.clear()
    if scope in ("all", "subsampled"):
        _sweep_trial.clear()
