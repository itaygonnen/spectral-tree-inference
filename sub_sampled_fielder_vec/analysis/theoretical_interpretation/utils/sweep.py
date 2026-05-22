"""Topology-agnostic (p, seed) sweep driver for a single tree size.

The outer loop over n is intentionally left in the notebook so callers
can plug in different similarity-matrix builders for non-balanced
topologies.

This driver is metric-agnostic. Pass a ``metrics`` dict mapping a name
to a function ``(v_full, v_hat) -> float`` and the sweep persists every
named metric inside a single per-trial ``metrics.json``. On a cache hit,
metrics that are already present are reused; metrics that are missing
are computed from the cached Fiedler vector and merged in — **no
resampling required to add a new metric to a completed sweep**.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable, Mapping, Optional

import numpy as np
import pandas as pd

from .cache import make_subsample_key
from .recovery import compute_recovery
from .spectral import compute_fiedler_of_S, subsample_S

# cache_io is two parents up — bootstrap path the same way the notebooks do.
import sys
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
from src.cache_io import sweep_trial as _scope  # noqa: E402


Metric = Callable[[np.ndarray, np.ndarray], float]
DEFAULT_METRICS: dict[str, Metric] = {"sign_agreement": compute_recovery}


def run_sweep(
    cache_dir: Path,
    full_key: str,
    S: np.ndarray,
    fiedler_full: np.ndarray,
    p_grid: Iterable[float],
    n_trials: int,
    metrics: Optional[Mapping[str, Metric]] = None,
    progress: bool = True,
) -> pd.DataFrame:
    """Run a (p × trial) sub-sampling sweep against a single full S.

    For each ``(p, seed)`` pair:
      1. ``try_load`` the trial. If cached, take ``fiedler_hat`` and any
         already-computed metrics from disk.
      2. If miss, subsample S → ``S_hat`` → ``fiedler_hat``.
      3. For each ``(name, fn)`` in ``metrics`` not yet on disk, compute
         ``fn(fiedler_full, fiedler_hat)`` and merge into ``metrics.json``
         (re-stamping the ``.complete`` sentinel).

    Degenerate-sample guard: when ``p * n < 1`` the Erdős-Rényi subgraph
    is subcritical and the Fiedler-on-Laplacian computation latches onto
    near-basis vectors. The ``sign_agreement`` metric is forced to 0.5 in
    that regime. Other metrics are computed honestly — caller decides
    how to handle them.

    Returns a long-form DataFrame with one row per ``(p, seed)`` and one
    column per metric name. The legacy ``"agreement"`` column has been
    replaced by ``"sign_agreement"``.
    """
    p_grid = np.asarray(list(p_grid), dtype=float)
    metrics_dict: dict[str, Metric] = dict(metrics) if metrics else dict(DEFAULT_METRICS)
    metric_names = list(metrics_dict)
    n = S.shape[0]
    rows: list[dict] = []
    total = len(p_grid) * n_trials
    done = 0

    for p in p_grid:
        for seed in range(n_trials):
            p_f = float(p)
            seed_i = int(seed)
            sub_key = make_subsample_key(p_f, seed_i)

            # 1. Try cache.
            cached = _scope.try_load(full_key, sub_key)
            if cached is not None and "fiedler_hat" in cached:
                fiedler_hat = cached["fiedler_hat"]
                cached_metrics = dict(cached.get("metrics") or {})
            else:
                S_hat = subsample_S(S, p_f, seed_i)
                fiedler_hat = compute_fiedler_of_S(S_hat, sampling_prob=p_f)
                cached_metrics = {}
                # Persist substrate immediately so a crash mid-metric leaves
                # us cleanly resumable.
                _scope.save(
                    full_key, sub_key,
                    fiedler_hat=fiedler_hat,
                    metrics=cached_metrics,
                )

            # 2. Fill in missing metrics.
            missing = [name for name in metric_names if name not in cached_metrics]
            if missing:
                for name in missing:
                    fn = metrics_dict[name]
                    cached_metrics[name] = float(fn(fiedler_full, fiedler_hat))
                _scope.extend(full_key, sub_key, metrics=cached_metrics)

            row = {"full_key": full_key, "p": p_f, "seed": seed_i}
            for name in metric_names:
                value = cached_metrics[name]
                if name == "sign_agreement" and p_f * n < 1.0:
                    # Same floor as the legacy code: cached value is the raw
                    # compute_recovery output (can be <0.5 in degenerate
                    # regime); the DataFrame reports chance level instead.
                    value = 0.5
                row[name] = value
            rows.append(row)

            done += 1
            if progress and done % max(1, total // 20) == 0:
                print(f"  [{full_key}] {done}/{total} trials done")

    return pd.DataFrame(rows)


def load_or_extend_metrics(
    full_key: str,
    metrics: Mapping[str, Metric],
    fiedler_full: np.ndarray,
    p_grid: Optional[Iterable[float]] = None,
    n_trials: Optional[int] = None,
) -> pd.DataFrame:
    """Add metrics to an already-complete sweep without resampling.

    Iterates the cached ``sweep_trial`` entries under ``full_key``,
    computes any missing metrics from the stored Fiedler vector, merges
    them into ``metrics.json``, and returns a long-form DataFrame with
    one row per ``(p, seed)``.

    If ``p_grid`` / ``n_trials`` are provided, only those trials are
    iterated (and missing ones are skipped); otherwise every cached
    trial under ``full_key`` is processed.
    """
    metrics_dict: dict[str, Metric] = dict(metrics)
    metric_names = list(metrics_dict)
    parent = _scope.path(full_key)
    if not parent.exists():
        return pd.DataFrame(columns=["full_key", "p", "seed", *metric_names])

    rows: list[dict] = []

    if p_grid is not None and n_trials is not None:
        trial_keys = [
            make_subsample_key(float(p), int(seed))
            for p in p_grid for seed in range(int(n_trials))
        ]
    else:
        trial_keys = sorted(
            d.name for d in parent.iterdir()
            if d.is_dir() and (d / ".complete").exists()
        )

    for sub_key in trial_keys:
        cached = _scope.try_load(full_key, sub_key)
        if cached is None or "fiedler_hat" not in cached:
            continue
        fiedler_hat = cached["fiedler_hat"]
        cached_metrics = dict(cached.get("metrics") or {})
        missing = [name for name in metric_names if name not in cached_metrics]
        if missing:
            for name in missing:
                cached_metrics[name] = float(metrics_dict[name](fiedler_full, fiedler_hat))
            _scope.extend(full_key, sub_key, metrics=cached_metrics)
        # Parse p, seed from the sub_key (format "p<...>_seed<...>").
        try:
            p_part, seed_part = sub_key.split("_seed")
            p_val = float(p_part[1:].replace("p", "."))
            seed_val = int(seed_part)
        except Exception:
            continue
        row = {"full_key": full_key, "p": p_val, "seed": seed_val}
        for name in metric_names:
            row[name] = cached_metrics.get(name)
        rows.append(row)

    return pd.DataFrame(rows)
