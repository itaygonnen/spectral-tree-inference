"""Topology-agnostic (p, seed) sweep driver for a single tree size.

The outer loop over n is intentionally left in the notebook so that callers
can plug in different similarity-matrix builders for non-balanced topologies.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .cache import get_or_compute_subsample
from .recovery import compute_recovery
from .spectral import compute_fiedler_of_S, subsample_S


def run_sweep(
    cache_dir: Path,
    full_key: str,
    S: np.ndarray,
    fiedler_full: np.ndarray,
    p_grid: Iterable[float],
    n_trials: int,
    progress: bool = True,
) -> pd.DataFrame:
    """Run a (p × trial) sub-sampling sweep against a single full S.

    For each ``(p, seed)`` pair, sub-samples S, recomputes the Fiedler vector,
    measures sign-recovery against ``fiedler_full``, and persists the trial
    result through ``get_or_compute_subsample``. Trials that hit the cache
    return immediately without recomputing.
    """
    p_grid = np.asarray(list(p_grid), dtype=float)
    rows = []
    total = len(p_grid) * n_trials
    done = 0

    for p in p_grid:
        for seed in range(n_trials):

            def compute_trial(p=float(p), seed=int(seed)):
                S_hat = subsample_S(S, p, seed)
                fiedler_hat = compute_fiedler_of_S(S_hat, sampling_prob=p)
                agreement = compute_recovery(fiedler_full, fiedler_hat)
                return fiedler_hat, agreement

            _, agreement = get_or_compute_subsample(
                cache_dir, full_key, float(p), int(seed), compute_trial
            )
            rows.append(
                {"full_key": full_key, "p": float(p), "seed": int(seed), "agreement": agreement}
            )
            done += 1
            if progress and done % max(1, total // 20) == 0:
                print(f"  [{full_key}] {done}/{total} trials done")

    return pd.DataFrame(rows)
