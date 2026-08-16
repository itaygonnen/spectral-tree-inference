"""Collect Griffing-on-D (B = HDH, sign threshold) sub-sampling sweeps over the
eta-binned pool, in the **same shape** as
:func:`src.utils.eta_pool_sweep.run_and_collect_sweeps`.

The Fiedler-on-S methods (``sign``/``sigma2``/``kmeans``) are collected by
``run_and_collect_sweeps``. Griffing-on-D lives in a *different* inner sweep
(:mod:`src.utils.bpart_sweep_cache`) and cache scope (``bpart_sweep``), but runs
on the **same pool samples** -- the distance matrix is recovered losslessly as
``D = -log(M)`` (see :mod:`src.runners.bpart_eta_pool`). This module is the thin
adapter that turns those per-sample bpart sweeps into the
``results[(n, eta)] = ndarray[samples x len(p_values)]`` array the notebooks
expect, so Griffing can be plotted as a third "method" alongside the others.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..runners.bpart_eta_pool import _D_from_M
from .bpart_sweep_cache import compute_or_load_bpart_sweep
from .griffing import DEFAULT_SOLVER
from .eta_pool_cache import bin_name, list_completed_samples, load_pool_entry, param_key
from .eta_pool_sweep import (
    ETA_SAMPLE_CAP, ETA_TARGETS, P_VALUES, _default_pool_params,
)

_METRICS = ("agreement", "dot", "ari", "nmi")


def collect_griffing_sweeps(
    ns: Sequence[int],
    eta_targets: Sequence[int] = ETA_TARGETS,
    *,
    p_values: Sequence[float] = P_VALUES,
    reps: int = 10,
    seed_base: int = 0,
    pool_params: Optional[Dict[str, object]] = None,
    max_per_eta: Optional[int] = None,
    use_cache: bool = True,
    metric: str = "nmi",
    eigsolver: str = DEFAULT_SOLVER,
    on_sample: Optional[Callable[[int, int, int, bool], None]] = None,
) -> Dict[Tuple[int, int], np.ndarray]:
    """Run (or load) the B-method (Griffing-on-D) sweep for every pool sample.

    Returns ``results[(n, eta)]`` = ndarray of shape
    ``(num_samples x len(p_values))`` -- one row per pool sample, each row the
    sample's per-p ``metric`` curve averaged over reps. Mirrors one method's
    sub-dict from ``run_and_collect_sweeps`` so the caller can splice it in as a
    third operator.

    ``eigsolver`` picks how the single leading eigenpair of ``B`` is computed;
    ``"lm_k1"`` is what makes a 25-point grid over the whole pool tractable (see
    :data:`src.utils.griffing.SOLVERS`). ``on_sample(n, eta, idx, was_cached)``
    fires per sample so a long run can report progress.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {_METRICS}, got {metric!r}")
    pp = pool_params if pool_params is not None else _default_pool_params()
    p_values = list(p_values)

    results: Dict[Tuple[int, int], np.ndarray] = {}
    for n in ns:
        key = param_key(n=n, **pp)
        for eta in eta_targets:
            cap = max_per_eta if max_per_eta is not None else ETA_SAMPLE_CAP.get(eta, 10)
            idxs = list_completed_samples(None, key, int(eta))[:cap]
            curves: List[List[float]] = []
            for idx in idxs:
                entry = load_pool_entry(None, key, int(eta), idx)
                if entry is None:
                    continue
                M = entry[0]
                outcome = compute_or_load_bpart_sweep(
                    key, bin_name(int(eta)), f"sample_{int(idx):04d}",
                    D_loader=lambda M=M: _D_from_M(M),
                    p_values=p_values, reps=reps, seed_base=seed_base,
                    imputation="mean", eigsolver=eigsolver, use_cache=use_cache,
                    # Figure 5's twin is Figure 3, which k-means-rounds the bootstrap
                    # AVERAGE of the Fiedler vectors. Scoring each replicate and
                    # averaging the metric (the historical default, and the right
                    # accounting for Figure 4 against Figure 2) is a different
                    # estimator, and paper_figures.py records the resulting mismatch
                    # as "NOT like-for-like with Figure 3". This makes it like-for-like.
                    aggregation="avg_vector",
                )
                if outcome is None:
                    continue
                result, was_cached = outcome
                if on_sample is not None:
                    on_sample(n, int(eta), int(idx), was_cached)
                curves.append([float(np.mean(pp_[metric])) for pp_ in result["per_p"]])
            results[(n, int(eta))] = (
                np.asarray(curves) if curves else np.empty((0, len(p_values)))
            )
    return results
