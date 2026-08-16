"""Collect Griffing-on-D (B = HDH, sign threshold) sub-sampling sweeps over a
synthetic flat-CBM ``(eta, m)`` grid, in the **same shape** as
:func:`src.utils.eta_pool_griffing.collect_griffing_sweeps`.

This is the distance-route twin of the sweep loop in
``analysis/sec5_empirical/synthesized/nonbalanced_flat_cbm.ipynb`` (the paper's
Fig 2), and the deliberate mirror of ``eta_pool_griffing`` on the generated side:
one function, one returned dict ``{(m, eta): ndarray[trials x len(p_grid)]}``,
which is exactly what ``analysis.utils.sweep_plots_two_panel.plot_pstar_pair``
consumes. Nothing here plots.

No simulation and no estimation noise: the two-clan structure is written down
directly as a flat CBM in *distance* space,

    D[i, j] = -log S_in    same clan
    D[i, j] = -log S_out   different clans
    D[i, i] = 0

so the reference split (the sign pattern of the leading eigenvector of
``B = HDH`` on the full ``D``) is the exact pair of true clans, and any loss of
agreement at ``p < 1`` is purely a sub-sampling artifact at controlled ``eta``.

Clan sizes use the notebook's exact ``n1 = m // (1 + eta)``, not
``bpart_synthetic._clans``' ``round(m / (1 + eta))``: ``m`` is chosen divisible by
``1 + eta`` precisely so the realized eta is the requested one, and rounding
would silently shift it. The matrix builder itself
(:func:`src.runners.bpart_synthetic._synthetic_clan_D`) is reused as-is.

The sweep and its disk cache are the shared
:mod:`src.utils.bpart_sweep_cache` mechanism, keyed per ``(m, eta, s_in, s_out)``
matrix plus the sweep config, so runs are resumable and a re-run with fewer
trials just reads the first seeds of an existing cache.
"""
from __future__ import annotations

from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..cache_io import make_key
from ..runners.bpart_synthetic import _synthetic_clan_D
from .bpart_sweep_cache import _METRICS, compute_or_load_bpart_sweep
from .griffing import DEFAULT_SOLVER


def clan_sizes(m: int, eta: int) -> Tuple[int, int]:
    """``(n1, n2)`` for total size ``m`` at imbalance ``eta``, notebook-exact.

    Mirrors ``nonbalanced_flat_cbm.ipynb`` cell 5. Exact only when ``1 + eta``
    divides ``m``; the caller's grid is responsible for that (each eta needs its
    own m grid, since no single grid is divisible by 2, 6, 11 and 16).
    """
    n1 = m // (1 + int(eta))
    return n1, m - n1


def collect_griffing_synth_sweeps(
    m_grids: Mapping[int, Sequence[int]],
    *,
    s_in: float = 0.9,
    s_out: float = 0.05,
    p_grid: Sequence[float],
    n_trials: int = 8,
    seed_base: int = 0,
    imputation: str = "mean",
    eigsolver: str = "lm_k1",
    metric: str = "nmi",
    use_cache: bool = True,
    on_cell: Optional[Callable[[int, int, bool], None]] = None,
) -> Dict[Tuple[int, int], np.ndarray]:
    """Run (or load) the B-method sweep for every ``(eta, m)`` in ``m_grids``.

    ``m_grids`` maps ``eta -> [m, ...]`` (each eta its own grid; see
    :func:`clan_sizes`). Returns ``results[(m, eta)]`` = ndarray of shape
    ``(n_trials x len(p_grid))``, one row per trial. ``n_trials`` is the bpart
    path's ``reps``; its seeds are ``seed_base + 10_000*p_index + rep``, so
    trials are independent within each ``p``.

    This path deliberately keeps ``aggregation="per_rep"`` (the default). Its twin
    is Figure 2, whose similarity-side sweep (``analysis/utils/sweep.run_sweep``)
    also scores ``n_trials`` independent single draws with no inner averaging --
    ``paper_figures.py`` records the pair as "LIKE-FOR-LIKE ... same per-trial
    accounting". Switching to ``"avg_vector"`` here would collapse every row to one
    and break that. Figure 5 is the opposite case; see
    :mod:`src.utils.eta_pool_griffing`.

    ``on_cell(eta, m, was_cached)`` fires once per grid cell for progress.
    """
    if metric not in _METRICS:
        raise ValueError(f"metric must be one of {_METRICS}, got {metric!r}")
    p_grid = [float(p) for p in p_grid]

    results: Dict[Tuple[int, int], np.ndarray] = {}
    for eta in sorted(m_grids):
        for m in sorted(m_grids[eta]):
            n1, n2 = clan_sizes(m, eta)
            if n1 < 1:
                raise ValueError(f"eta={eta} m={m} gives an empty clan (n1={n1})")

            identity = make_key("synthetic_flat_cbm", n=m, eta=float(eta),
                                s_in=s_in, s_out=s_out)
            result, was_cached = compute_or_load_bpart_sweep(
                identity,
                D_loader=lambda n1=n1, n2=n2: _synthetic_clan_D(n1, n2, s_in, s_out),
                p_values=p_grid, reps=n_trials, seed_base=seed_base,
                imputation=imputation, eigsolver=eigsolver, use_cache=use_cache,
            )
            if on_cell is not None:
                on_cell(int(eta), int(m), bool(was_cached))

            # per_p is ordered as p_grid; each entry holds n_trials raw values.
            # Transpose to (trials x p) so a row is one trial's whole curve.
            per_p = [list(entry[metric]) for entry in result["per_p"]]
            results[(int(m), int(eta))] = np.asarray(per_p, dtype=float).T
    return results
