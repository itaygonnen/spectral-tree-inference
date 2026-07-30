"""B-matrix clan-partition sweep on a SYNTHETIC, η-controlled flat-CBM tree.

Parallel of ``bpart_eta_pool`` but with **no simulation**: the two-clan
structure is built directly as a flat constant-block-model distance matrix

    D[i, j] = d_in   if i, j in the same clan      (d_in = -log S_in)
    D[i, j] = d_out  if i, j in different clans     (d_out = -log S_out)
    D[i, i] = 0

with clan sizes ``(n1, n2)`` chosen so the top-split imbalance
``η = max(n1, n2) / min(n1, n2)`` is set *exactly* by the caller. This is the
distance-space analogue of the synthetic flat-CBM similarity used in
``analysis/paper/fig02_pstar_synth_cbm.ipynb``
(``D = -log S``, JC convention), so the reference split ``v_ref`` is exact (the
two true clans) and any degradation in agreement is purely a sub-sampling
artifact at controlled η — no topology-resolution confound.

The sweep + its disk cache are the shared ``utils.bpart_sweep_cache`` mechanism
(keyed per ``(n, η, s_in, s_out)`` matrix), so re-running is an instant cache
hit. Writes ``bpart_eta_pool_n{n}.json`` in the SAME schema as ``bpart_eta_pool``
so ``src.plots.plot_bpart_eta_grid`` renders synthetic and Kingman sweeps identically.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from ..cache_io import make_key
from ..utils.bpart_sweep_cache import aggregate_per_p, compute_or_load_bpart_sweep
from ..utils.logging import log_info
from .bpart_eta_pool import _DEFAULT_P
from .nj_sweep import _SIM_FLOOR


def _clans(n: int, eta: float) -> Tuple[int, int]:
    """Split ``n`` leaves into ``(n1, n2)`` with ``max/min`` closest to ``eta``."""
    n1 = int(round(n / (1.0 + float(eta))))
    n1 = max(1, min(n - 1, n1))
    n2 = n - n1
    return (n1, n2) if n1 <= n2 else (n2, n1)


def _synthetic_clan_D(n1: int, n2: int, s_in: float, s_out: float) -> np.ndarray:
    """Flat two-clan JC distance matrix ``D = -log S`` (diag 0)."""
    d_in = -np.log(max(float(s_in), _SIM_FLOOR))
    d_out = -np.log(max(float(s_out), _SIM_FLOOR))
    n = n1 + n2
    D = np.full((n, n), d_out, dtype=np.float64)
    D[:n1, :n1] = d_in
    D[n1:, n1:] = d_in
    np.fill_diagonal(D, 0.0)
    return D


def bpart_synthetic_for_n(
    n: int, run_dir: str, *,
    eta_targets: Sequence[float] = (1, 5),
    s_in: float = 0.9, s_out: float = 0.05,
    reps: int = 5,
    p_values: Sequence[float] = tuple(_DEFAULT_P),
    seed_base: int = 0,
) -> Dict:
    """Synthetic flat-CBM B-method subsampling sweep at fixed η bins for one ``n``.

    For each η target the clans ``(n1, n2)`` are sized to realize that η, the
    exact distance matrix is built, and the (cached) B-method sweep over
    ``p_values`` × ``reps`` measures agreement / dot / ARI / NMI against the
    full-data reference split. Writes ``bpart_eta_pool_n{n}.json``.
    """
    p_values = list(p_values)

    bins: Dict[str, Dict] = {}
    for target in eta_targets:
        n1, n2 = _clans(n, float(target))
        realized_eta = max(n1, n2) / min(n1, n2)
        log_info("bpart_synth",
                 f"  n={n} eta~{target}: n1={n1} n2={n2} "
                 f"(eta={realized_eta:.2f}) x {reps} reps", force=True)

        identity = make_key("synthetic_flat_cbm", n=n, eta=float(target),
                            s_in=s_in, s_out=s_out)
        outcome = compute_or_load_bpart_sweep(
            identity,
            D_loader=lambda n1=n1, n2=n2: _synthetic_clan_D(n1, n2, s_in, s_out),
            p_values=p_values, reps=reps,
            seed_base=seed_base, imputation="mean",
        )
        result, was_cached = outcome
        log_info("bpart_synth", f"    cache {'HIT' if was_cached else 'MISS'}", force=True)

        bins[str(int(round(target)))] = {
            "eta_target": float(target),
            "eta_mean": float(realized_eta), "eta_std": 0.0,
            "n_samples": 1, "reps": reps, "n1": n1, "n2": n2,
            "per_p": aggregate_per_p([result["per_p"]]),
        }

    out = {"n": n, "s_in": s_in, "s_out": s_out,
           "tree_model": "synthetic_flat_cbm", "bins": bins}
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(run_dir) / f"bpart_eta_pool_n{n}.json"
    out_path.write_text(json.dumps(out, indent=2))
    log_info("bpart_synth", f"Wrote {out_path}", force=True)
    return out
