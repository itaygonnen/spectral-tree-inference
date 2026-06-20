"""B-matrix clan-partition sweep on the η-binned Kingman pool.

The plain ``bpart_sweep`` lets the top-split imbalance η = max(n1,n2)/min(n1,n2)
fall wherever a random Kingman tree puts it. This runner instead reuses the
**η-controlled** trees already built in ``cache/pool_sample/`` (see
``eta_pool``/``build_eta_pool``), so we can report the B-method's subsampling
recovery at fixed η bins (η≈1, η≈5, ...) and stay aligned with the
``kingman_threshold_vs_theory`` analysis.

Each pool sample stores the JC *similarity* matrix ``M`` (∈ (0,1], diag 1); the
B-method distance is recovered losslessly as ``D = -log(M)`` — no re-simulation.
Each sample's full ``p × reps`` sweep is disk-cached **per sample** via
``utils.bpart_sweep_cache`` (same granularity + mechanism as
``sweep_cache.compute_or_load_sweep``), so re-running — or adding an ``n`` /
widening an η bin — only computes the genuinely new samples; everything else is
a cache hit. We average the cached raw metrics across the samples in a bin and
report the bin's realized η as mean ± std of the samples' metadata η.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from ..utils.bpart_sweep_cache import aggregate_per_p, compute_or_load_bpart_sweep
from ..utils.eta_pool_cache import (
    param_key, bin_name, list_completed_samples, load_pool_entry,
)
from ..utils.logging import log_info
from .nj_sweep import _SIM_FLOOR

_DEFAULT_P = [1.0, 0.9, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]


def _D_from_M(M: np.ndarray) -> np.ndarray:
    """JC distance from a stored similarity matrix: D = -log(clip(M)), diag 0."""
    R = np.clip(np.asarray(M, dtype=np.float64), _SIM_FLOOR, 1.0)
    D = -np.log(R)
    np.fill_diagonal(D, 0.0)
    return D


def bpart_eta_pool_for_n(
    n: int, run_dir: str, *,
    eta_targets: Sequence[int] = (1, 5),
    reps: int = 1,
    p_values: Sequence[float] = tuple(_DEFAULT_P),
    seq_len: int = 10000, mu: float = 0.1, pop_size: float = 1.0,
    seq_model: str = "JC69", max_samples: int | None = None,
    seed_base: int = 0,
) -> Dict:
    """Run the B-method subsampling sweep on every pooled Kingman tree at the
    given η bins for one ``n``. Writes ``bpart_eta_pool_n{n}.json`` and returns it.
    """
    key = param_key(n, seq_len, mu, "kingman", pop_size, seq_model)
    p_values = list(p_values)

    bins: Dict[str, Dict] = {}
    for target in eta_targets:
        idxs = list_completed_samples(None, key, int(target))
        if max_samples is not None:
            idxs = idxs[:max_samples]
        if not idxs:
            log_info("bpart_eta", f"  n={n} eta~{target}: no pool samples, skipping", force=True)
            continue
        log_info("bpart_eta", f"  n={n} eta~{target}: {len(idxs)} samples × {reps} reps", force=True)

        etas: List[float] = []
        raws: List[List[Dict]] = []
        hits = 0
        for idx in idxs:
            entry = load_pool_entry(None, key, int(target), idx)
            if entry is None:
                continue
            M, _v_pop, _part, _tree, meta = entry

            outcome = compute_or_load_bpart_sweep(
                key, bin_name(int(target)), f"sample_{int(idx):04d}",
                D_loader=lambda M=M: _D_from_M(M),
                p_values=p_values, reps=reps,
                seed_base=seed_base, imputation="mean",
            )
            if outcome is None:
                continue
            result, was_cached = outcome
            etas.append(float(meta["eta"]))
            raws.append(result["per_p"])
            hits += int(was_cached)

        if not raws:
            continue
        log_info("bpart_eta",
                 f"    {hits}/{len(raws)} samples from cache", force=True)
        bins[str(int(target))] = {
            "eta_target": int(target),
            "eta_mean": float(np.mean(etas)), "eta_std": float(np.std(etas)),
            "n_samples": len(etas), "reps": reps,
            "per_p": aggregate_per_p(raws),
        }

    out = {"n": n, "seq_len": seq_len, "mutation_rate": mu, "pop_size": pop_size,
           "tree_model": "kingman", "bins": bins}
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    out_path = Path(run_dir) / f"bpart_eta_pool_n{n}.json"
    out_path.write_text(json.dumps(out, indent=2))
    log_info("bpart_eta", f"Wrote {out_path}", force=True)
    return out
