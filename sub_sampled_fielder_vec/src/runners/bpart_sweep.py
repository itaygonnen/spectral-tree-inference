"""Distance-based spectral clan-partitioning sweep over ``p`` (Fig. 9 variant).

Replaces the classical-NJ path (``nj_sweep``) with the distance-based
partitioner of Spectral-Top-Down-Recovery, Appendix D (p. 2329, from [20]):
the terminal nodes split into two clans by the **sign pattern of the leading
eigenvector** of the doubly-centred distance matrix

    B = (I - 11^T/m) D (I - 11^T/m).

We measure how well that bipartition survives uniform sub-sampling of ``D``.
All metrics compare the sub-sampled partition (sign of the leading eigenvector
of ``B̂ = (I-11^T/m) D̂ (I-11^T/m)``) against the full-data reference partition
(sign of the leading eigenvector of ``B``):

  1. ``agreement`` — % of taxa with matching clan label. Orientation-invariant.
  2. ``dot``       — |⟨v̂₁, v₁⟩| of the sign-aligned unit-norm eigenvectors.
  3. ``ari``       — Adjusted Rand Index of the two sign bipartitions.
  4. ``nmi``       — Normalized Mutual Information of the two sign bipartitions.

The inner sweep + its disk cache live in ``utils.bpart_sweep_cache`` (the
canonical ``CacheScope`` mechanism, shared by all three bpart runners); this
module only builds the distance matrix and renders ``bpart_meta.json``. The
bipartition-comparison helpers (``_align_sign`` etc.) are re-exported from there
for backward compatibility.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List

from ..cache_io import make_key
from ..config import StructuredConfig
from ..utils.bpart_sweep_cache import (  # re-exported helpers + cached sweep
    _align_sign, _partition_agreement, _binary_ari, _binary_nmi,
    aggregate_per_p, compute_or_load_bpart_sweep,
)
from ..utils.logging import log_info
from ..utils.summaries import save_json
from .nj_sweep import cached_distance_matrix

__all__ = [
    "_align_sign", "_partition_agreement", "_binary_ari", "_binary_nmi",
    "bpart_sweep_for_params",
]


def bpart_sweep_for_params(cfg: StructuredConfig, n_taxa: int, seq_len: int,
                           run_dir: str, *, imputation: str = "mean") -> Dict:
    """Run the B-matrix clan-partition sweep at fixed (n_taxa, seq_len).

    Writes ``bpart_meta.json`` into ``run_dir`` and returns the in-memory dict.
    The sweep itself is disk-cached (``utils.bpart_sweep_cache``), so re-running
    with the same ``(tree, n, L, mu, seq, tree_params)`` and the same p-grid /
    reps / seed / imputation is an instant cache hit.

    Parameters
    ----------
    imputation : {"mean", "zero"}
        How to handle unsampled entries before double-centering (see
        ``nj_sweep`` for the rationale). Default ``"mean"``.
    """
    if imputation not in ("mean", "zero"):
        raise ValueError(f"imputation must be 'mean' or 'zero', got {imputation!r}")
    os.makedirs(run_dir, exist_ok=True)
    log_info("bpart", f"Building ground truth (n={n_taxa}, L={seq_len})...", force=True)

    p_values: List[float] = list(cfg.experiment.p_values)
    bootstrap_reps: int = int(cfg.experiment.bootstrap_reps)
    seed_base: int = int(cfg.experiment.seed)

    # Identity = the same (tree, n, L, mu, seq, *tree_params) tuple the distance
    # cache uses; the p-grid / reps / seed / imputation are folded into the
    # sweep key inside compute_or_load_bpart_sweep.
    identity = make_key(
        f"bpart_{cfg.tree.model}",
        n=n_taxa, L=seq_len,
        mu=float(cfg.sequence.params.get("mutation_rate", 1.0)),
        seq=cfg.sequence.model,
        **{k: cfg.tree.params[k] for k in sorted(cfg.tree.params)},
    )

    t0 = time.time()
    result, was_cached = compute_or_load_bpart_sweep(
        identity,
        D_loader=lambda: cached_distance_matrix(cfg, n_taxa, seq_len),
        p_values=p_values, reps=bootstrap_reps,
        seed_base=seed_base, imputation=imputation,
    )
    per_p = aggregate_per_p([result["per_p"]], include_per_rep=True)
    elapsed = time.time() - t0

    meta = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": [float(p) for p in p_values],
        "bootstrap_reps": bootstrap_reps,
        "tree_model": cfg.tree.model,
        "seq_model": cfg.sequence.model,
        "mutation_rate": float(cfg.sequence.params.get("mutation_rate", 1.0)),
        "n1": result["n1"], "n2": result["n2"], "eta": result["eta"],
        "sampling_method": f"uniform_imputation={imputation}",
        "method": "distance_double_centering",
        "distance_kind": "jc_minus_log",
        "imputation": imputation,
        "per_p": per_p,
        "elapsed_seconds": elapsed,
        "was_cached": was_cached,
    }
    meta_path = Path(run_dir) / "bpart_meta.json"
    save_json(meta, str(meta_path))
    log_info("bpart", f"Wrote {meta_path} ({elapsed:.1f}s, "
                      f"cache {'HIT' if was_cached else 'MISS'})", force=True)
    return meta
