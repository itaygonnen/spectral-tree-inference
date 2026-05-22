"""Griffing-distance-partition sweep over ``p`` for a fixed (n_taxa, seq_len).

Distance-matrix analogue of figure_3's Fiedler-recovery sweep. Mirrors the
shape of ``nj_sweep.nj_sweep_for_params`` but the inner metric is
*sign-agnostic recovery of the leading eigenvector of* ``J D̂ J`` against
the reference eigenvector of ``J D J`` (the full, unsampled matrix).

  1. ``||D - D̂||_2``  (spec-norm error; same dimensional artifact as NJ)
  2. ``recovery``     (sign-agnostic agreement, [0.5, 1.0])

Sampling defaults to **zero-fill IPW** (faithful to the paper's
no-completion baseline) because the eigenvector of B is scale-invariant
under the mean-imputation bias, so we don't need the band-aid that NJ
required. ``imputation="mean"`` is exposed as a parameter for direct
comparison.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from ..config import StructuredConfig
from ..core.similarity_builder import SimilarityMatrixBuilder
from ..utils.logging import log_info
from ..utils.metrics import estimate_operator_norm_diff, compute_sign_agreement
from ..utils.griffing import griffing_leading_eigvec
from ..utils.summaries import save_json
from .nj_sweep import _build_truth, _impute_mean


def griffing_sweep_for_params(cfg: StructuredConfig, n_taxa: int, seq_len: int,
                              run_dir: str, *, imputation: str = "zero",
                              save_eigvecs: bool = False) -> Dict:
    """Run the Griffing-partition sweep at fixed ``(n_taxa, seq_len)``.

    Writes ``griffing_meta.json`` (and optionally ``griffing_arrays.npz``)
    into ``run_dir`` and returns the in-memory dict.

    Parameters
    ----------
    imputation : {"zero", "mean"}
        Default ``"zero"`` matches the paper's no-completion baseline.
        ``"mean"`` reuses ``nj_sweep._impute_mean`` for direct comparison
        with figures 8/9.
    save_eigvecs : bool
        If True, also writes the reference and per-(p, rep) eigenvectors to
        ``griffing_arrays.npz`` (one ``(m,)`` vector for v_ref, and
        ``(n_p, n_reps, m)`` for v_hat). Off by default — at n=2048 that's
        ~80MB per sweep, useful for offline analysis but not for the headline
        recovery plot.
    """
    if imputation not in ("zero", "mean"):
        raise ValueError(f"imputation must be 'zero' or 'mean', got {imputation!r}")
    os.makedirs(run_dir, exist_ok=True)
    log_info("griffing", f"Building ground truth (n={n_taxa}, L={seq_len})...", force=True)
    tree, observations, D, taxa_metadata = _build_truth(cfg, n_taxa, seq_len)

    # Reference eigenvector — leading by |eigenvalue| of J D J
    log_info("griffing", "Computing reference v_ref = leading eigvec of J·D·J...", force=True)
    v_ref = griffing_leading_eigvec(D)

    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="distance")
    sampler = builder.sampler  # UniformSampler with self_value=0.0

    p_values: List[float] = list(cfg.experiment.p_values)
    bootstrap_reps: int = int(cfg.experiment.bootstrap_reps)
    seed_base: int = int(cfg.experiment.seed)

    per_p: List[Dict] = []
    v_hat_stack: List[np.ndarray] = [] if save_eigvecs else []
    t0 = time.time()
    for p_idx, p in enumerate(p_values):
        log_info("griffing", f"  p = {p:.4g} ({bootstrap_reps} reps)", force=True)
        spec_norms: List[float] = []
        recoveries: List[float] = []
        v_hat_reps: List[np.ndarray] = [] if save_eigvecs else []

        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            D_hat = sampler.sample(D, p, seed=seed)
            if imputation == "mean":
                D_hat = _impute_mean(D_hat)

            spec_norm = float(estimate_operator_norm_diff(D_hat, D, n_iter=20))
            spec_norms.append(spec_norm)

            v_hat = griffing_leading_eigvec(D_hat)
            # Sign-agnostic agreement in [0.5, 1.0] — matches utils/recovery.py
            s = compute_sign_agreement(v_ref, v_hat) / 100.0
            recovery = float(max(s, 1.0 - s))
            recoveries.append(recovery)

            if save_eigvecs:
                v_hat_reps.append(v_hat)

        per_p.append({
            "p": float(p),
            "spec_norm_mean": float(np.mean(spec_norms)),
            "spec_norm_std": float(np.std(spec_norms)),
            "spec_norm_per_rep": [float(x) for x in spec_norms],
            "recovery_mean": float(np.mean(recoveries)),
            "recovery_std": float(np.std(recoveries)),
            "recovery_per_rep": [float(x) for x in recoveries],
        })
        if save_eigvecs:
            v_hat_stack.append(np.stack(v_hat_reps, axis=0))

    elapsed = time.time() - t0
    meta = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": [float(p) for p in p_values],
        "bootstrap_reps": bootstrap_reps,
        "tree_model": cfg.tree.model,
        "seq_model": cfg.sequence.model,
        "mutation_rate": float(cfg.sequence.params.get("mutation_rate", 1.0)),
        "tree_params": dict(cfg.tree.params),
        "sampling_method": f"uniform_imputation={imputation}",
        "method": "griffing_distance_partition",
        "imputation": imputation,
        "per_p": per_p,
        "elapsed_seconds": elapsed,
    }
    run_path = Path(run_dir)
    meta_path = run_path / "griffing_meta.json"
    save_json(meta, str(meta_path))
    if save_eigvecs:
        npz_path = run_path / "griffing_arrays.npz"
        np.savez_compressed(npz_path,
                            v_ref=v_ref,
                            v_hat_stack=np.stack(v_hat_stack, axis=0))
        log_info("griffing", f"Wrote {meta_path} + {npz_path.name} ({elapsed:.1f}s)",
                 force=True)
    else:
        log_info("griffing", f"Wrote {meta_path} ({elapsed:.1f}s)", force=True)
    return meta
