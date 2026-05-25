"""Reusable bootstrap p-sweep over a precomputed (M, fiedler_ref) pair.

Narrow extraction of the inner loop in ``bootstrap_sweep.py::sweep_for_params``
(lines ~498–722) for the **uniform sampling** path, focused on
``partition_agreement_M``. Drops the metric_composer, leveraged-sampler
diagnostics, persistent caching, guardrails, and partition-validity checks
that the full runner needs but a notebook does not.

Use ``sweep_for_params`` for full experiments; use this helper when you
already have ``(M, v_pop)`` (e.g. from the eta-pool cache) and want a
quick agreement curve.
"""
from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional

import numpy as np

from ..core.utils import compute_laplacian
from ..utils.logging import log_warning
from ..utils.metrics import (
    _normalize_vector,
    compute_fiedler_dot_product,
    compute_reference_partition_and_quality,
    compute_sign_agreement,
)
from ..utils.random_entries import _subsample_matrix_entries, compute_fiedler_from_laplacian


def _bipartition_agreement(partition_ref: np.ndarray, partition_avg: np.ndarray) -> float:
    matches_direct = int(np.sum(partition_ref == partition_avg))
    matches_flipped = int(np.sum(partition_ref != partition_avg))
    return 100.0 * max(matches_direct, matches_flipped) / len(partition_ref)


def align_fiedler_by_dot_product(
    fiedler_vector: np.ndarray, reference_vector: np.ndarray
) -> np.ndarray:
    """Normalize and sign-align a Fiedler vector against a reference."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            v_normalized = _normalize_vector(fiedler_vector)
            u_normalized = _normalize_vector(reference_vector)
        except ValueError as e:
            log_warning('align', f"Normalization failed: {e}")
            return fiedler_vector
        dot_product = np.dot(v_normalized, u_normalized)
        if not np.isfinite(dot_product):
            log_warning('align', f"Non-finite dot product: {dot_product}")
            return fiedler_vector
        return -v_normalized if dot_product < 0 else v_normalized


def bootstrap_p_sweep_simple(
    M: np.ndarray,
    fiedler_ref: np.ndarray,
    p_values: List[float],
    bootstrap_reps: int = 50,
    seed: int = 0,
    num_gaps: int = 10,
    min_split: int = 1,
    partition_ref: Optional[np.ndarray] = None,
    early_stop_consecutive_100: int = 0,
    partition_method: str = "sigma2",
) -> Dict[str, Any]:
    """Bootstrap p-sweep with precomputed (M, fiedler_ref). Uniform sampling.

    Parameters
    ----------
    M : (n, n) similarity matrix.
    fiedler_ref : (n,) reference Fiedler vector of L(M).
    p_values : sampling rates to sweep.
    bootstrap_reps : bootstrap iterations per p (skipped when p ≈ 1.0).
    seed : base seed; rep i uses seed+i.
    num_gaps, min_split : passed to partition_taxa (sigma2 mode only).
    partition_ref : optional precomputed reference partition. If None,
        derived from fiedler_ref according to ``partition_method``.
    early_stop_consecutive_100 : if > 0, stop the sweep after this many
        consecutive p's hit partition_agreement_M == 100.0, and fill the
        remaining p's with (100.0, 100.0, 1.0). Assumes p_values are
        ordered such that agreement is monotone-ish at the high-p tail.
    partition_method : how to derive partitions from Fiedler vectors.
        - ``"sigma2"`` (default): use ``partition_taxa`` (gap search by σ₂,
          honours ``num_gaps``/``min_split``). Can return degenerate splits
          when an outlier taxon dominates.
        - ``"sign"``: split by sign of the Fiedler vector (``v > 0`` vs
          ``v <= 0``). No σ₂ search; matches the eta-binning convention.

    Returns
    -------
    dict with keys: ``p_values``, ``partition_agreement_M``,
    ``partition_ari_M``, ``partition_nmi_M``, ``sign_agreement``,
    ``dot_product``, ``partition_split_ref``, ``reference_partition_quality``.
    ``partition_ari_M`` is the Adjusted Rand Index and ``partition_nmi_M``
    is the Normalized Mutual Information between the reference bipartition
    and the bootstrap-averaged bipartition (both label-permutation invariant).
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    from spectraltree.spectral_tree_reconstruction import partition_taxa

    if partition_method not in ("sigma2", "sign"):
        raise ValueError(f"partition_method must be 'sigma2' or 'sign', got {partition_method!r}")

    if partition_ref is None:
        if partition_method == "sigma2":
            partition_ref, ref_quality, ref_split = compute_reference_partition_and_quality(
                fiedler_ref, M, num_gaps=num_gaps, min_split=min_split
            )
        else:
            partition_ref = fiedler_ref > 0
            ref_quality = float('nan')
            n_true = int(np.sum(partition_ref))
            n_false = len(partition_ref) - n_true
            ref_split = (min(n_true, n_false), max(n_true, n_false))
    else:
        ref_quality = float('nan')
        n_true = int(np.sum(partition_ref))
        n_false = len(partition_ref) - n_true
        ref_split = (min(n_true, n_false), max(n_true, n_false))

    partition_ref_int = partition_ref.astype(int)

    partition_agreement_M: List[float] = []
    partition_ari_M: List[float] = []
    partition_nmi_M: List[float] = []
    sign_agreement: List[float] = []
    dot_product: List[float] = []
    consecutive_100 = 0

    for idx, p in enumerate(p_values):
        if early_stop_consecutive_100 > 0 and consecutive_100 >= early_stop_consecutive_100:
            remaining = len(p_values) - idx
            partition_agreement_M.extend([100.0] * remaining)
            partition_ari_M.extend([1.0] * remaining)
            partition_nmi_M.extend([1.0] * remaining)
            sign_agreement.extend([100.0] * remaining)
            dot_product.extend([1.0] * remaining)
            break

        if p >= 0.9999:
            partition_agreement_M.append(100.0)
            partition_ari_M.append(1.0)
            partition_nmi_M.append(1.0)
            sign_agreement.append(100.0)
            dot_product.append(1.0)
            consecutive_100 += 1
            continue

        aligned: List[np.ndarray] = []
        for i in range(bootstrap_reps):
            S = _subsample_matrix_entries(M, p, seed=seed + i, builder=None)
            try:
                L_S = compute_laplacian(S)
                f_est = compute_fiedler_from_laplacian(L_S)
            except Exception as e:
                log_warning('p_sweep_inner', f"Fiedler failed at p={p:.4g}, rep={i}: {e}")
                continue
            aligned.append(align_fiedler_by_dot_product(f_est, fiedler_ref))

        if not aligned:
            partition_agreement_M.append(float('nan'))
            partition_ari_M.append(float('nan'))
            partition_nmi_M.append(float('nan'))
            sign_agreement.append(float('nan'))
            dot_product.append(float('nan'))
            continue

        try:
            v_avg = _normalize_vector(np.mean(aligned, axis=0))
        except ValueError:
            v_avg = np.zeros_like(fiedler_ref)

        try:
            if partition_method == "sigma2":
                partition_avg = partition_taxa(v_avg, M, num_gaps, min_split)
            else:
                partition_avg = v_avg > 0
            agr_M = _bipartition_agreement(partition_ref, partition_avg)
            ari_M = float(adjusted_rand_score(partition_ref_int, partition_avg.astype(int)))
            nmi_M = float(normalized_mutual_info_score(partition_ref_int, partition_avg.astype(int)))
        except Exception as e:
            log_warning('p_sweep_inner', f"partition_agreement failed at p={p:.4g}: {e}")
            agr_M = float('nan')
            ari_M = float('nan')
            nmi_M = float('nan')
        partition_agreement_M.append(float(agr_M))
        partition_ari_M.append(ari_M)
        partition_nmi_M.append(nmi_M)
        sign_agreement.append(float(compute_sign_agreement(fiedler_ref, v_avg)))
        try:
            dot_product.append(float(compute_fiedler_dot_product(fiedler_ref, v_avg)))
        except ValueError as e:
            log_warning('p_sweep_inner', f"dot_product failed at p={p:.4g}: {e}")
            dot_product.append(float('nan'))
        consecutive_100 = consecutive_100 + 1 if agr_M == 100.0 else 0

    return {
        "p_values": list(p_values),
        "partition_agreement_M": partition_agreement_M,
        "partition_ari_M": partition_ari_M,
        "partition_nmi_M": partition_nmi_M,
        "sign_agreement": sign_agreement,
        "dot_product": dot_product,
        "partition_split_ref": ref_split,
        "reference_partition_quality": ref_quality,
    }
