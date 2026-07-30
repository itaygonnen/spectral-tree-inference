"""Disk cache for the B-method (distance double-centering) subsampling sweep.

The bpart family (``bpart_sweep``, ``bpart_eta_pool``, ``bpart_synthetic``)
recovers a clan bipartition from the sign of the leading eigenvector of
``B = (I-11ᵀ/m) D (I-11ᵀ/m)`` and measures how it survives uniform
sub-sampling of ``D``. This module is the **single** place that runs that
inner sweep and the **single** place that caches it — the direct analogue of
:mod:`src.utils.sweep_cache` (which does the same for the Fiedler-on-S
``bootstrap_p_sweep_simple`` path).

Mechanism is the canonical one used everywhere in the repo: a
:class:`~src.cache_io.CacheScope` (here ``bpart_sweep``), a ``sha1(config)``
sweep key, a ``{config, result}`` payload, and a ``(result, was_cached)``
return. One cached entry == one full ``p × reps`` sweep on one distance matrix,
addressed by caller-supplied ``identity_parts`` (which matrix) plus the sweep
key (which p-grid / reps / seed / imputation).
"""
from __future__ import annotations

import hashlib
import json
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..cache_io import bpart_sweep as _scope
from ..core.similarity_builder import SimilarityMatrixBuilder
from ..utils.griffing import DEFAULT_SOLVER, griffing_leading_eigvec
from ..runners.nj_sweep import _impute_mean

SCHEMA_VERSION = "v1"
_METRICS = ("agreement", "dot", "ari", "nmi")

# UniformSampler with self_value=0.0 (distance convention). Stateless → reuse.
_SAMPLER = SimilarityMatrixBuilder(method="uniform", matrix_kind="distance").sampler


# --- bipartition comparison helpers (lower layer: the runners import these) ---
def _align_sign(v: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Flip ``v`` so its sign pattern best matches ``ref`` (sign-vector dot)."""
    s = np.dot(np.sign(v), np.sign(ref))
    return -v if s < 0 else v


def _partition_agreement(v_ref: np.ndarray, v_hat: np.ndarray) -> float:
    """Orientation-invariant % of matching clan labels (sign >= 0)."""
    a = v_ref >= 0
    b = v_hat >= 0
    matches = int(np.sum(a == b))
    matches = max(matches, len(a) - matches)
    return 100.0 * matches / len(a)


def _binary_ari(v_ref: np.ndarray, v_hat: np.ndarray) -> float:
    """Adjusted Rand Index of the two sign-based bipartitions (label-permutation
    invariant, so robust to a global sign flip). 1 = perfect, ~0 = random."""
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(np.asarray(v_ref) >= 0, np.asarray(v_hat) >= 0))


def _binary_nmi(v_ref: np.ndarray, v_hat: np.ndarray) -> float:
    """Normalized Mutual Information of the two sign-based bipartitions
    (label-permutation invariant). 1 = perfect, 0 = independent."""
    from sklearn.metrics import normalized_mutual_info_score
    return float(normalized_mutual_info_score(np.asarray(v_ref) >= 0, np.asarray(v_hat) >= 0))


def _build_config(
    p_values: Sequence[float], reps: int, seed_base: int, imputation: str,
    method: str = "griffing_double_centering",
    eigsolver: str = DEFAULT_SOLVER,
) -> Dict[str, Any]:
    """Config dict that IS the cache key (sha1'd by :func:`compute_sweep_key`).

    ``eigsolver`` is emitted **only when it differs from the historical default**.
    Adding it unconditionally would change every hash and orphan the sweeps
    already on disk under ``cache/bpart_sweep/``; omitting it on ``"eigh_full"``
    keeps those keys byte-identical while ``"lm_k1"`` runs get their own
    namespace, so the two implementations can never be silently averaged
    together.
    """
    config: Dict[str, Any] = {
        "p_values": [round(float(p), 9) for p in p_values],
        "reps": int(reps),
        "seed_base": int(seed_base),
        "imputation": str(imputation),
        "method": method,
        "schema_version": SCHEMA_VERSION,
    }
    if eigsolver != DEFAULT_SOLVER:
        config["eigsolver"] = str(eigsolver)
    return config


def compute_sweep_key(
    p_values: Sequence[float], reps: int, seed_base: int, imputation: str,
    eigsolver: str = DEFAULT_SOLVER,
) -> Tuple[str, Dict[str, Any]]:
    config = _build_config(p_values, reps, seed_base, imputation,
                           eigsolver=eigsolver)
    blob = json.dumps(config, sort_keys=True).encode("utf-8")
    return hashlib.sha1(blob).hexdigest()[:12], config


def bpart_sweep_raw(
    D: np.ndarray, p_values: Sequence[float], reps: int,
    seed_base: int = 0, imputation: str = "mean",
    eigsolver: str = DEFAULT_SOLVER,
) -> Dict[str, Any]:
    """Run one B-method subsampling sweep on a single distance matrix ``D``.

    Returns the reference split ``(n1, n2, eta)`` and, for each ``p``, the
    **raw per-rep** value of each metric (no aggregation — callers average as
    they wish, e.g. across pool samples)::

        {"p_values":[...], "n1":int, "n2":int, "eta":float,
         "per_p":[{"p":p, "agreement":[…reps], "dot":[…], "ari":[…], "nmi":[…]}, …]}

    ``eigsolver`` is used for the reference vector **and** every sub-sampled one,
    so a sweep never compares across two implementations.
    """
    if imputation not in ("mean", "zero"):
        raise ValueError(f"imputation must be 'mean' or 'zero', got {imputation!r}")
    v_ref = griffing_leading_eigvec(D, eigsolver)
    n1 = int(np.sum(v_ref >= 0))
    n2 = int(v_ref.size - n1)
    eta = float(max(n1, n2)) / float(max(min(n1, n2), 1))

    per_p: List[Dict[str, Any]] = []
    for p_idx, p in enumerate(p_values):
        acc: Dict[str, List[float]] = {m: [] for m in _METRICS}
        for rep in range(reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            D_hat = _SAMPLER.sample(D, float(p), seed=seed)
            if imputation == "mean":
                D_hat = _impute_mean(D_hat)
            v_hat = _align_sign(griffing_leading_eigvec(D_hat, eigsolver), v_ref)
            acc["agreement"].append(_partition_agreement(v_ref, v_hat))
            acc["dot"].append(float(np.abs(np.dot(v_hat, v_ref))))
            acc["ari"].append(_binary_ari(v_ref, v_hat))
            acc["nmi"].append(_binary_nmi(v_ref, v_hat))
        per_p.append({"p": float(p), **acc})

    return {
        "p_values": [float(p) for p in p_values],
        "n1": n1, "n2": n2, "eta": eta, "per_p": per_p,
    }


def compute_or_load_bpart_sweep(
    *identity_parts: str,
    D_loader: Callable[[], Optional[np.ndarray]],
    p_values: Sequence[float],
    reps: int,
    seed_base: int = 0,
    imputation: str = "mean",
    eigsolver: str = DEFAULT_SOLVER,
    use_cache: bool = True,
) -> Optional[Tuple[Dict[str, Any], bool]]:
    """Return ``(result, was_cached)`` for one matrix's sweep, or ``None``.

    ``identity_parts`` name *which matrix* (e.g. a synthetic param key, or
    ``param_key, bin_name(eta), "sample_0001"`` for a pool sample); the sweep
    key (p-grid / reps / seed / imputation / non-default eigsolver) is appended
    automatically. The distance matrix is built lazily via ``D_loader`` only on a
    cache miss; ``None`` is returned when ``D_loader`` yields ``None`` (e.g. a
    missing pool sample), so callers can skip cleanly.
    """
    sweep_key, config = compute_sweep_key(p_values, reps, seed_base, imputation,
                                          eigsolver=eigsolver)
    parts = tuple(str(x) for x in identity_parts) + (sweep_key,)

    if use_cache and _scope.is_complete(*parts):
        data = _scope.try_load(*parts)
        if data is not None and "result" in data:
            return data["result"], True

    D = D_loader()
    if D is None:
        return None
    raw = bpart_sweep_raw(D, p_values, reps, seed_base, imputation, eigsolver)
    _scope.save(*parts, result=raw, config=config)
    return raw, False


def aggregate_per_p(
    raws: Sequence[List[Dict[str, Any]]], include_per_rep: bool = False,
) -> List[Dict[str, Any]]:
    """Aggregate one-or-more raw ``per_p`` lists into mean/std per metric.

    ``raws`` is a list of ``per_p`` arrays (each from :func:`bpart_sweep_raw`),
    all sharing the same ``p`` ordering. Values are pooled across every raw
    (across pool samples *and* reps) before computing mean/std — so a single
    matrix passes ``[raw["per_p"]]`` and an η-bin passes one entry per sample.
    ``include_per_rep`` additionally emits the pooled raw list as
    ``{metric}_per_rep`` (used by ``bpart_sweep``'s ``bpart_meta.json``).
    """
    base = raws[0]
    out: List[Dict[str, Any]] = []
    for i, e0 in enumerate(base):
        rec: Dict[str, Any] = {"p": float(e0["p"])}
        for m in _METRICS:
            vals: List[float] = []
            for r in raws:
                vals.extend(r[i][m])
            rec[f"{m}_mean"] = float(np.mean(vals))
            rec[f"{m}_std"] = float(np.std(vals))
            if include_per_rep:
                rec[f"{m}_per_rep"] = [float(x) for x in vals]
        out.append(rec)
    return out
