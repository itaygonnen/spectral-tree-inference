"""Shared utilities for theoretical-interpretation notebooks.

Topology-specific code lives in the ``balanced_binary`` and ``block_model``
modules (and future siblings). Everything else is generic over the similarity
matrix.
"""
from .balanced_binary import build_balanced_binary_S, balanced_binary_population_fiedler
from .block_model import build_flat_cbm_S, build_decay_cbm_S, flat_cbm_population_fiedler
from .tree_features import (
    imbalance_eta, n_min, structural_margin_rho, estimate_features_from_M,
    hbm_d_max, hbm_s_in_min,
)
from .linalg_features import (
    coherence_mu, cross_clan_variance,
    compute_top_eigenpairs, spectral_gap, compute_lemma04_row,
)
from .spectral import (
    compute_fiedler_of_S, subsample_S,
    uniform_mask, ipw_from_mask, nnm_from_mask,
)
from .recovery import compute_recovery, compute_ari, find_threshold_p_star
from .cache import (
    make_full_key,
    make_subsample_key,
    load_full,
    save_full,
    get_or_compute_full,
    load_subsample,
    save_subsample,
    get_or_compute_subsample,
    clear_theoretical_cache,
)
from .sweep import run_sweep, load_or_extend_metrics
from .plotting import compute_C_constant, plot_recovery_figure

__all__ = [
    "build_balanced_binary_S",
    "balanced_binary_population_fiedler",
    "build_flat_cbm_S",
    "build_decay_cbm_S",
    "flat_cbm_population_fiedler",
    "imbalance_eta",
    "n_min",
    "structural_margin_rho",
    "estimate_features_from_M",
    "hbm_d_max",
    "hbm_s_in_min",
    "coherence_mu",
    "cross_clan_variance",
    "compute_top_eigenpairs",
    "spectral_gap",
    "compute_lemma04_row",
    "compute_fiedler_of_S",
    "subsample_S",
    "uniform_mask",
    "ipw_from_mask",
    "nnm_from_mask",
    "compute_recovery",
    "compute_ari",
    "find_threshold_p_star",
    "make_full_key",
    "make_subsample_key",
    "load_full",
    "save_full",
    "get_or_compute_full",
    "load_subsample",
    "save_subsample",
    "get_or_compute_subsample",
    "clear_theoretical_cache",
    "run_sweep",
    "load_or_extend_metrics",
    "compute_C_constant",
    "plot_recovery_figure",
]
