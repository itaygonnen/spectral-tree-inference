"""Shared utilities for theoretical-interpretation notebooks.

Topology-specific code lives in the ``balanced_binary`` and ``block_model``
modules (and future siblings). Everything else is generic over the similarity
matrix.

The figure modules ``sweep_plots``, ``sweep_plots_two_panel`` and the loader
``generated_data`` are NOT re-exported here on purpose: notebooks import them by
module path so it stays visible at the call site which figure family a plot belongs
to (``sweep_plots`` -> appendix Figs 8-9, ``sweep_plots_two_panel`` -> main Figs
2-3). Those two modules estimate the constant ``C`` differently and their values are
not comparable, so collapsing them into one flat namespace would invite mixing them.

``distance_similarity`` and ``distance_features`` used to live here and were
documented as the shared engine behind the two distance_vs_similarity notebooks. They
had zero importers -- the notebooks inline their own logic -- and were removed;
recover them from git history if that extraction is ever revived.
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
from .perturbation import (
    rank2_subspace, laplacian_perturbation_norm,
    subspace_sin_theta, davis_kahan_bound, two_mode_decomposition,
)
from .recovery import compute_recovery, compute_ari, compute_nmi, find_threshold_p_star
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
    "rank2_subspace",
    "laplacian_perturbation_norm",
    "subspace_sin_theta",
    "davis_kahan_bound",
    "two_mode_decomposition",
    "compute_fiedler_of_S",
    "subsample_S",
    "uniform_mask",
    "ipw_from_mask",
    "nnm_from_mask",
    "compute_recovery",
    "compute_ari",
    "compute_nmi",
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
