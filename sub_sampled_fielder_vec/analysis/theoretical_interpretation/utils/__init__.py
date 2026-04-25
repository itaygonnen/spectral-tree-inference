"""Shared utilities for theoretical-interpretation notebooks.

Topology-specific code lives in the ``balanced_binary`` module (and future
siblings). Everything else is generic over the similarity matrix.
"""
from .balanced_binary import build_balanced_binary_S, balanced_binary_population_fiedler
from .spectral import compute_fiedler_of_S, subsample_S
from .recovery import compute_recovery, find_threshold_p_star
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
from .sweep import run_sweep
from .plotting import compute_C_constant, plot_recovery_figure

__all__ = [
    "build_balanced_binary_S",
    "balanced_binary_population_fiedler",
    "compute_fiedler_of_S",
    "subsample_S",
    "compute_recovery",
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
    "compute_C_constant",
    "plot_recovery_figure",
]
