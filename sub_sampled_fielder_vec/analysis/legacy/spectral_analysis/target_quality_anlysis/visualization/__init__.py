"""Visualization modules for plotting diagnostic results."""

from .scree_plots import plot_scree, generate_all_scree_plots, plot_comparison_scree
from .coherence_plots import plot_coherence_comparison
from .heatmaps import plot_metric_heatmaps
# tree_plots.py moved to src/utils/tree_plots.py: src/ imports it (pipeline A
# diagnostics) and src must never import from analysis/.
from src.utils.tree_plots import plot_tree_with_partition, plot_combined_tree_and_fiedler

__all__ = [
    'plot_scree',
    'generate_all_scree_plots',
    'plot_comparison_scree',
    'plot_coherence_comparison',
    'plot_metric_heatmaps',
    'plot_tree_with_partition',
    'plot_combined_tree_and_fiedler'
]

