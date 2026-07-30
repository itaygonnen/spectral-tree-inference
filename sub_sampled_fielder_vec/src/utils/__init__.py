"""Utilities for the sub-sampled-Fiedler experiments.

Deliberately does NOT re-export anything. ``src/utils`` holds ~30 loosely related
modules -- caching shims, metrics, plotting, per-method IO, logging -- and several are
expensive to import (``tree_plots`` pulls toytree/toyplot/PIL, ``plotting`` is 1,150
lines of matplotlib). A flat re-export here would make every ``from src.utils.x import
y`` pay for all of them. Import the specific module instead:

    from src.utils.threshold_utils import find_discrete_threshold
    from src.utils.eta_pool_cache import param_key, load_pool_entry

This file exists because ``src/utils`` was the only subpackage without an
``__init__.py``, which left it an implicit namespace package -- it worked, but it made
``src.utils`` behave subtly differently from its siblings.

Caching lives here: ``persistent_cache`` (experiment_data), ``sweep_cache``
(bootstrap_sweep), ``eta_pool_cache`` (pool_sample) and ``bpart_sweep_cache``
(bpart_sweep) are thin shims over the scopes declared in ``src/cache_io.py``. Each
IGNORES its ``cache_root`` argument -- the real root is always ``cache_io.CACHE_ROOT``.
See ``docs/CACHE_AND_RESULTS.md``.
"""
