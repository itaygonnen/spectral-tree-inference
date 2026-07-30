"""Compatibility shim — the generated-tree data source moved into ``src/``.

``simulation_distance_vs_similarity.ipynb`` imports ``make_generated`` and
``build_ids`` from here, so the names stay put. Everything now lives where the
rest of the pipeline can reach it without ``src/`` importing ``analysis/``:

* id grammar (``build_ids``, ``parse_id``, ``parse_spec``, ``make_tree_id``,
  ``GenSpec``, ``MODEL_OF_TAG``) -> ``src/utils/tree_ids.py``
* generation (``make_generated``, ``TREE_BUILDERS``, ``GENERATORS``)
  -> ``src/models/generated_trees.py``

New code should import from those two modules directly.
"""
from __future__ import annotations

from src.models.generated_trees import (  # noqa: F401
    GENERATORS, TREE_BUILDERS, make_generated,
)
from src.utils.tree_ids import (  # noqa: F401
    MODEL_OF_TAG, GenSpec, build_ids, make_tree_id, parse_id, parse_spec,
)

__all__ = ["GENERATORS", "TREE_BUILDERS", "make_generated", "MODEL_OF_TAG",
           "GenSpec", "build_ids", "make_tree_id", "parse_id", "parse_spec"]
