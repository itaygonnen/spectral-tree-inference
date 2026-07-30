"""Generic operator screen for the distance-vs-similarity notebooks.

Collapses the three near-identical screen cells (Fiedler-on-S, Griffing-on-D,
Fiedler-on-L_sym) into one helper: loop the trees, build a partition with the
given operator+threshold, test whether it is a real single-edge clan, and record
the partition imbalance eta. Lightweight cache (just the per-tree rows).
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, List, Optional, Tuple

import numpy as np

from .partition_validity import check_partition_valid_in_tree
# Re-export so notebooks can `from src.utils.screening import run_screen, _kmeans_bipartition`.
from ..runners.p_sweep_inner import _kmeans_bipartition  # noqa: F401


def _eta(partition: np.ndarray) -> float:
    n1 = int(np.sum(partition))
    n2 = len(partition) - n1
    return max(n1, n2) / max(min(n1, n2), 1)


def _tree_leaf_index(tree, labels) -> np.ndarray:
    idx = {l: i for i, l in enumerate(labels)}
    return np.array([idx[str(l.taxon.label)] for l in tree.leaf_node_iter()])


def run_screen(
    screen_ids: List[str],
    load_all: Callable[[str], Optional[Tuple]],
    partition_fn: Callable[[np.ndarray, np.ndarray], np.ndarray],
    cache_path,
    *,
    label: str = "screen",
    progress: int = 50,
) -> Tuple[List[str], np.ndarray, List[dict]]:
    """Screen ``screen_ids`` with one operator+threshold.

    Parameters
    ----------
    load_all(tid) -> (S, labels, tree, D)  with ``D`` aligned to ``labels`` order
        (or ``None`` to skip a missing tree).
    partition_fn(S, D) -> boolean partition in ``labels`` order.

    Returns
    -------
    (valid_ids, etas_valid, rows) where ``rows`` is a list of
    ``{'tree': tid, 'eta': float, 'valid': bool}`` for every screened tree.
    The (id, eta, valid) rows are cached so re-runs skip the loop.
    """
    cache_path = Path(cache_path)
    if cache_path.exists():
        try:
            z = np.load(cache_path, allow_pickle=True)
            rows = list(z["rows"])
            if (list(z["screen_ids"]) == list(screen_ids)
                    and rows and isinstance(rows[0], dict) and "valid" in rows[0]):
                print(f"loaded {label} from cache ({cache_path})")
                return _split(rows)
        except Exception:
            pass

    rows: List[dict] = []
    for k, tid in enumerate(screen_ids, 1):
        loaded = load_all(tid)
        if loaded is None:
            continue
        S, labels, tree, D = loaded
        part = np.asarray(partition_fn(S, D)).astype(bool)
        sidx = _tree_leaf_index(tree, labels)
        valid = bool(check_partition_valid_in_tree(tree, part[sidx]))
        rows.append({"tree": tid, "eta": float(_eta(part)), "valid": valid})
        if k % progress == 0:
            print(f"  {label}: screened {k}/{len(screen_ids)} "
                  f"(valid so far: {sum(r['valid'] for r in rows)})")
    np.savez(cache_path, screen_ids=np.array(list(screen_ids), dtype=object),
             rows=np.array(rows, dtype=object))
    print(f"saved {label} to cache ({cache_path})")
    return _split(rows)


def _split(rows: List[dict]) -> Tuple[List[str], np.ndarray, List[dict]]:
    valid_ids = [r["tree"] for r in rows if r["valid"]]
    etas_valid = np.array([r["eta"] for r in rows if r["valid"]], float)
    return valid_ids, etas_valid, rows
