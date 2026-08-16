"""Shared spine of the two distance-vs-similarity notebooks.

``analysis/supporting/distance_vs_similarity_{real,generated}.ipynb`` run the
same pipeline over different data. ``screening.run_screen`` already collapsed
their three per-operator screen cells; this module collapses what came after —
the operator x threshold sweep and the two per-tree helpers that both notebooks
had defined inline, byte for byte.

The operator and threshold callables are **parameters**, not constants defined
here: the sigma2 rule is ``spectraltree.partition_taxa``, and ``src/`` does not
import ``spectraltree`` (see CLAUDE.md, cross-project boundary). The notebooks
build the two dicts, where that import already lives.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .partition_validity import check_partition_valid_in_tree


def tree_leaf_index(tree, labels) -> np.ndarray:
    """Permutation taking ``tree``'s leaf order to the order of ``labels``."""
    idx = {l: i for i, l in enumerate(labels)}
    return np.array([idx[str(l.taxon.label)] for l in tree.leaf_node_iter()])


def patristic_distance_matrix(tree) -> np.ndarray:
    """Dense leaf-by-leaf patristic distances, in ``tree``'s own leaf order."""
    leaf_taxa = [leaf.taxon for leaf in tree.leaf_node_iter()]
    n = len(leaf_taxa)
    pdm = tree.phylogenetic_distance_matrix()
    Dp = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = pdm.distance(leaf_taxa[i], leaf_taxa[j])
            Dp[i, j] = Dp[j, i] = d
    return Dp


def p_star_per_tree(nmi_arr, p_values, thresh: float = 0.95) -> np.ndarray:
    """Smallest ``p`` whose NMI clears ``thresh``, per row; ``inf`` if none does.

    ``inf`` is a NON-DETECTION, not a threshold — the operator never recovered
    on this grid. Do not fit a rate through those points.
    """
    out = []
    for row in np.asarray(nmi_arr, float):
        ok = np.where(row >= thresh)[0]
        out.append(float(p_values[ok[0]]) if ok.size else np.inf)
    return np.array(out)


def run_threshold_sweep(
    ids: Sequence[str],
    load_all: Callable[[str], Optional[Tuple]],
    operators: Dict[str, Callable[[np.ndarray, np.ndarray], np.ndarray]],
    thresholds: Dict[str, Callable[[np.ndarray, np.ndarray], np.ndarray]],
    cache_path,
    meta: dict,
    *,
    progress: int = 50,
) -> Tuple[Dict[Tuple[str, str], np.ndarray], int]:
    """Partition imbalance eta per (operator, threshold), valid partitions only.

    For every tree in ``ids``: build a vector with each operator, cut it with
    each threshold rule, and keep the resulting eta only when the partition is a
    real single-edge clan of the true tree. Trees that fail to load are skipped.

    Parameters
    ----------
    load_all(tid) -> (S, labels, tree, D) or None
        ``D`` aligned to ``labels`` order. ``None`` skips a missing tree.
    operators / thresholds
        ``op(S, D) -> vector`` and ``thr(vector, S) -> boolean partition``.
    meta
        Cache-key dict. A cached file whose ``meta`` differs is ignored, so
        changing ``min_split`` or the id list recomputes rather than silently
        reusing a stale sweep.

    Returns ``({(op, thr): etas}, n_screened)``.
    """
    op_names, thr_names = list(operators), list(thresholds)
    cache_path = Path(cache_path)

    if cache_path.exists():
        z = np.load(cache_path, allow_pickle=True)
        if dict(z["meta"].item()) == meta:
            etas = {tuple(k.split("|")): np.asarray(v, float)
                    for k, v in z["etas"].item().items()}
            # The generated notebook's cache predates n_screened; it loads every id,
            # so len(ids) is the value it would have written.
            n_screened = int(z["n_screened"]) if "n_screened" in z.files else len(ids)
            print(f"loaded threshold-eta cache ({cache_path})")
            return etas, n_screened

    acc: Dict[Tuple[str, str], List[float]] = {
        (op, thr): [] for op in op_names for thr in thr_names
    }
    n_screened = 0
    for k, tid in enumerate(ids, 1):
        loaded = load_all(tid)
        if loaded is None:
            continue
        S, labels, tree, D = loaded
        n_screened += 1
        sidx = tree_leaf_index(tree, labels)
        for op in op_names:
            v = operators[op](S, D)
            for thr in thr_names:
                try:
                    part = np.asarray(thresholds[thr](v, S)).astype(bool)
                except Exception:
                    continue
                if not bool(check_partition_valid_in_tree(tree, part[sidx])):
                    continue
                n1 = int(part.sum())
                n2 = len(part) - n1
                acc[(op, thr)].append(max(n1, n2) / max(min(n1, n2), 1))
        if progress and k % progress == 0:
            print(f"  screened {k}/{len(ids)}")

    etas = {key: np.asarray(val, float) for key, val in acc.items()}
    np.savez(
        cache_path,
        meta=meta,
        n_screened=n_screened,
        etas={f"{op}|{thr}": etas[(op, thr)] for op in op_names for thr in thr_names},
    )
    print(f"saved threshold-eta cache ({cache_path})")
    return etas, n_screened


def print_validity_table(
    etas: Dict[Tuple[str, str], np.ndarray],
    op_names: Sequence[str],
    thr_names: Sequence[str],
    n_screened: int,
) -> None:
    """Valid-partition counts as an operator-by-threshold table."""
    print(f"\nvalid partitions / {n_screened} screened:")
    print(f"  {'op':>6} | " + " | ".join(f"{t:>8}" for t in thr_names))
    for op in op_names:
        print(f"  {op:>6} | " + " | ".join(f"{len(etas[(op, t)]):>8d}" for t in thr_names))
