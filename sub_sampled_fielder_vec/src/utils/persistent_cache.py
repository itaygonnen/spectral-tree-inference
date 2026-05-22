"""Persistent disk-based caching for experiment data.

Thin shim over :mod:`src.cache_io`. Backing scope is ``experiment_data``;
all keys go through :func:`src.cache_io.make_key`. The public API is
preserved so existing callers keep working unchanged.

Stored artifacts (filenames unchanged for byte-compat with migrated entries):
    ``tree.npz``  (multi-array NPZ — ``adjacency_matrix`` or ``newick`` keys)
    ``observations.npz`` (single-array NPZ, key ``observations``)
    ``similarity_matrix.npz`` (single-array NPZ, key ``similarity_matrix``)
    ``fiedler_ref.npz`` (single-array NPZ, key ``fiedler_ref``)
    ``metadata.json``
    ``.complete`` (sentinel)
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .logging import log_info, log_warning
from ..cache_io import experiment_data as _scope, make_key as _make_key, SENTINEL


def _get_cache_key(
    n_taxa: int,
    seq_len: int,
    mutation_rate: float,
    tree_model_name: str,
    seq_model_name: str,
    tree_params: Optional[Dict[str, Any]] = None,
    seq_params: Optional[Dict[str, Any]] = None,
    matrix_kind: str = "similarity",
    distance_alpha: float = 1.0,
) -> str:
    """Build the canonical experiment_data key for these parameters."""
    kwargs: Dict[str, Any] = {
        "n":    int(n_taxa),
        "L":    int(seq_len),
        "mu":   float(mutation_rate),
        "tree": str(tree_model_name),
        "seq":  str(seq_model_name),
    }
    if tree_params:
        for k, v in tree_params.items():
            if k == "num_taxa":
                continue
            kwargs[k] = v
    if seq_params:
        for k, v in seq_params.items():
            if k == "mutation_rate":
                continue
            kwargs[k] = v
    if matrix_kind == "distance":
        kwargs["matrix_kind"] = "distance"
        kwargs["distance_alpha"] = float(distance_alpha)
    return _make_key(**kwargs)


def _get_cache_dir(cache_key: str) -> Path:
    """Path under experiment_data scope for this key. Creates it if missing."""
    d = _scope.path(cache_key)
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_experiment_data(
    cache_key: str,
    tree: Any,
    observations: np.ndarray,
    similarity_matrix: np.ndarray,
    fiedler_ref: np.ndarray,
    metadata: Dict[str, Any],
) -> None:
    """Save (tree, observations, similarity_matrix, fiedler_ref, metadata) atomically."""
    d = _scope.path(cache_key)
    d.mkdir(parents=True, exist_ok=True)
    sentinel = d / SENTINEL
    if sentinel.exists():
        sentinel.unlink()
    try:
        tree_data = _serialize_tree(tree)
        np.savez_compressed(d / "tree.npz", **tree_data)
        np.savez_compressed(d / "observations.npz", observations=observations)
        np.savez_compressed(d / "similarity_matrix.npz", similarity_matrix=similarity_matrix)
        np.savez_compressed(d / "fiedler_ref.npz", fiedler_ref=fiedler_ref)
        (d / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str))
        sentinel.touch()
        log_info("cache", f"Saved experiment data to cache: {cache_key}")
    except Exception as e:
        log_warning("cache", f"Failed to save cache {cache_key}: {e}")


def load_experiment_data(cache_key: str) -> Optional[Dict[str, Any]]:
    """Load all five artifacts. Returns None on miss or corruption."""
    d = _scope.path(cache_key)
    required = [".complete", "tree.npz", "observations.npz",
                "similarity_matrix.npz", "fiedler_ref.npz", "metadata.json"]
    for f in required:
        if not (d / f).exists():
            log_info("cache", f"Cache miss: {cache_key} (missing {f})")
            return None
    try:
        tree_npz = np.load(d / "tree.npz", allow_pickle=True)
        tree = _deserialize_tree(tree_npz)
        observations = np.load(d / "observations.npz")["observations"]
        similarity_matrix = np.load(d / "similarity_matrix.npz")["similarity_matrix"]
        fiedler_ref = np.load(d / "fiedler_ref.npz")["fiedler_ref"]
        metadata = json.loads((d / "metadata.json").read_text())
        log_info("cache", f"Cache hit: {cache_key}")
        return {
            "tree": tree,
            "observations": observations,
            "similarity_matrix": similarity_matrix,
            "fiedler_ref": fiedler_ref,
            "metadata": metadata,
        }
    except Exception as e:
        log_warning("cache", f"Failed to load cache {cache_key}: {e}")
        return None


def list_cached_experiments() -> List[Dict[str, Any]]:
    """List all complete entries under the experiment_data scope."""
    out: List[Dict[str, Any]] = []
    if not _scope.root.exists():
        return out
    for d in sorted(_scope.root.iterdir()):
        if not d.is_dir() or not (d / SENTINEL).exists():
            continue
        meta_path = d / "metadata.json"
        if not meta_path.exists():
            continue
        try:
            metadata = json.loads(meta_path.read_text())
            out.append({"cache_key": d.name, "metadata": metadata, "path": str(d)})
        except Exception as e:
            log_warning("cache", f"Failed to read metadata for {d.name}: {e}")
    return out


def clear_cache(cache_key: Optional[str] = None) -> None:
    """Remove a specific entry or the whole experiment_data scope."""
    if cache_key is None:
        _scope.clear()
        log_info("cache", "Cleared all caches")
        return
    d = _scope.path(cache_key)
    if d.exists():
        shutil.rmtree(d)
        log_info("cache", f"Cleared cache: {cache_key}")
    else:
        log_info("cache", f"Cache not found: {cache_key}")


def clean_incomplete_caches() -> int:
    """Remove every cache entry under experiment_data that has no sentinel."""
    if not _scope.root.exists():
        return 0
    removed = 0
    for d in list(_scope.root.iterdir()):
        if not d.is_dir():
            continue
        if not (d / SENTINEL).exists():
            try:
                shutil.rmtree(d)
                log_info("cache", f"Removed incomplete cache: {d.name}")
                removed += 1
            except Exception as e:
                log_warning("cache", f"Failed to remove incomplete cache {d.name}: {e}")
    if removed > 0:
        log_info("cache", f"Cleaned up {removed} incomplete cache entries")
    return removed


def _serialize_tree(tree: Any) -> Dict[str, np.ndarray]:
    if hasattr(tree, "adjacency_matrix"):
        return {"adjacency_matrix": tree.adjacency_matrix}
    if hasattr(tree, "newick"):
        newick_str = tree.newick() if callable(tree.newick) else str(tree.newick)
        return {"newick": np.array([newick_str], dtype=object)}
    return {"tree_str": np.array([str(tree)], dtype=object)}


def _deserialize_tree(tree_npz) -> Any:
    import dendropy
    if "adjacency_matrix" in tree_npz:
        return tree_npz["adjacency_matrix"]
    if "newick" in tree_npz:
        newick_str = str(tree_npz["newick"][0])
        return dendropy.Tree.get(data=newick_str, schema="newick")
    if "tree_str" in tree_npz:
        tree_str = str(tree_npz["tree_str"][0])
        try:
            return dendropy.Tree.get(data=tree_str, schema="newick")
        except Exception:
            return tree_str
    return None
