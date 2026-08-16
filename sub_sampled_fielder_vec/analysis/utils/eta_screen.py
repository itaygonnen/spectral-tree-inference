"""Per-tree partition imbalance eta, for two operators, on generated trees.

For each tree: build S once, then read a bipartition off each of two operators and
record its imbalance eta = larger clan / smaller clan together with whether the
partition is a genuine single-edge bipartition of the true tree.

    L(S)  Fiedler of Deg(S) - S        -> k-means bipartition
    B     leading-|lambda| of H D H    -> sign rule

Both partitions come from ONE ``make_generated`` call per tree. ``run_screen`` in
:mod:`src.utils.screening` would re-simulate every tree once per operator, which at
seq_len=10^4 doubles a ~50 min job for nothing.

Cache: one ``.npz`` per configuration, named for ``(model, n, seq_len)``. Do NOT reuse
``run_screen``'s cache convention -- its hit test keys on the id list alone, so two
different ``n`` with the same ids silently return the first run's rows.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

MIN_SPLIT = 5  # matches the STDR/test convention used throughout the repo


def _eta(part: np.ndarray) -> float:
    """Imbalance of a boolean bipartition: larger clan / smaller clan."""
    n1 = int(np.sum(part))
    n2 = int(part.size) - n1
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


def screen_one(tid: str, seq_len: int) -> Optional[Dict]:
    """Both operators on one tree. Returns None when the tree cannot be built."""
    # imported here so each worker process pays the import once, not the parent
    from src.core.utils import (compute_fiedler_from_laplacian, compute_laplacian)
    from src.models.generated_trees import make_generated
    from src.runners.p_sweep_inner import _kmeans_bipartition
    from src.utils.griffing import griffing_leading_eigvec
    from src.utils.partition_validity import check_partition_valid_in_tree
    from src.utils.screening import _tree_leaf_index

    out = make_generated(tid, seq_len=seq_len)
    if out is None:
        return None
    S, labels, tree, D = out
    sidx = _tree_leaf_index(tree, labels)          # labels order -> tree-leaf order

    rec: Dict = {"tree": tid, "m": int(S.shape[0])}
    for tag, part in (
        ("S", _kmeans_bipartition(
            compute_fiedler_from_laplacian(compute_laplacian(S)), MIN_SPLIT)),
        ("B", griffing_leading_eigvec(D, solver="lm_k1") >= 0),
    ):
        part = np.asarray(part).astype(bool)
        rec[f"eta_{tag}"] = _eta(part)
        # the validity check needs the mask in TREE-LEAF order, not labels order
        rec[f"valid_{tag}"] = bool(check_partition_valid_in_tree(tree, part[sidx]))
    return rec


def _worker(args):
    tid, seq_len = args
    try:
        return screen_one(tid, seq_len)
    except Exception as exc:                       # one bad tree must not kill the run
        return {"tree": tid, "error": repr(exc)}


def run_eta_screen(
    ids: Sequence[str],
    cache_path: Path,
    *,
    seq_len: int = 10_000,
    workers: int = 6,
    force: bool = False,
    progress: int = 100,
) -> List[Dict]:
    """Screen ``ids``, caching to ``cache_path``. Resumable: only missing ids are run."""
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    done: Dict[str, Dict] = {}
    if cache_path.exists() and not force:
        z = np.load(cache_path, allow_pickle=True)
        done = {r["tree"]: r for r in z["rows"]}
        print(f"cache: {len(done)} rows in {cache_path.name}")

    todo = [t for t in ids if t not in done]
    if not todo:
        print("nothing to do -- all ids cached")
        return [done[t] for t in ids]

    print(f"computing {len(todo)} trees on {workers} workers (seq_len={seq_len})")
    n_done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_worker, (t, seq_len)): t for t in todo}
        for f in as_completed(futs):
            rec = f.result()
            if rec is not None:
                done[rec["tree"]] = rec
            n_done += 1
            if n_done % progress == 0 or n_done == len(todo):
                rows = [done[t] for t in ids if t in done]
                np.savez(cache_path, rows=np.array(rows, dtype=object))
                print(f"  [{n_done}/{len(todo)}] saved {len(rows)} rows", flush=True)

    rows = [done[t] for t in ids if t in done]
    bad = [r for r in rows if "error" in r]
    if bad:
        print(f"WARNING: {len(bad)} trees failed, e.g. {bad[0]['error'][:120]}")
    return rows


def main() -> None:
    import sys
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root))
    from src.utils.tree_ids import build_ids

    model = os.environ.get("ETA_MODEL", "kingman")
    n = int(os.environ.get("ETA_N", "1000"))
    seq_len = int(os.environ.get("ETA_L", "10000"))
    count = int(os.environ.get("ETA_COUNT", "1000"))
    workers = int(os.environ.get("ETA_WORKERS", "6"))

    ids = build_ids(models=[model], n_values=[n], per_cell=count)
    cache = (root / "analysis" / "notebooks_cache" / "eta_by_operator"
             / f"{model}_n{n}_L{seq_len}.npz")
    run_eta_screen(ids, cache, seq_len=seq_len, workers=workers)
    print("done ->", cache)


if __name__ == "__main__":
    main()
