"""Per-tree partition imbalance eta, for two operators, on a real cohort.

Same screen as :mod:`analysis.utils.eta_screen`, and the same record schema, but the
trees are the real FASTA/newick pairs of a :class:`analysis.utils.real_cohorts.Cohort`
rather than simulated ones. For each tree: load once, read a bipartition off each of two operators
and record its imbalance eta = larger clan / smaller clan together with whether the
partition is a genuine single-edge bipartition of the true tree.

    L(S)  Fiedler of Deg(S) - S        -> k-means bipartition
    B     leading-|lambda| of H D H    -> sign rule

Both arms come from ONE ``load_all`` call per tree; the load alone is ~1 min at
m=6000, so re-loading per operator (which ``src.utils.screening.run_screen`` does)
would double a ~2 h job for nothing. Unlike ``run_threshold_sweep`` this keeps eta for
*invalid* partitions too -- that is what the stacked valid/not-an-edge histogram needs.

Cache: one ``.npz`` of per-tree rows, resumable, keyed by the row list itself (rows are
looked up by tree id), so a longer id list extends an existing cache.
"""
from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from src.utils.logging import create_progress_bar, log_info, log_warning

MIN_SPLIT = 5  # matches the STDR/test convention used throughout the repo


def _eta(part: np.ndarray) -> float:
    """Imbalance of a boolean bipartition: larger clan / smaller clan."""
    n1 = int(np.sum(part))
    n2 = int(part.size) - n1
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


def screen_one(tid: str, cohort_name: str) -> Optional[Dict]:
    """Both operators on one tree. Returns None when the tree cannot be loaded."""
    # imported here so each worker process pays the import once, not the parent
    from analysis.utils.real_cohorts import get_cohort
    from src.core.utils import (compute_fiedler_from_laplacian, compute_laplacian)
    from src.runners.p_sweep_inner import _kmeans_bipartition
    from src.utils.griffing import griffing_leading_eigvec
    from src.utils.partition_validity import check_partition_valid_in_tree
    from src.utils.screening import _tree_leaf_index

    out = get_cohort(cohort_name).load_all(tid)
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
    tid, cohort_name = args
    try:
        return screen_one(tid, cohort_name)
    except Exception as exc:                       # one bad tree must not kill the run
        return {"tree": tid, "error": repr(exc)}


def run_real_eta_screen(
    ids: Sequence[str],
    cache_path,
    *,
    cohort_name: str = "6000 taxa",
    workers: int = 4,
    force: bool = False,
    progress: int = 5,
) -> List[Dict]:
    """Screen ``ids``, caching to ``cache_path``. Resumable: only missing ids are run.

    ``workers`` defaults to 4: each worker holds the (6000, 5000) alignment plus S, D
    and L at 288 MB apiece, so ~1.2 GB of resident memory per process.
    """
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    done: Dict[str, Dict] = {}
    if cache_path.exists() and not force:
        z = np.load(cache_path, allow_pickle=True)
        done = {r["tree"]: r for r in z["rows"]}
        log_info("screen", f"cache: {len(done)} rows in {cache_path}")

    todo = [t for t in ids if t not in done]
    if not todo:
        print(f"  all {len(ids)} trees already cached", flush=True)
        log_info("screen", f"{cohort_name}: nothing to do, all {len(ids)} ids cached")
        return [done[t] for t in ids]

    print(f"  {len(done)}/{len(ids)} trees cached, {len(todo)} to run "
          f"on {workers} workers", flush=True)
    log_info("screen", f"{cohort_name}: computing {len(todo)} trees on {workers} "
                       f"workers ({len(done)} already cached)")
    # One BLAS thread per worker: the per-tree work is already parallel across trees, and
    # on a many-core Linux box the default (every worker opening a full thread pool)
    # oversubscribes the machine and runs slower than serial. "spawn" so the children
    # inherit this env rather than a numpy already initialised with the old thread count
    # -- it also makes Linux and macOS behave identically here.
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(_v, "1")
    n_done = 0
    bar = create_progress_bar(len(todo), f"  screening {cohort_name}", unit="tree",
                              leave=False)
    with ProcessPoolExecutor(max_workers=workers,
                             mp_context=mp.get_context("spawn")) as ex:
        futs = {ex.submit(_worker, (t, cohort_name)): t for t in todo}
        for f in as_completed(futs):
            rec = f.result()
            if rec is not None:
                done[rec["tree"]] = rec
                if "error" in rec:
                    log_info("screen", f"{rec['tree']}: FAILED {rec['error'][:100]}")
                else:
                    log_info("screen",
                             f"{rec['tree']}: L(S) eta={rec['eta_S']:.2f} "
                             f"{'edge' if rec['valid_S'] else 'not-an-edge'}, "
                             f"B eta={rec['eta_B']:.2f} "
                             f"{'edge' if rec['valid_B'] else 'not-an-edge'}")
            n_done += 1
            bar.update(1)
            if n_done % progress == 0 or n_done == len(todo):
                rows = [done[t] for t in ids if t in done]
                np.savez(cache_path, rows=np.array(rows, dtype=object))
                log_info("screen", f"[{n_done}/{len(todo)}] saved {len(rows)} rows")
    bar.close()

    rows = [done[t] for t in ids if t in done]
    n_s = sum(1 for r in rows if r.get("valid_S"))
    n_b = sum(1 for r in rows if r.get("valid_B"))
    print(f"  {len(rows)} tree(s) screened: L(S) cuts a real edge on {n_s}, B on {n_b}",
          flush=True)
    bad = [r for r in rows if "error" in r]
    if bad:
        print(f"  WARNING: {len(bad)} trees failed, e.g. {bad[0]['error'][:120]}")
        log_warning("screen", f"{len(bad)} trees failed")
    return rows


def main() -> None:
    import sys
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root))
    from analysis.utils.real_cohorts import get_cohort, screen_cache_path

    name = os.environ.get("REAL_COHORT", "6000 taxa")
    workers = int(os.environ.get("REAL_WORKERS", "4"))
    limit = int(os.environ.get("REAL_LIMIT", "0")) or None
    cohort = get_cohort(name)
    cache = screen_cache_path(name)
    run_real_eta_screen(cohort.ids(limit), cache,
                        cohort_name=name, workers=workers)
    print("done ->", cache)


if __name__ == "__main__":
    main()
