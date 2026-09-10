"""Per-tree partition imbalance eta, for two operators, on a real cohort.

Same screen as :mod:`analysis.utils.eta_screen`, and the same record schema, but the
trees are the real FASTA/newick pairs of a :class:`analysis.utils.real_cohorts.Cohort`
rather than simulated ones. For each tree: load once, read a bipartition off each of two operators
and record its imbalance eta = larger clan / smaller clan together with whether the
partition is a genuine single-edge bipartition of the true tree.

    L(S)  Fiedler of Deg(S) - S        -> k-means AND sign, whichever cuts more evenly
    B     leading-|lambda| of H D H    -> sign rule

The L arm is scored under both threshold rules and the more balanced one wins. k-means
alone routinely isolates a single taxon on real data (1/999, eta=999) -- a real pendant
edge, so a validity gate passes it, but a split with no recovery signal. Both rules'
eta and validity are recorded, so the choice is auditable and either can be re-read.

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
    fr = compute_fiedler_from_laplacian(compute_laplacian(S))
    candidates = {
        "kmeans": np.asarray(_kmeans_bipartition(fr, MIN_SPLIT)).astype(bool),
        "sign": np.asarray(fr >= 0).astype(bool),
    }
    for rule, part in candidates.items():
        rec[f"eta_S_{rule}"] = _eta(part)
        rec[f"valid_S_{rule}"] = bool(check_partition_valid_in_tree(tree, part[sidx]))
    # the more balanced rule wins; ties go to k-means, the historical default
    rec["rule_S"] = min(candidates, key=lambda r: (rec[f"eta_S_{r}"], r != "kmeans"))
    rec["eta_S"] = rec[f"eta_S_{rec['rule_S']}"]
    rec["valid_S"] = rec[f"valid_S_{rec['rule_S']}"]

    part_b = np.asarray(griffing_leading_eigvec(D, solver="lm_k1") >= 0).astype(bool)
    rec["eta_B"] = _eta(part_b)
    # the validity check needs the mask in TREE-LEAF order, not labels order
    rec["valid_B"] = bool(check_partition_valid_in_tree(tree, part_b[sidx]))
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
    progress: int = 1,
) -> List[Dict]:
    """Screen ``ids``, caching to ``cache_path``. Resumable: only missing ids are run.

    ``workers`` defaults to 4: each worker holds the (6000, 5000) alignment plus S, D
    and L at 288 MB apiece, so ~1.2 GB of resident memory per process -- drop to 2 on a
    machine with other work on it, or the OS starts killing workers.

    ``progress`` is how many completed trees to accumulate before rewriting the cache.
    It is 1 because the cost of a rewrite (a few KB) is nothing beside a tree, and a run
    killed for memory otherwise loses everything since the last checkpoint.
    """
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    done: Dict[str, Dict] = {}
    if cache_path.exists() and not force:
        z = np.load(cache_path, allow_pickle=True)
        done = {r["tree"]: r for r in z["rows"]}
        log_info("screen", f"cache: {len(done)} rows in {cache_path}")

    # rows written before the L arm gained the sign/k-means comparison lack rule_S, and
    # silently keeping them would mix two definitions of eta_S in one screen
    stale = [t for t, r in done.items() if "error" not in r and "rule_S" not in r]
    if stale:
        log_info("screen", f"{len(stale)} cached row(s) predate the two-rule L arm "
                           f"-- recomputing them", force=True)
        for t in stale:
            done.pop(t, None)
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
                             f"{rec['tree']}: L(S) k-means eta="
                             f"{rec['eta_S_kmeans']:.2f}"
                             f"{'/edge' if rec['valid_S_kmeans'] else '/not-an-edge'}, "
                             f"sign eta={rec['eta_S_sign']:.2f}"
                             f"{'/edge' if rec['valid_S_sign'] else '/not-an-edge'} "
                             f"-> {rec['rule_S']}; B eta={rec['eta_B']:.2f}"
                             f"{'/edge' if rec['valid_B'] else '/not-an-edge'}")
            n_done += 1
            bar.update(1)
            if n_done % progress == 0 or n_done == len(todo):
                rows = [done[t] for t in ids if t in done]
                np.savez(cache_path, rows=np.array(rows, dtype=object))
                log_info("screen", f"[{n_done}/{len(todo)}] saved {len(rows)} rows")
    bar.close()

    rows = [done[t] for t in ids if t in done]
    bad = [r for r in rows if "error" in r]
    ok = [r for r in rows if "error" not in r]
    n_s = sum(1 for r in ok if r.get("valid_S"))
    n_b = sum(1 for r in ok if r.get("valid_B"))
    if ok:
        print(f"  {len(ok)} tree(s) screened: L(S) cuts a real edge on {n_s}, "
              f"B on {n_b}", flush=True)
    if bad:
        # one full message, not a truncated one repeated per tree: these failures are
        # nearly always one cause (unaligned FASTA, missing tree) shared by every tree
        print(f"  {len(bad)} of {len(rows)} tree(s) FAILED. First error in full:",
              flush=True)
        print(f"    {bad[0]['tree']}: {bad[0]['error']}", flush=True)
        log_warning("screen", f"{len(bad)}/{len(rows)} trees failed: "
                              f"{bad[0]['error']}", force=True)
        if not ok:
            print("  nothing was screened -- fix the cause above before sweeping",
                  flush=True)
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
