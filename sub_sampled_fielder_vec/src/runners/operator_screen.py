"""Per-tree partition imbalance and validity, for any operator, from any source.

For each tree: load once, read a bipartition off every selected operator, and record
its imbalance ``eta`` = larger clan / smaller clan together with whether the partition
is a genuine single-edge bipartition of the true tree. Both cut rules are scored on
every Fiedler operator and the more balanced one wins -- see :mod:`src.runners.operators`.

Where the tree comes from is the ``loader``'s business (``RealLoader`` for a FASTA
directory, ``GeneratedLoader`` for a simulated id), so this module has no idea whether
it is screening real alignments or simulated ones. That is the whole point: the screen
used to exist twice, once per data source, and the two drifted apart.

One load per tree serves every operator. The load alone is ~1 min at m=6000, so
re-loading per operator (which ``src.utils.screening.run_screen`` does) would multiply
a 2 h job by the number of operators for nothing.

Cache: one ``.npz`` of per-tree rows, resumable, looked up by tree id. A row written
before an operator existed is **topped up**, not thrown away -- it is re-run and its
other operators land on the same values, because nothing about them changed.
"""
from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..utils.logging import create_progress_bar, log_info, log_warning
from .operators import MIN_SPLIT, OPERATORS, best_cut, resolve

_LOADER = None                     # set per worker process by _init


def _init(loader) -> None:
    global _LOADER
    _LOADER = loader


def screen_one(tree_id: str, loader, operators: Sequence[str] = (),
               min_split: int = MIN_SPLIT) -> Optional[Dict]:
    """Every selected operator on one tree. None when the tree cannot be loaded."""
    from ..utils.partition_metrics import eta as _eta
    from ..utils.partition_validity import check_partition_valid_in_tree
    from ..utils.screening import _tree_leaf_index

    loaded = loader(tree_id)
    if loaded is None:
        return None
    S, labels, tree, D = loaded
    # the validity check needs the mask in TREE-LEAF order, not labels order
    sidx = _tree_leaf_index(tree, labels) if tree is not None else None

    rec: Dict = {"tree": tree_id, "m": int(S.shape[0])}
    for key in resolve(operators):
        op = OPERATORS[key]
        try:
            vec = op.vector(S, D)
        except Exception as exc:                  # one operator must not cost the tree
            log_warning("screen", f"{tree_id}: {key} failed: {exc!r}")
            continue
        rule, part, etas = best_cut(op, vec, min_split)
        for r, e in etas.items():
            rec[f"eta_{key}_{r}"] = float(e)
            if sidx is not None:
                from .operators import cut
                rec[f"valid_{key}_{r}"] = bool(
                    check_partition_valid_in_tree(tree, cut(vec, r, min_split)[sidx]))
        rec[f"rule_{key}"] = rule
        rec[f"eta_{key}"] = float(_eta(part))
        if sidx is not None:
            rec[f"valid_{key}"] = bool(check_partition_valid_in_tree(tree, part[sidx]))
    return rec


def _worker(args):
    tree_id, operators, min_split = args
    try:
        return screen_one(tree_id, _LOADER, operators, min_split)
    except Exception as exc:                       # one bad tree must not kill the run
        return {"tree": tree_id, "error": repr(exc)}


def _needs_rerun(row: Dict, operators: Sequence[str]) -> bool:
    """True when a cached row predates one of the operators (or the two-rule cut)."""
    if "error" in row:
        return False                               # keep the failure; re-run is futile
    return any(f"rule_{k}" not in row for k in operators)


def run_screen(
    ids: Sequence[str],
    cache_path,
    loader,
    *,
    operators: Sequence[str] = (),
    source: str = "",
    workers: int = 4,
    min_split: int = MIN_SPLIT,
    force: bool = False,
    progress: int = 1,
) -> List[Dict]:
    """Screen ``ids`` through ``loader``, caching to ``cache_path``. Resumable.

    ``workers`` defaults to 4: each worker holds the alignment plus S, D and L, which
    at m=6000 is ~1.2 GB of resident memory per process -- drop to 2 on a machine with
    other work on it, or the OS starts killing workers.

    ``progress`` is how many completed trees to accumulate before rewriting the cache.
    It is 1 because the cost of a rewrite (a few KB) is nothing beside a tree, and a run
    killed for memory otherwise loses everything since the last checkpoint.
    """
    ops = resolve(operators)
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    done: Dict[str, Dict] = {}
    if cache_path.exists() and not force:
        z = np.load(cache_path, allow_pickle=True)
        done = {r["tree"]: r for r in z["rows"]}
        log_info("screen", f"cache: {len(done)} rows in {cache_path}")

    # A row that predates an operator is re-run, but NOT dropped first: the old
    # verdicts stay in the cache until their replacement arrives, so a top-up killed
    # half-way leaves the screen exactly as complete as it was before it started.
    stale = [t for t in ids if t in done and _needs_rerun(done[t], ops)]
    if stale:
        missing = sorted({k for t in stale for k in ops if f"rule_{k}" not in done[t]})
        print(f"  {len(stale)} cached row(s) predate {', '.join(missing)} "
              f"-- topping them up (the existing verdicts stay until then)", flush=True)
        log_info("screen", f"{len(stale)} row(s) re-run for {missing}", force=True)
    todo = [t for t in ids if t not in done or t in set(stale)]
    if not todo:
        print(f"  all {len(ids)} trees already cached", flush=True)
        log_info("screen", f"{source}: nothing to do, all {len(ids)} ids cached")
        return [done[t] for t in ids if t in done]

    print(f"  {len(done)}/{len(ids)} trees cached, {len(todo)} to run "
          f"on {workers} workers ({', '.join(ops)})", flush=True)
    log_info("screen", f"{source}: computing {len(todo)} trees on {workers} "
                       f"workers, operators {ops} ({len(done)} already cached)")
    # One BLAS thread per worker: the per-tree work is already parallel across trees,
    # and on a many-core Linux box the default (every worker opening a full thread
    # pool) oversubscribes the machine and runs slower than serial. "spawn" so the
    # children inherit this env rather than a numpy already initialised with the old
    # thread count -- it also makes Linux and macOS behave identically here.
    for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(_v, "1")

    n_done = 0
    bar = create_progress_bar(len(todo), f"  screening {source}".rstrip(),
                              unit="tree", leave=False)
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"),
                             initializer=_init, initargs=(loader,)) as ex:
        futs = {ex.submit(_worker, (t, ops, min_split)): t for t in todo}
        for f in as_completed(futs):
            rec = f.result()
            if rec is not None:
                done[rec["tree"]] = rec
                _log_row(rec, ops)
            n_done += 1
            bar.update(1)
            if n_done % progress == 0 or n_done == len(todo):
                rows = [done[t] for t in ids if t in done]
                np.savez(cache_path, rows=np.array(rows, dtype=object))
                log_info("screen", f"[{n_done}/{len(todo)}] saved {len(rows)} rows")
    bar.close()

    rows = [done[t] for t in ids if t in done]
    _report(rows, ops)
    return rows


def _log_row(rec: Dict, ops: Sequence[str]) -> None:
    if "error" in rec:
        log_info("screen", f"{rec['tree']}: FAILED {rec['error'][:100]}")
        return
    parts = []
    for k in ops:
        if f"eta_{k}" not in rec:
            continue
        edge = "edge" if rec.get(f"valid_{k}") else "not-an-edge"
        parts.append(f"{k} eta={rec[f'eta_{k}']:.2f}/{edge} ({rec.get(f'rule_{k}')})")
    log_info("screen", f"{rec['tree']}: " + ", ".join(parts))


def _report(rows: List[Dict], ops: Sequence[str]) -> None:
    bad = [r for r in rows if "error" in r]
    ok = [r for r in rows if "error" not in r]
    if ok:
        counts = ", ".join(f"{k} on {sum(1 for r in ok if r.get(f'valid_{k}'))}"
                           for k in ops)
        print(f"  {len(ok)} tree(s) screened; cuts a real edge: {counts}", flush=True)
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
