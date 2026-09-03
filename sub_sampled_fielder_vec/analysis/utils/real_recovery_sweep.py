"""Sub-sampling recovery sweep on a real cohort, two operators.

For one tree the sweep asks: replace the full matrix by a uniform sub-sample at rate
``p`` and re-read the bipartition -- how much of the full-matrix partition survives?
Both arms are scored by NMI against their OWN full-matrix reference:

    L   Fiedler of Deg(S) - S, k-means cut      (``bootstrap_p_sweep_simple``)
    B   leading-|lambda| of H D H, sign cut     (``bpart_sweep_raw``)

Two settings are not the library defaults and both are load-bearing here:

* ``eigsolver="lm_k1"`` -- ARPACK for the single eigenpair B needs. Same vector as the
  default ``eigh_full`` (verified in :mod:`src.utils.griffing`) at 0.39 s vs 17.2 s per
  solve at m=6000, which is what makes a 20-point x 10-rep grid finishable at all.
* ``aggregation="avg_vector"`` -- bootstrap-average the eigenvector, THEN partition
  once. That is what ``bootstrap_p_sweep_simple`` does on the L arm, so the two curves
  are comparable; the ``per_rep`` default would score B under a different rule.

One ``.npz`` per tree, so the run is resumable and interruptible at any tree boundary
(the full 100-tree grid is a ~50 h job dominated by the L arm at ~9 s per solve).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Cache-key fields. A tree whose stored meta differs is recomputed rather than mixed
# into a figure with a different grid.
META_KEYS = ("reps", "num_gaps", "min_split", "partition_method",
             "eigsolver", "aggregation", "m")


def sweep_meta(reps: int, num_gaps: int, min_split: int, m: int = 6000) -> Dict:
    return dict(reps=int(reps), num_gaps=int(num_gaps), min_split=int(min_split),
                partition_method="kmeans", eigsolver="lm_k1",
                aggregation="avg_vector", m=int(m))


def _cache_file(cache_dir: Path, tid: str) -> Path:
    return Path(cache_dir) / f"{tid}.npz"


def seed_for(tid: str, stride: int = 1000) -> int:
    """Seed from the tree's own number, not its position in some cohort list, so a
    pilot run, a partial run and the full run all produce the same curve for a tree."""
    return stride * int(tid.rsplit("_", 1)[-1])


def load_tree_sweep(cache_dir, tid: str, p_values: Sequence[float],
                    meta: Dict) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Return ``(nmi_L, nmi_B)`` for one tree if a matching cache exists, else None."""
    path = _cache_file(Path(cache_dir), tid)
    if not path.exists():
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if list(np.round(z["p_values"], 6)) != list(np.round(p_values, 6)):
            return None
        if {k: dict(z["meta"].item())[k] for k in META_KEYS} != {k: meta[k] for k in META_KEYS}:
            return None
        return np.asarray(z["nmi_L"], float), np.asarray(z["nmi_B"], float)
    except Exception:
        return None


def sweep_one_tree(tid: str, seed: int, p_values: Sequence[float], *, reps: int,
                   num_gaps: int, min_split: int,
                   cohort_name: str = "6000 taxa") -> Tuple[np.ndarray, np.ndarray]:
    """Run both arms on one tree. Returns ``(nmi_L, nmi_B)``, one value per ``p``."""
    from analysis.utils.real_cohorts import get_cohort
    from src.core.utils import compute_fiedler_from_laplacian, compute_laplacian
    from src.runners.p_sweep_inner import bootstrap_p_sweep_simple
    from src.utils.bpart_sweep_cache import bpart_sweep_raw

    loaded = get_cohort(cohort_name).load_all(tid)
    if loaded is None:
        raise FileNotFoundError(f"{tid}: alignment or tree missing")
    S, _labels, _tree, D = loaded

    fr_L = compute_fiedler_from_laplacian(compute_laplacian(S))
    out_L = bootstrap_p_sweep_simple(
        S, fr_L, list(p_values), bootstrap_reps=reps, seed=seed,
        num_gaps=num_gaps, min_split=min_split,
        partition_method="kmeans", laplacian="unnormalized")
    nmi_L = np.asarray(out_L["partition_nmi_M"], float)

    raw = bpart_sweep_raw(D, list(p_values), reps=reps, seed_base=seed,
                          eigsolver="lm_k1", aggregation="avg_vector")
    nmi_B = np.array([float(np.mean(pp["nmi"])) for pp in raw["per_p"]], float)
    return nmi_L, nmi_B


def run_sweep(ids: Sequence[str], cache_dir, p_values: Sequence[float], *,
              reps: int = 10, num_gaps: int = 10, min_split: int = 5,
              cohort_name: str = "6000 taxa", m: int = 6000,
              seed_stride: int = 1000) -> List[str]:
    """Sweep every id in ``ids``, skipping trees already cached. Returns ids done."""
    import time

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = sweep_meta(reps, num_gaps, min_split, m)
    pv = np.round(np.asarray(p_values, float), 6)

    finished: List[str] = []
    todo = [t for t in ids if load_tree_sweep(cache_dir, t, pv, meta) is None]
    print(f"{len(ids) - len(todo)}/{len(ids)} trees already cached; "
          f"{len(todo)} to run ({len(pv)} p x {reps} reps)", flush=True)

    t0 = time.time()
    for k, tid in enumerate(todo, 1):
        t_tree = time.time()
        nmi_L, nmi_B = sweep_one_tree(
            tid, seed=seed_for(tid, seed_stride), p_values=pv,
            reps=reps, num_gaps=num_gaps, min_split=min_split,
            cohort_name=cohort_name)
        np.savez(_cache_file(cache_dir, tid), p_values=pv, nmi_L=nmi_L, nmi_B=nmi_B,
                 meta=np.array(meta, dtype=object))
        finished.append(tid)
        el, per = time.time() - t0, (time.time() - t0) / k
        print(f"  [{k}/{len(todo)}] {tid}  {time.time() - t_tree:.0f} s  "
              f"(elapsed {el/3600:.1f} h, ETA {per*(len(todo)-k)/3600:.1f} h)", flush=True)
    return finished


def collect(ids: Sequence[str], cache_dir, p_values: Sequence[float],
            meta: Dict) -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    """Stack the cached curves. Returns ``(nmi_L, nmi_B, have, pending)``."""
    pv = np.round(np.asarray(p_values, float), 6)
    have, pending, L, B = [], [], [], []
    for tid in ids:
        got = load_tree_sweep(cache_dir, tid, pv, meta)
        if got is None:
            pending.append(tid)
            continue
        have.append(tid)
        L.append(got[0])
        B.append(got[1])
    nmi_L = np.vstack(L) if L else np.empty((0, len(pv)))
    nmi_B = np.vstack(B) if B else np.empty((0, len(pv)))
    return nmi_L, nmi_B, have, pending
