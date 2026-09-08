"""Sub-sampling recovery sweep on a real cohort, two operators.

For one tree the sweep asks: replace the full matrix by a uniform sub-sample at rate
``p`` and re-read the bipartition -- how much of the full-matrix partition survives?
Both arms are scored by NMI against their OWN full-matrix reference:

    L   Fiedler of Deg(S) - S, k-means cut      (``bootstrap_p_sweep_simple``)
    B   leading-|lambda| of H D H, sign cut     (``bpart_sweep_raw``)

Every metric both arms produce is stored per tree -- NMI, ARI, partition agreement, sign
agreement and the dot product with the reference vector. NMI is what the recovery figure
plots; the others are free once the eigenvector exists and answer the next question
without a 50 h recompute.

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
# into a figure with a different grid. ``metrics`` is part of the key so caches written
# before every metric was stored are recomputed instead of read back half-empty.
META_KEYS = ("reps", "num_gaps", "min_split", "partition_method",
             "eigsolver", "aggregation", "m", "metrics")

# Metrics BOTH arms produce per p, so a figure can be drawn on any of them. NMI is what the
# recovery figure plots; the rest cost nothing once the eigenvector is in hand, so they are
# stored rather than thrown away. ``agreement`` and ``dot`` are the orientation-invariant
# clan-label match (%) and the dot product with the reference vector.
METRICS = ("nmi", "ari", "agreement", "dot")

# The same three scores against a SECOND reference: the true tree's own top bipartition
# (its root's two leaf sets), which is what the earlier real-data benchmark reported as
# ``*_gt``. Recovery towards the full matrix says the sub-sample kept what the full matrix
# saw; recovery towards the tree says the full matrix was seeing the right thing. They can
# disagree sharply -- 0.68 vs 0.013 on the m=1000 run -- so both are stored.
GT_METRICS = ("nmi_gt", "ari_gt", "agreement_gt")

# L-only extra: ``bootstrap_p_sweep_simple`` reports sign agreement separately, while on the
# B arm the sign pattern IS the partition, so its ``agreement`` already is that number.
EXTRA_L = ("sign",)


def sweep_meta(reps: int, num_gaps: int, min_split: int, m: int = 6000) -> Dict:
    return dict(reps=int(reps), num_gaps=int(num_gaps), min_split=int(min_split),
                partition_method="kmeans", eigsolver="lm_k1",
                aggregation="avg_vector", m=int(m),
                metrics=",".join(METRICS + GT_METRICS))


def _gt_partition(tree, labels) -> np.ndarray:
    """The true tree's top bipartition, in ``labels`` order.

    The root's two child subtrees, matching ``src.runners.real_data_bpart`` so the numbers
    are comparable with the earlier benchmark. Falls back to the first internal node with
    two children when the seed node is degenerate.
    """
    root = tree.seed_node
    children = list(root.child_nodes())
    if len(children) < 2:
        for node in tree.preorder_node_iter():
            if len(list(node.child_nodes())) >= 2:
                children = list(node.child_nodes())
                break
    left = {str(leaf.taxon.label) for leaf in children[0].leaf_iter()}
    part = np.array([l in left for l in labels], dtype=bool)
    n_left = int(part.sum())
    if n_left in (0, part.size):            # degenerate: keep NMI defined
        part[:part.size // 2] = True
    return part


def _score(ref: np.ndarray, pred) -> Dict[str, float]:
    """NMI, ARI and clan agreement of one partition against a reference."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    if pred is None:
        return dict(nmi=float("nan"), ari=float("nan"), agreement=float("nan"))
    a, b = np.asarray(ref).astype(int), np.asarray(pred).astype(int)
    agr = max((a == b).mean(), (a != b).mean()) * 100.0   # orientation-invariant
    return dict(nmi=float(normalized_mutual_info_score(a, b)),
                ari=float(adjusted_rand_score(a, b)), agreement=float(agr))


def _cache_file(cache_dir: Path, tid: str) -> Path:
    return Path(cache_dir) / f"{tid}.npz"


def seed_for(tid: str, stride: int = 1000) -> int:
    """Seed from the tree's own number, not its position in some cohort list, so a
    pilot run, a partial run and the full run all produce the same curve for a tree."""
    return stride * int(tid.rsplit("_", 1)[-1])


def load_tree_sweep(cache_dir, tid: str, p_values: Sequence[float], meta: Dict,
                    metric: str = "nmi") -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Return ``(<metric>_L, <metric>_B)`` for one tree if a matching cache exists."""
    path = _cache_file(Path(cache_dir), tid)
    if not path.exists():
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if list(np.round(z["p_values"], 6)) != list(np.round(p_values, 6)):
            return None
        stored = dict(z["meta"].item())
        if any(stored.get(k) != meta[k] for k in META_KEYS):
            return None
        return (np.asarray(z[f"{metric}_L"], float),
                np.asarray(z[f"{metric}_B"], float))
    except Exception:
        return None


def sweep_one_tree(tid: str, seed: int, p_values: Sequence[float], *, reps: int,
                   num_gaps: int, min_split: int,
                   cohort_name: str = "6000 taxa") -> Dict[str, np.ndarray]:
    """Run both arms on one tree. Returns ``{"<metric>_L"/"_B": one value per p}``."""
    from analysis.utils.real_cohorts import get_cohort
    from src.core.utils import compute_fiedler_from_laplacian, compute_laplacian
    from src.runners.p_sweep_inner import bootstrap_p_sweep_simple
    from src.utils.bpart_sweep_cache import bpart_sweep_raw

    loaded = get_cohort(cohort_name).load_all(tid)
    if loaded is None:
        raise FileNotFoundError(f"{tid}: alignment or tree missing")
    S, labels, tree, D = loaded
    gt = _gt_partition(tree, labels)

    fr_L = compute_fiedler_from_laplacian(compute_laplacian(S))
    out_L = bootstrap_p_sweep_simple(
        S, fr_L, list(p_values), bootstrap_reps=reps, seed=seed,
        num_gaps=num_gaps, min_split=min_split,
        partition_method="kmeans", laplacian="unnormalized")
    # the two arms name the same quantities differently; map both onto METRICS
    out = {f"{m}_L": np.asarray(out_L[k], float) for m, k in (
        ("nmi", "partition_nmi_M"), ("ari", "partition_ari_M"),
        ("agreement", "partition_agreement_M"), ("dot", "dot_product"),
        ("sign", "sign_agreement"))}

    raw = bpart_sweep_raw(D, list(p_values), reps=reps, seed_base=seed,
                          eigsolver="lm_k1", aggregation="avg_vector")
    for met in METRICS:
        out[f"{met}_B"] = np.array(
            [float(np.mean(pp[met])) for pp in raw["per_p"]], float)

    # second reference: the true tree's own split, scored from the partitions the two
    # sweeps just produced -- no extra eigensolve
    for arm, parts in (("L", out_L["partitions"]),
                       ("B", [pp["partitions"][0] for pp in raw["per_p"]])):
        scored = [_score(gt, part) for part in parts]
        for met in ("nmi", "ari", "agreement"):
            out[f"{met}_gt_{arm}"] = np.array([sc[met] for sc in scored], float)
    return out


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
        curves = sweep_one_tree(
            tid, seed=seed_for(tid, seed_stride), p_values=pv,
            reps=reps, num_gaps=num_gaps, min_split=min_split,
            cohort_name=cohort_name)
        np.savez(_cache_file(cache_dir, tid), p_values=pv,
                 meta=np.array(meta, dtype=object), **curves)
        finished.append(tid)
        el, per = time.time() - t0, (time.time() - t0) / k
        print(f"  [{k}/{len(todo)}] {tid}  {time.time() - t_tree:.0f} s  "
              f"(elapsed {el/3600:.1f} h, ETA {per*(len(todo)-k)/3600:.1f} h)", flush=True)
    return finished


def collect(ids: Sequence[str], cache_dir, p_values: Sequence[float], meta: Dict,
            metric: str = "nmi") -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    """Stack one metric's cached curves. Returns ``(arr_L, arr_B, have, pending)``."""
    if metric not in METRICS + GT_METRICS:
        raise ValueError(f"metric must be one of {METRICS + GT_METRICS}, got {metric!r}")
    pv = np.round(np.asarray(p_values, float), 6)
    have, pending, L, B = [], [], [], []
    for tid in ids:
        got = load_tree_sweep(cache_dir, tid, pv, meta, metric)
        if got is None:
            pending.append(tid)
            continue
        have.append(tid)
        L.append(got[0])
        B.append(got[1])
    arr_L = np.vstack(L) if L else np.empty((0, len(pv)))
    arr_B = np.vstack(B) if B else np.empty((0, len(pv)))
    return arr_L, arr_B, have, pending
