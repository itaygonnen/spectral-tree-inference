"""Sub-sampling recovery sweep on a real dataset, two operators.

For one tree the sweep asks: replace the full matrix by a uniform sub-sample at rate
``p`` and re-read the bipartition -- how much of the full-matrix partition survives?
Both arms are scored by NMI against their OWN full-matrix reference:

    L   Fiedler of Deg(S) - S, k-means OR sign cut -- whichever cuts the reference more
        evenly, decided per tree and recorded in the tree's own .npz
        (``bootstrap_p_sweep_simple``)
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

from src.utils.logging import create_progress_bar, log_info

# Cache-key fields. A tree whose stored meta differs is recomputed rather than mixed
# into a figure with a different grid. ``metrics`` is part of the key so caches written
# before every metric was stored are recomputed instead of read back half-empty.
META_KEYS = ("reps", "num_gaps", "min_split", "partition_method",
             "eigsolver", "aggregation", "m", "metrics")

# Metrics BOTH arms produce per p, so a figure can be drawn on any of them. NMI is what the
# recovery figure plots; the rest cost nothing once the eigenvector is in hand, so they are
# stored rather than thrown away. ``agreement`` and ``dot`` are the orientation-invariant
# clan-label match (%) and the dot product with the reference vector.
# Column names spell out WHAT is compared to WHAT: <metric>_<operator>_vs_<reference>.
# "gt" meant the true tree and nobody could tell.
METRICS = ("nmi", "ari", "agreement", "dot")
REF_FULL = "vs_fullmatrix"      # that operator's own split of the complete matrix
REF_TREE = "vs_truetree"        # the true tree's own top bipartition

# The same three scores against a SECOND reference: the true tree's own top bipartition
# (its root's two leaf sets), which is what the earlier real-data benchmark reported as
# ``*_gt``. Recovery towards the full matrix says the sub-sample kept what the full matrix
# saw; recovery towards the tree says the full matrix was seeing the right thing. They can
# disagree sharply -- 0.68 vs 0.013 on the m=1000 run -- so both are stored.
GT_METRICS = tuple(f"{m}_{REF_TREE}" for m in ("nmi", "ari", "agreement"))

# L-only extra: ``bootstrap_p_sweep_simple`` reports sign agreement separately, while on the
# B arm the sign pattern IS the partition, so its ``agreement`` already is that number.
EXTRA_L = ("sign",)


def sweep_meta(reps: int, num_gaps: int, min_split: int, m: int = 6000) -> Dict:
    return dict(reps=int(reps), num_gaps=int(num_gaps), min_split=int(min_split),
                partition_method="kmeans_or_sign", eigsolver="lm_k1",
                aggregation="avg_vector", m=int(m),
                metrics=",".join(tuple(f"{m}_{REF_FULL}" for m in METRICS)
                                 + GT_METRICS + ("split", "eta")))


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


def _eta_of(part: np.ndarray) -> float:
    """Imbalance of a boolean bipartition: larger clan / smaller clan."""
    n1 = int(np.sum(part))
    n2 = int(part.size) - n1
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


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
    """Seed from the tree's own number, not its position in some dataset list, so a
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
        # metric already carries its reference, e.g. "nmi_vs_fullmatrix"
        met, ref = metric.rsplit("_vs_", 1)
        return (np.asarray(z[f"{met}_L_vs_{ref}"], float),
                np.asarray(z[f"{met}_B_vs_{ref}"], float))
    except Exception:
        return None


def sweep_one_tree(tid: str, seed: int, p_values: Sequence[float], *, reps: int,
                   num_gaps: int, min_split: int,
                   dataset: str = "6000 taxa",
                   progress_cb=None) -> Dict[str, np.ndarray]:
    """Run both arms on one tree. Returns ``{"<metric>_L"/"_B": one value per p}``.

    ``progress_cb(stage, p_index, p)`` is called after every p of every arm, so a caller
    can show movement inside a tree that takes half an hour.
    """
    from analysis.utils.real_datasets import get_dataset
    from src.core.utils import compute_fiedler_from_laplacian, compute_laplacian
    from src.runners.p_sweep_inner import bootstrap_p_sweep_simple
    from src.utils.bpart_sweep_cache import bpart_sweep_raw

    loaded = get_dataset(dataset).load_all(tid)
    if loaded is None:
        raise FileNotFoundError(f"{tid}: alignment or tree missing")
    S, labels, tree, D = loaded
    gt = _gt_partition(tree, labels)
    log_info("bootstrap", f"{tid}: loaded m={S.shape[0]}, "
                          f"true-tree top split {int(gt.sum())}/{int((~gt).sum())}")

    def _cb(stage):
        if progress_cb is None:
            return None
        return lambda i, p: progress_cb(stage, i, p)

    fr_L = compute_fiedler_from_laplacian(compute_laplacian(S))
    # Both threshold rules on the reference vector; the more balanced one drives the
    # sweep. k-means alone routinely isolates one taxon on real data (eta ~ m), which has
    # no recovery signal at any p. Both etas are recorded so the choice is auditable.
    from src.runners.p_sweep_inner import _kmeans_bipartition
    cand = {"kmeans": np.asarray(_kmeans_bipartition(fr_L, min_split)).astype(bool),
            "sign": np.asarray(fr_L >= 0).astype(bool)}
    etas = {r: _eta_of(part) for r, part in cand.items()}
    rule_L = min(cand, key=lambda r: (etas[r], r != "kmeans"))
    log_info("bootstrap", f"{tid}: L(S) reference k-means eta={etas['kmeans']:.2f}, "
                          f"sign eta={etas['sign']:.2f} -> using {rule_L}")

    out_L = bootstrap_p_sweep_simple(
        S, fr_L, list(p_values), bootstrap_reps=reps, seed=seed,
        num_gaps=num_gaps, min_split=min_split,
        partition_method=rule_L, laplacian="unnormalized",
        progress_cb=_cb("L"))
    n1, n2 = out_L["partition_split_ref"]
    log_info("bootstrap", f"{tid}: L(S) reference split {n1}/{n2} ({rule_L})")
    # the two arms name the same quantities differently; map both onto METRICS
    out = {f"{m}_L_{REF_FULL}": np.asarray(out_L[k], float) for m, k in (
        ("nmi", "partition_nmi_M"), ("ari", "partition_ari_M"),
        ("agreement", "partition_agreement_M"), ("dot", "dot_product"),
        ("signagreement", "sign_agreement"))}

    raw = bpart_sweep_raw(D, list(p_values), reps=reps, seed_base=seed,
                          eigsolver="lm_k1", aggregation="avg_vector",
                          progress_cb=_cb("B"))
    log_info("bootstrap", f"{tid}: B reference split {raw['n1']}/{raw['n2']} "
                          f"(eta={raw['eta']:.2f})")
    for met in METRICS:
        out[f"{met}_B_{REF_FULL}"] = np.array(
            [float(np.mean(pp[met])) for pp in raw["per_p"]], float)

    # second reference: the true tree's own split, scored from the partitions the two
    # sweeps just produced -- no extra eigensolve
    for arm, parts in (("L", out_L["partitions"]),
                       ("B", [pp["partitions"][0] for pp in raw["per_p"]])):
        scored = [_score(gt, part) for part in parts]
        for met in ("nmi", "ari", "agreement"):
            out[f"{met}_{arm}_{REF_TREE}"] = np.array([sc[met] for sc in scored], float)
        # the split the sweep actually recovered at each p, not just its score
        out[f"split_small_{arm}"] = np.array(
            [np.nan if part is None else min(int(np.sum(part)),
                                             int(part.size - np.sum(part)))
             for part in parts], float)
        out[f"eta_{arm}"] = np.array(
            [np.nan if part is None else _eta_of(np.asarray(part))
             for part in parts], float)

    # provenance of the L arm's threshold choice, per tree
    out["eta_ref_L_kmeans"] = np.array([etas["kmeans"]], float)
    out["eta_ref_L_sign"] = np.array([etas["sign"]], float)
    out["eta_ref_B"] = np.array([raw["eta"]], float)
    out["rule_L"] = np.array([rule_L], dtype=object)
    return out


def run_sweep(ids: Sequence[str], cache_dir, p_values: Sequence[float], *,
              reps: int = 10, num_gaps: int = 10, min_split: int = 5,
              dataset: str = "6000 taxa", m: int = 6000,
              seed_stride: int = 1000) -> List[str]:
    """Sweep every id in ``ids``, skipping trees already cached. Returns ids done."""
    import time

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = sweep_meta(reps, num_gaps, min_split, m)
    pv = np.round(np.asarray(p_values, float), 6)

    finished: List[str] = []
    todo = [t for t in ids if load_tree_sweep(cache_dir, t, pv, meta) is None]
    print(f"  {len(ids) - len(todo)}/{len(ids)} trees cached, {len(todo)} to run "
          f"({len(pv)} p x {reps} reps, p={pv[0]:g}..{pv[-1]:g})", flush=True)
    log_info("bootstrap", f"{dataset}: {len(ids) - len(todo)}/{len(ids)} cached, "
                          f"{len(todo)} to run, p={pv[0]:g}..{pv[-1]:g}, reps={reps}")

    t0 = time.time()
    # ONE bar, not a nested pair: nested bars repaint over each other whenever the output
    # is piped to a log. The tree counter lives in the description, the ETA in the
    # postfix, and each bar covers the 2x|p| arm-steps of its own tree -- so it moves
    # several times a minute even at m=6000.
    for k, tid in enumerate(todo, 1):
        t_tree = time.time()
        inner = create_progress_bar(2 * len(pv),
                                    f"  [{k}/{len(todo)}] {tid}", unit="p", leave=False)
        curves = sweep_one_tree(
            tid, seed=seed_for(tid, seed_stride), p_values=pv,
            reps=reps, num_gaps=num_gaps, min_split=min_split,
            dataset=dataset,
            progress_cb=lambda stage, i, p: inner.update(1))
        inner.close()
        np.savez(_cache_file(cache_dir, tid), p_values=pv,
                 meta=np.array(meta, dtype=object), **curves)
        log_info("bootstrap", f"{tid}: L rule {curves['rule_L'][0]}, "
                              f"eta_ref L={curves['eta_ref_L_' + curves['rule_L'][0]][0]:.2f} "
                              f"B={curves['eta_ref_B'][0]:.2f}")
        finished.append(tid)
        el = time.time() - t0
        eta_h = (el / k) * (len(todo) - k) / 3600.0
        log_info("bootstrap",
                 f"[{k}/{len(todo)}] {tid} done in {time.time() - t_tree:.0f} s "
                 f"(elapsed {el/3600:.2f} h, ETA {eta_h:.2f} h)")
    print(f"  {len(finished)} tree(s) swept in {(time.time() - t0)/3600:.2f} h",
          flush=True)
    return finished


def collect(ids: Sequence[str], cache_dir, p_values: Sequence[float], meta: Dict,
            metric: str = "nmi") -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    """Stack one metric's cached curves. Returns ``(arr_L, arr_B, have, pending)``."""
    known = tuple(f"{m}_{REF_FULL}" for m in METRICS) + GT_METRICS
    if metric not in known:
        raise ValueError(f"metric must be one of {known}, got {metric!r}")
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
