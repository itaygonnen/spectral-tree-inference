"""Sub-sampling recovery sweep for any operator, on a tree from any source.

For one tree the sweep asks: replace the full matrix by a uniform sub-sample at rate
``p`` and re-read the bipartition -- how much of the full-matrix partition survives?
Every selected operator is scored against its OWN full-matrix reference and against
the true tree's top bipartition, so the two questions stay separable:

    vs_fullmatrix   did sub-sampling keep what the complete matrix saw?
    vs_truetree     was the complete matrix seeing the right thing?

Which operators exist, and how each one is cut, is :mod:`src.runners.operators`.
Where the tree comes from is the ``loader``'s business -- this module never touches a
FASTA file or a simulator, which is why one sweep now serves both.

Two settings are not the library defaults and both are load-bearing:

* ``eigsolver="lm_k1"`` for B -- ARPACK for the single eigenpair it needs. Same vector
  as the dense default at 0.39 s vs 17.2 s per solve at m=6000.
* ``aggregation="avg_vector"`` -- bootstrap-average the eigenvector, THEN partition
  once, which is what ``bootstrap_p_sweep_simple`` does on the Fiedler arms. The
  ``per_rep`` default would score B under a different estimator.

One ``.npz`` per tree, so a run is resumable and interruptible at any tree boundary.
Columns are **additive**: they are not part of the cache key, so a tree swept before an
operator or a diagnostic existed still counts as swept and simply lacks those columns.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..utils.logging import create_progress_bar, log_info
from ..utils.partition_metrics import eta as _eta_of
from ..utils.partition_metrics import score as _score
from .operators import ARM_OF, MIN_SPLIT, OPERATORS, choose_cut, resolve
from .p_sweep_inner import EXTRA_METRIC_KEYS

# Cache-key fields: a tree whose stored meta differs on one of these is recomputed
# rather than mixed into a figure with a different grid. Neither the column list nor
# the operator list is here -- both are additive, and putting them in the key would
# make every cached tree a miss the moment one was added (a 50 h recompute at m=6000
# to gain a column).
META_KEYS = ("reps", "num_gaps", "min_split", "partition_method",
             "eigsolver", "aggregation", "m")

METRICS = ("nmi", "ari", "agreement", "dot")
REF_FULL = "vs_fullmatrix"      # that operator's own split of the complete matrix
REF_TREE = "vs_truetree"        # the true tree's own top bipartition
# ``dot`` needs the reference eigenvector, which the true tree does not have
TREE_METRICS = ("nmi", "ari", "agreement")

# The linear-algebra diagnostics the original pipeline-A runs recorded. The sweep keeps
# the quantities but not the names: "empirical_rank_L_S" says nothing about which arm or
# which matrix once three operators share one table. Here "full" is the complete matrix
# an arm reads (S on the Fiedler arms, D on B), "sub" is its bootstrap-averaged
# sub-sample, and "op" is the operator built from it.
DIAGNOSTICS = {
    "sigma2_avg_M":               "sigma2_full",
    "sigma2_avg_S":               "sigma2_sub",
    "mean_operator_norm_error":   "opnorm_err_mean",
    "median_operator_norm_error": "opnorm_err_median",
    "std_operator_norm_error":    "opnorm_err_std",
    "mean_empirical_rank_S":      "rank_sub_mean",
    "median_empirical_rank_S":    "rank_sub_median",
    "std_empirical_rank_S":       "rank_sub_std",
    "mean_empirical_rank_L_S":    "rank_op_sub_mean",
    "median_empirical_rank_L_S":  "rank_op_sub_median",
    "std_empirical_rank_L_S":     "rank_op_sub_std",
    "empirical_rank_M":           "rank_full",
    "empirical_rank_L_M":         "rank_op_full",
}
# a key added to the sweep but never mapped to a column name would be collected and
# then silently dropped on the way to the CSV, so fail here instead
assert tuple(DIAGNOSTICS) == tuple(EXTRA_METRIC_KEYS), (
    "DIAGNOSTICS is out of step with p_sweep_inner.EXTRA_METRIC_KEYS: "
    f"{set(EXTRA_METRIC_KEYS) ^ set(DIAGNOSTICS)}")


def columns(operators: Sequence[str] = ()) -> Tuple[str, ...]:
    """Every per-p column a tree's ``.npz`` can hold, in the order the CSVs print them."""
    arms = [ARM_OF[k] for k in resolve(operators)]
    fiedler = [ARM_OF[k] for k in resolve(operators) if OPERATORS[k].kind == "fiedler"]
    return (tuple(f"{m}_{a}_{REF_FULL}" for m in METRICS for a in arms)
            + tuple(f"{m}_{a}_{REF_TREE}" for m in TREE_METRICS for a in arms)
            + tuple(f"{x}_{a}" for x in ("eta", "split_small") for a in arms)
            + tuple(f"signagreement_{a}_{REF_FULL}" for a in fiedler)
            + tuple(f"{name}_{a}" for name in DIAGNOSTICS.values() for a in arms))


COLUMNS = columns()


def sweep_meta(reps: int, num_gaps: int, min_split: int, m: int = 6000) -> Dict:
    """The cache key for one tree's sweep. Unchanged since the first cluster run."""
    return dict(reps=int(reps), num_gaps=int(num_gaps), min_split=int(min_split),
                partition_method="kmeans_or_sign", eigsolver="lm_k1",
                aggregation="avg_vector", m=int(m),
                metrics=",".join(COLUMNS))


def tree_top_bipartition(tree, labels) -> Optional[np.ndarray]:
    """The true tree's top bipartition, in ``labels`` order, or None without a tree.

    The root's two child subtrees, matching ``src.runners.real_data_bpart`` so the
    numbers are comparable with the earlier benchmark. Falls back to the first internal
    node with two children when the seed node is degenerate.
    """
    if tree is None:
        return None
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


def _cache_file(cache_dir: Path, tree_id: str) -> Path:
    return Path(cache_dir) / f"{tree_id}.npz"


def seed_for(tree_id: str, stride: int = 1000) -> int:
    """Seed from the tree's own number, not its position in some list, so a pilot run,
    a partial run and the full run all produce the same curve for a tree."""
    tail = tree_id.rsplit("_", 1)[-1]
    return stride * (int(tail) if tail.isdigit() else abs(hash(tree_id)) % 100000)


def load_tree_sweep(cache_dir, tree_id: str, p_values: Sequence[float], meta: Dict,
                    metric: str = f"nmi_{REF_FULL}", arms: Sequence[str] = ("L", "B")
                    ) -> Optional[Tuple[np.ndarray, ...]]:
    """One metric's curve per arm for one tree, if a matching cache exists.

    ``metric`` names its reference, e.g. ``"nmi_vs_fullmatrix"``; a bare ``"nmi"`` is
    read as the full-matrix reference. A column the cached tree predates comes back as
    NaN rather than a cache miss -- the tree IS swept, it just does not carry that
    column, and recomputing it costs ~30 min at m=6000.
    """
    if "_vs_" not in metric:
        metric = f"{metric}_{REF_FULL}"
    path = _cache_file(Path(cache_dir), tree_id)
    if not path.exists():
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if list(np.round(z["p_values"], 6)) != list(np.round(p_values, 6)):
            return None
        stored = dict(z["meta"].item())
        if any(stored.get(k) != meta[k] for k in META_KEYS):
            return None
        met, ref = metric.rsplit("_vs_", 1)
        nan = np.full(len(z["p_values"]), np.nan)
        return tuple(np.asarray(z[f"{met}_{a}_vs_{ref}"], float)
                     if f"{met}_{a}_vs_{ref}" in z.files else nan for a in arms)
    except Exception:
        return None


def sweep_one_tree(tree_id: str, loader, seed: int, p_values: Sequence[float], *,
                   reps: int, num_gaps: int = 10, min_split: int = MIN_SPLIT,
                   operators: Sequence[str] = (), extra_metrics: bool = False,
                   rule_policy: str = "best", eigsolver: str = "lm_k1",
                   aggregation: str = "avg_vector",
                   progress_cb=None) -> Dict[str, np.ndarray]:
    """Every selected operator on one tree. Returns ``{column: one value per p}``.

    ``progress_cb(arm, p_index, p)`` is called after every p of every arm, so a caller
    can show movement inside a tree that takes half an hour.
    """
    from ..utils.bpart_sweep_cache import bpart_sweep_raw
    from .p_sweep_inner import bootstrap_p_sweep_simple

    loaded = loader(tree_id)
    if loaded is None:
        raise FileNotFoundError(f"{tree_id}: alignment or tree missing")
    S, labels, tree, D = loaded
    gt = tree_top_bipartition(tree, labels)
    log_info("bootstrap", f"{tree_id}: loaded m={S.shape[0]}" + (
        f", true-tree top split {int(gt.sum())}/{int((~gt).sum())}" if gt is not None
        else ", no true tree -- vs_truetree columns will be NaN"))

    out: Dict[str, np.ndarray] = {}
    n_p = len(list(p_values))
    for key in resolve(operators):
        op = OPERATORS[key]
        arm = op.arm
        cb = (None if progress_cb is None
              else (lambda i, p, _a=arm: progress_cb(_a, i, p)))
        vec = op.vector(S, D)
        rule, _part, etas = choose_cut(op, vec, min_split, S, num_gaps,
                                       rule_policy)
        for r, e in etas.items():
            out[f"eta_ref_{arm}_{r}"] = np.array([e], float)
        out[f"eta_ref_{arm}"] = np.array([etas[rule]], float)
        out[f"rule_{arm}"] = np.array([rule], dtype=object)
        log_info("bootstrap", f"{tree_id}: {key} reference "
                              + ", ".join(f"{r} eta={e:.2f}" for r, e in etas.items())
                              + f" -> using {rule}")

        if op.kind == "fiedler":
            res = bootstrap_p_sweep_simple(
                S, vec, list(p_values), bootstrap_reps=reps, seed=seed,
                num_gaps=num_gaps, min_split=min_split, partition_method=rule,
                laplacian=op.laplacian, extra_metrics=extra_metrics, progress_cb=cb)
            for metric, src in (("nmi", "partition_nmi_M"), ("ari", "partition_ari_M"),
                                ("agreement", "partition_agreement_M"),
                                ("dot", "dot_product")):
                out[f"{metric}_{arm}_{REF_FULL}"] = np.asarray(res[src], float)
            out[f"signagreement_{arm}_{REF_FULL}"] = np.asarray(res["sign_agreement"],
                                                                float)
            parts = res["partitions"]
            diagnostics = {k: np.asarray(res[k], float) for k in DIAGNOSTICS
                           if k in res}
        else:
            raw = bpart_sweep_raw(D, list(p_values), reps=reps, seed_base=seed,
                                  eigsolver=eigsolver, aggregation=aggregation,
                                  extra_metrics=extra_metrics, progress_cb=cb)
            for metric in METRICS:
                out[f"{metric}_{arm}_{REF_FULL}"] = np.array(
                    [float(np.mean(pp[metric])) for pp in raw["per_p"]], float)
            # under per_rep aggregation each p holds one partition per replicate; the
            # recovered-split columns describe the first, which is the only one a
            # single column can describe
            parts = [pp["partitions"][0] for pp in raw["per_p"]]
            diagnostics = {k: np.array([float(pp[k]) for pp in raw["per_p"]], float)
                           for k in DIAGNOSTICS if k in raw["per_p"][0]}

        # the split the sweep actually recovered at each p, not just its score
        out[f"split_small_{arm}"] = np.array(
            [np.nan if p is None else min(int(np.sum(p)), int(p.size - np.sum(p)))
             for p in parts], float)
        out[f"eta_{arm}"] = np.array(
            [np.nan if p is None else _eta_of(np.asarray(p)) for p in parts], float)
        # second reference: the true tree's own split, scored from the partitions the
        # sweep just produced -- no extra eigensolve
        scored = [_score(gt, p) if gt is not None else _score(None, None)
                  for p in parts]
        for metric in TREE_METRICS:
            out[f"{metric}_{arm}_{REF_TREE}"] = np.array(
                [sc[metric] for sc in scored], float)
        for src, name in DIAGNOSTICS.items():
            if src in diagnostics:
                out[f"{name}_{arm}"] = diagnostics[src]
        assert all(v.size in (1, n_p) for k, v in out.items()
                   if isinstance(v, np.ndarray) and v.dtype.kind == "f"), \
            f"{tree_id}: a column came back with the wrong length"
    return out


def run_sweep(ids: Sequence[str], cache_dir, loader, p_values: Sequence[float], *,
              reps: int = 10, num_gaps: int = 10, min_split: int = MIN_SPLIT,
              operators: Sequence[str] = (), source: str = "", m: int = 6000,
              extra_metrics: bool = False, rule_policy: str = "best",
              eigsolver: str = "lm_k1", aggregation: str = "avg_vector",
              seed_stride: int = 1000) -> List[str]:
    """Sweep every id in ``ids``, skipping trees already cached. Returns ids done."""
    import time

    ops = resolve(operators)
    arms = [ARM_OF[k] for k in ops]
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = dict(sweep_meta(reps, num_gaps, min_split, m),
                partition_method=("kmeans_or_sign" if rule_policy == "best"
                                  else rule_policy),
                eigsolver=eigsolver, aggregation=aggregation)
    pv = np.round(np.asarray(p_values, float), 6)

    todo = [t for t in ids if load_tree_sweep(cache_dir, t, pv, meta) is None]
    print(f"  {len(ids) - len(todo)}/{len(ids)} trees cached, {len(todo)} to run "
          f"({len(pv)} p x {reps} reps x {len(ops)} operators, "
          f"p={pv[0]:g}..{pv[-1]:g})", flush=True)
    log_info("bootstrap", f"{source}: {len(ids) - len(todo)}/{len(ids)} cached, "
                          f"{len(todo)} to run, operators {ops}, "
                          f"p={pv[0]:g}..{pv[-1]:g}, reps={reps}")

    finished: List[str] = []
    t0 = time.time()
    # ONE bar, not a nested pair: nested bars repaint over each other whenever the
    # output is piped to a log. The tree counter lives in the description, the ETA in
    # the postfix, and each bar covers the |arms| x |p| steps of its own tree -- so it
    # moves several times a minute even at m=6000.
    for k, tree_id in enumerate(todo, 1):
        t_tree = time.time()
        inner = create_progress_bar(len(arms) * len(pv),
                                    f"  [{k}/{len(todo)}] {tree_id}", unit="p",
                                    leave=False)
        curves = sweep_one_tree(
            tree_id, loader, seed=seed_for(tree_id, seed_stride), p_values=pv,
            reps=reps, num_gaps=num_gaps, min_split=min_split, operators=ops,
            extra_metrics=extra_metrics, rule_policy=rule_policy,
            eigsolver=eigsolver, aggregation=aggregation,
            progress_cb=lambda arm, i, p: inner.update(1))
        inner.close()
        np.savez(_cache_file(cache_dir, tree_id), p_values=pv,
                 meta=np.array(meta, dtype=object), **curves)
        finished.append(tree_id)
        el = time.time() - t0
        log_info("bootstrap",
                 f"[{k}/{len(todo)}] {tree_id} done in {time.time() - t_tree:.0f} s "
                 f"(elapsed {el/3600:.2f} h, ETA {(el/k)*(len(todo)-k)/3600:.2f} h)")
    print(f"  {len(finished)} tree(s) swept in {(time.time() - t0)/3600:.2f} h",
          flush=True)
    return finished


def collect(ids: Sequence[str], cache_dir, p_values: Sequence[float], meta: Dict,
            metric: str = f"nmi_{REF_FULL}", arms: Sequence[str] = ("L", "B")):
    """Stack one metric's cached curves. Returns ``(arrays_per_arm, have, pending)``."""
    pv = np.round(np.asarray(p_values, float), 6)
    have, pending = [], []
    stacks: List[List[np.ndarray]] = [[] for _ in arms]
    for tree_id in ids:
        got = load_tree_sweep(cache_dir, tree_id, pv, meta, metric, arms)
        if got is None:
            pending.append(tree_id)
            continue
        have.append(tree_id)
        for i, curve in enumerate(got):
            stacks[i].append(curve)
    out = tuple(np.vstack(s) if s else np.empty((0, len(pv))) for s in stacks)
    return out, have, pending
