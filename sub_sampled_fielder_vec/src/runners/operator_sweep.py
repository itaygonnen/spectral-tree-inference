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

import re
import sys
import time

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
# to gain a column). ``early_stop`` is not here either, for a different reason: it is
# stored and compared, but under a partial order rather than for equality -- see
# ``load_tree_sweep``.
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
            + tuple(f"{x}_{a}" for x in ("eta", "split_small", "measured") for a in arms)
            + tuple(f"signagreement_{a}_{REF_FULL}" for a in fiedler)
            + tuple(f"{name}_{a}" for name in DIAGNOSTICS.values() for a in arms))


COLUMNS = columns()


def sweep_meta(reps: int, num_gaps: int, min_split: int, m: int = 6000,
               early_stop: int = 0) -> Dict:
    """The cache key for one tree's sweep. Unchanged since the first cluster run.

    ``early_stop`` is stored but NOT in :data:`META_KEYS` -- it is compared by
    :func:`load_tree_sweep` under a partial order rather than for equality, because a
    tree swept with a smaller (or no) early stop holds strictly more measurement than
    one swept with a larger one.
    """
    return dict(reps=int(reps), num_gaps=int(num_gaps), min_split=int(min_split),
                partition_method="kmeans_or_sign", eigsolver="lm_k1",
                aggregation="avg_vector", m=int(m), early_stop=int(early_stop),
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
        # Asymmetric on the early stop, unlike every key above. A tree swept with a
        # SMALLER stop measured at least every p this run wants measured, so it is
        # strictly better data and is reused; 0 (measured throughout) serves any run.
        # The converse is refused: a run that asked for measurement is never handed an
        # extrapolated tail.
        want = int(meta.get("early_stop", 0) or 0)
        have = int(stored.get("early_stop", 0) or 0)
        if have and (not want or have > want):
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
                   aggregation: str = "avg_vector", early_stop: int = 0,
                   progress_cb=None) -> Dict[str, np.ndarray]:
    """Every selected operator on one tree. Returns ``{column: one value per p}``.

    ``progress_cb(arm, p_index, p, nmi)`` is called after every p of every arm, so a caller
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

    def _log_p(arm: str, idx: int, p: float, row: Dict) -> None:
        """One line per (tree, arm, p), the moment that p finishes.

        The pre-refactor pipeline logged this and the refactor dropped it, so a tree that
        takes half an hour at m=6000 said nothing between "loaded" and "done". force=True
        because it must reach experiment.log whether or not the terminal is drawing bars.

        A p the early stop extrapolated is logged too, and says so: reading the log, the
        gap between two measured p's would otherwise look like a crash.
        """
        part = row.get("partition")
        if part is None:
            # an extrapolated low-side row carries no partition at all, not even a None
            parts = row.get("partitions") or []
            part = parts[0] if len(parts) else None
        if part is not None:
            part = np.asarray(part)
            a = int(part.sum())
            split = f"{min(a, part.size - a)}-{max(a, part.size - a)}"
            eta = f"{_eta_of(part):.2f}"
        else:
            split = eta = "N/A"

        def _num(key):
            v = row.get(key)
            if v is None:
                return None
            v = float(np.mean(v)) if isinstance(v, (list, tuple, np.ndarray)) else float(v)
            return v

        bits = []
        for key, label in (("nmi", "nmi"), ("ari", "ari"), ("agreement", "agr"),
                           ("dot", "dot"), ("sign_agreement", "sign")):
            v = _num(key)
            if v is not None:
                bits.append(f"{label}={v:.4f}" if label in ("nmi", "ari", "dot")
                            else f"{label}={v:.2f}%")
        for key, label in (("sigma2_avg_M", "s2_full"), ("sigma2_avg_S", "s2_sub"),
                           ("mean_operator_norm_error", "opnorm")):
            v = _num(key)
            if v is not None:
                bits.append(f"{label}={v:.4f}")
        kind = "results" if row.get("measured", 1.0) else "EXTRAPOLATED"
        # NOT force=True: that flag is for important final messages, and this is one
        # line per (tree, p) -- 160 of them on an 8-tree run at 20 p, which buried the
        # progress display it now sits behind. It still reaches run.log either way, and
        # --display-mode debug still prints it.
        log_info("bootstrap",
                 f"{tree_id} [{arm}] p={p:.4g} {kind}: " + ", ".join(bits)
                 + f", eta={eta}, split={split}")

    out: Dict[str, np.ndarray] = {}
    n_p = len(list(p_values))
    for key in resolve(operators):
        op = OPERATORS[key]
        arm = op.arm
        cb = (None if progress_cb is None
              else (lambda i, p, nmi, _a=arm: progress_cb(_a, i, p, nmi)))
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
                laplacian=op.laplacian, extra_metrics=extra_metrics,
                early_stop=early_stop, progress_cb=cb,
                on_p_result=lambda i, pp, row, _a=arm: _log_p(_a, i, pp, row))
            for metric, src in (("nmi", "partition_nmi_M"), ("ari", "partition_ari_M"),
                                ("agreement", "partition_agreement_M"),
                                ("dot", "dot_product")):
                out[f"{metric}_{arm}_{REF_FULL}"] = np.asarray(res[src], float)
            out[f"signagreement_{arm}_{REF_FULL}"] = np.asarray(res["sign_agreement"],
                                                                float)
            out[f"measured_{arm}"] = np.asarray(res["measured"], float)
            # one partition per p: bootstrap_p_sweep_simple averages the vectors first
            parts = [[q] if q is not None else [] for q in res["partitions"]]
            diagnostics = {k: np.asarray(res[k], float) for k in DIAGNOSTICS
                           if k in res}
        else:
            raw = bpart_sweep_raw(D, list(p_values), reps=reps, seed_base=seed,
                                  eigsolver=eigsolver, aggregation=aggregation,
                                  extra_metrics=extra_metrics, early_stop=early_stop,
                                  progress_cb=cb,
                                  on_p_result=lambda i, pp, row, _a=arm:
                                  _log_p(_a, i, pp, row))
            for metric in METRICS:
                out[f"{metric}_{arm}_{REF_FULL}"] = np.array(
                    [float(np.mean(pp[metric])) for pp in raw["per_p"]], float)
            # under avg_vector each p holds one partition; under per_rep it holds one
            # per replicate, and every column below averages over them -- the same
            # accounting the metrics above already use. Taking only the first would have
            # put two different accountings in one row.
            parts = [list(pp["partitions"]) for pp in raw["per_p"]]
            out[f"measured_{arm}"] = np.array(
                [float(pp.get("measured", 1.0)) for pp in raw["per_p"]], float)
            diagnostics = {k: np.array([float(pp[k]) for pp in raw["per_p"]], float)
                           for k in DIAGNOSTICS if k in raw["per_p"][0]}

        # the split the sweep actually recovered at each p, not just its score, averaged
        # over whatever partitions that p produced (one, or one per replicate)
        def _mean(vals):
            vals = [v for v in vals if v is not None and not np.isnan(v)]
            return float(np.mean(vals)) if vals else float("nan")

        out[f"split_small_{arm}"] = np.array(
            [_mean([min(int(np.sum(q)), int(q.size - np.sum(q)))
                    for q in ps if q is not None]) for ps in parts], float)
        out[f"eta_{arm}"] = np.array(
            [_mean([_eta_of(np.asarray(q)) for q in ps if q is not None])
             for ps in parts], float)
        # second reference: the true tree's own split, scored from the partitions the
        # sweep just produced -- no extra eigensolve
        for metric in TREE_METRICS:
            out[f"{metric}_{arm}_{REF_TREE}"] = np.array(
                [_mean([_score(gt, q)[metric] for q in ps if q is not None])
                 if gt is not None else float("nan") for ps in parts], float)
        for src, name in DIAGNOSTICS.items():
            if src in diagnostics:
                out[f"{name}_{arm}"] = diagnostics[src]
        assert all(v.size in (1, n_p) for k, v in out.items()
                   if isinstance(v, np.ndarray) and v.dtype.kind == "f"), \
            f"{tree_id}: a column came back with the wrong length"
    return out


# The arm the progress shades. DEFAULT_OPERATORS is ("S",) whose arm is "L"; a run
# sweeping several arms still draws one grid, because the figure reads one curve.
REF_ARM = "L"


class SweepProgress:
    """The p-grid filling in, one row per tree, NMI shaded into the cell.

    The sibling of ``sailfish_pool.Progress``, and for the same reason: the sweep used to
    emit one INFO line per (tree, p), so a modest 8-tree run at 20 p scrolled 160 lines of
    numbers nobody reads while it ran. The recovery curve is the point of the sweep, so
    the progress IS the curve: cells go dark as NMI rises, and the transition appears from
    the left as it is measured.

    Two modes, as elsewhere -- repaint on a TTY, one line per finished tree in a log,
    because a repainting block in a cluster log is noise.
    """

    # NMI -> glyph. The ramp is deliberately coarse: it reads as a curve, not a table.
    RAMP = ((0.10, "\u00b7"), (0.40, "\u2591"), (0.70, "\u2592"),
            (0.95, "\u2593"), (1.01, "\u2588"))
    PENDING = " "

    def __init__(self, tree_ids, p_values, threshold: float = 0.90,
                 tty: Optional[bool] = None):
        self.ids = list(tree_ids)
        self.p = list(p_values)
        self.threshold = threshold
        self.cells: Dict[str, List[str]] = {t: [self.PENDING] * len(self.p)
                                            for t in self.ids}
        self.nmi: Dict[str, List[float]] = {t: [float("nan")] * len(self.p)
                                            for t in self.ids}
        self.done: List[str] = []
        self.current: Optional[str] = None
        self.t0 = time.time()
        # getattr, not sys.stdout.isatty(): a run tees stdout into run.log through a
        # wrapper with no isatty, and treating that as "not a terminal" is also the right
        # answer -- repainting would write escape codes into the log file.
        self.tty = (getattr(sys.stdout, "isatty", lambda: False)()
                    if tty is None else tty)
        self._drawn = 0

    @classmethod
    def _glyph(cls, nmi: float) -> str:
        if nmi != nmi:                       # NaN
            return cls.PENDING
        for bound, ch in cls.RAMP:
            if nmi < bound:
                return ch
        return cls.RAMP[-1][1]

    def p_star(self, tree_id: str) -> Optional[float]:
        """Smallest measured p whose NMI clears the threshold -- the run's read-off."""
        for value, nmi in zip(self.p, self.nmi[tree_id]):
            if nmi == nmi and nmi >= self.threshold:
                return value
        return None

    @staticmethod
    def _bucket_of(tree_id: str) -> Optional[int]:
        """The eta bucket from the stem ``eta05_003``, or None for an unbucketed set.

        Read off the name rather than plumbed through, because ``run_sweep`` is given ids
        and a loader and deliberately knows nothing about where they came from -- a
        downloaded set named ``random_tree_001`` simply groups flat.
        """
        m = re.match(r"^eta(\d+)_", str(tree_id))
        return int(m.group(1)) if m else None

    def _groups(self) -> List[Tuple[Optional[int], List[str]]]:
        """``[(bucket, tree ids)]`` in bucket order; one flat group when unbucketed."""
        by: Dict[Optional[int], List[str]] = {}
        for t in self.ids:
            by.setdefault(self._bucket_of(t), []).append(t)
        if list(by) == [None]:
            return [(None, self.ids)]
        return sorted(by.items(), key=lambda kv: (kv[0] is None, kv[0] or 0))

    def _lines(self) -> List[str]:
        width = max((len(t) for t in self.ids), default=4) + 4
        head = (f"  {'tree':<{width}}  p {self.p[0]:.0e} "
                + "-" * max(len(self.p) - 12, 1) + f"> {self.p[-1]:g}"
                + f"   p*>={self.threshold:g}")
        out = [head]
        for bucket, members in self._groups():
            stars = [self.p_star(t) for t in members]
            got = [v for v in stars if v is not None]
            if bucket is not None:
                out.append(f"  eta{bucket:02d}  {len(members)} tree(s)"
                           + (f"   median p* {float(np.median(got)):.4g}"
                              if len(got) == len(members) else
                              f"   {len(got)}/{len(members)} past threshold"))
            for t, star in zip(members, stars):
                mark = "<" if t == self.current else " "
                indent = "    " if bucket is not None else "  "
                out.append(f"{indent}{t:<{width - len(indent) + 2}}{mark}"
                           + "".join(self.cells[t])
                           + "  " + (f"{star:.4g}" if star is not None else "-"))
        el = time.time() - self.t0
        eta = (el / len(self.done) * (len(self.ids) - len(self.done))
               if self.done else 0.0)
        out.append(f"  {len(self.done)}/{len(self.ids)} trees, {_clock(el)} elapsed"
                   + (f", ETA {_clock(eta)}" if self.done and
                      len(self.done) < len(self.ids) else ""))
        return out

    def _redraw(self) -> None:
        if self._drawn:
            sys.stdout.write(f"\033[{self._drawn}A")
        lines = self._lines()
        for line in lines:
            sys.stdout.write("\033[2K" + line + "\n")
        sys.stdout.flush()
        self._drawn = len(lines)

    def start(self, tree_id: str) -> None:
        self.current = tree_id
        if self.tty:
            self._redraw()

    def point(self, tree_id: str, index: int, nmi: float) -> None:
        if 0 <= index < len(self.p):
            self.nmi[tree_id][index] = nmi
            self.cells[tree_id][index] = self._glyph(nmi)
        if self.tty:
            self._redraw()

    def finish(self, tree_id: str, secs: float) -> None:
        self.done.append(tree_id)
        self.current = None
        if self.tty:
            self._redraw()
        else:
            star = self.p_star(tree_id)
            print(f"  [{len(self.done)}/{len(self.ids)}] {tree_id} "
                  f"{secs:.0f}s  p*={star:.4g}" if star is not None
                  else f"  [{len(self.done)}/{len(self.ids)}] {tree_id} "
                       f"{secs:.0f}s  p*=-", flush=True)

    def close(self) -> None:
        if self.tty:
            self._redraw()
        print(f"  shading: {chr(183)} NMI<0.1  {chr(9617)} <0.4  {chr(9618)} <0.7  "
              f"{chr(9619)} <0.95  {chr(9608)} ~1.0", flush=True)


def _clock(seconds: float) -> str:
    m, s = divmod(int(max(seconds, 0)), 60)
    h, m = divmod(m, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"

def run_sweep(ids: Sequence[str], cache_dir, loader, p_values: Sequence[float], *,
              reps: int = 10, num_gaps: int = 10, min_split: int = MIN_SPLIT,
              operators: Sequence[str] = (), source: str = "", m: int = 6000,
              extra_metrics: bool = False, rule_policy: str = "best",
              eigsolver: str = "lm_k1", aggregation: str = "avg_vector",
              early_stop: int = 0, seed_stride: int = 1000) -> List[str]:
    """Sweep every id in ``ids``, skipping trees already cached. Returns ids done."""
    import time

    ops = resolve(operators)
    arms = [ARM_OF[k] for k in ops]
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    meta = dict(sweep_meta(reps, num_gaps, min_split, m, early_stop),
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
    bar = SweepProgress(todo, pv)
    bar.start(todo[0]) if todo else None
    for k, tree_id in enumerate(todo, 1):
        t_tree = time.time()
        bar.start(tree_id)
        curves = sweep_one_tree(
            tree_id, loader, seed=seed_for(tree_id, seed_stride), p_values=pv,
            reps=reps, num_gaps=num_gaps, min_split=min_split, operators=ops,
            extra_metrics=extra_metrics, rule_policy=rule_policy,
            eigsolver=eigsolver, aggregation=aggregation, early_stop=early_stop,
            progress_cb=lambda arm, i, p, nmi, _t=tree_id: (
                bar.point(_t, i, nmi) if arm == REF_ARM else None))
        np.savez(_cache_file(cache_dir, tree_id), p_values=pv,
                 meta=np.array(meta, dtype=object), **curves)
        finished.append(tree_id)
        el = time.time() - t0
        bar.finish(tree_id, time.time() - t_tree)
        log_info("bootstrap",
                 f"[{k}/{len(todo)}] {tree_id} done in {time.time() - t_tree:.0f} s "
                 f"(elapsed {el/3600:.2f} h, ETA {(el/k)*(len(todo)-k)/3600:.2f} h)")
    if todo:
        bar.close()
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
