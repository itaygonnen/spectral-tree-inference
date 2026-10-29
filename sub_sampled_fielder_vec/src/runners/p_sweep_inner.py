"""Reusable bootstrap p-sweep over a precomputed (M, fiedler_ref) pair.

Narrow extraction of the inner loop in ``bootstrap_sweep.py::sweep_for_params``,
focused on ``partition_agreement_M``. Drops the metric_composer, persistent
caching, guardrails, and partition-validity checks that the full runner needs
but a notebook does not.

Use ``sweep_for_params`` for full experiments; use this helper when you
already have ``(M, v_pop)`` (e.g. from the eta-pool cache) and want a
quick agreement curve.
"""
from __future__ import annotations

import warnings
from typing import Any, Callable, Dict, List, Optional

import numpy as np

from ..core.utils import compute_laplacian, compute_normalized_laplacian
from ..utils.logging import log_warning
from ..utils.metrics import (
    _normalize_vector,
    compute_fiedler_dot_product,
    compute_reference_partition_and_quality,
    compute_sign_agreement,
    numerical_rank,
    operator_norm_error,
    sigma2_of_partition,
)
from ..utils.random_entries import _subsample_matrix_entries, compute_fiedler_from_laplacian


# The linear-algebra diagnostics the original pipeline-A ``results.json`` carried beside
# the partition metrics. Collected only when ``extra_metrics=True``; every one is O(n^2)
# per replicate (Lanczos / power iteration / randomised SVD), never a dense factorisation.
EXTRA_METRIC_KEYS = (
    "sigma2_avg_M", "sigma2_avg_S",
    "mean_operator_norm_error", "median_operator_norm_error", "std_operator_norm_error",
    "mean_empirical_rank_S", "median_empirical_rank_S", "std_empirical_rank_S",
    "mean_empirical_rank_L_S", "median_empirical_rank_L_S", "std_empirical_rank_L_S",
    "empirical_rank_M", "empirical_rank_L_M",
)

# --- middle-out traversal ----------------------------------------------------------
# The grid is swept from its middle outwards, one step per side, so a sweep can stop on
# either end the moment the curve has said what it is going to say:
#
#   high side   NMI == 1.0 twice running        -> every larger p recovers perfectly too
#   low  side   NMI <= LOW_FLOOR twice running  -> every smaller p sits on the floor
#
# Both inner sweeps seed their RNG from the p INDEX and never from visit order
# (``seed + i`` here, ``seed_base + 10_000*p_idx + rep`` in ``bpart_sweep_cache``), so
# visiting the grid middle-out is bit-exact on every point still computed.
#
# Serial and single-step on purpose. ``middle_out_runner`` expands in parallel blocks and
# widens ``workers_per_side`` once one side finishes; ``offset`` then advances past where
# the other side stopped and the indices in between are never computed AND never filled
# (n_p=20, workers=8, high guardrail in round 1 leaves 2,3,4,5 as None). Two cursors
# stepping one index at a time cannot do that.
HIGH_PERFECT = 1.0      # NMI of a partition identical to the reference -- exactly 1.0,
                        # not approximately: sklearn's NMI of a labelling against itself
                        # or its complement is bit-exact 1.0.
LOW_FLOOR = 0.10        # NMI at or below which a point counts as "no recovery"


class MiddleOut:
    """Which p index to compute next, and which ends of the grid to extrapolate.

        trav = MiddleOut(len(p_values), early_stop)
        for idx in trav.indices():
            row = compute(p_values[idx])
            trav.observe(idx, row["nmi"])
        # trav.high_fill / trav.low_fill are the indices to extrapolate, and
        # trav.high_trigger / trav.low_trigger the measured points that decided it

    The middle seed counts towards both sides, so a grid whose centre already recovers
    perfectly stops the high side one step later rather than measuring it twice.

    BOTH endpoints are ALWAYS measured, even when a side stops before reaching them.
    ``p_min`` is the bottom of the axis the recovery figure is drawn to show, and one
    extra solve there is the cheapest on the grid (sparse below density 0.1, dense
    above). ``p_max`` is measured for a different reason: at p=1 the sampler is the
    identity, so the arms disagree about what "exact" means -- the Fiedler arms
    short-circuit to 1.0 while the distance arm's ARPACK solve lands ~1e-4 away. Keeping
    the endpoint measured means a ``--no-early-stop`` run and an early-stopped one agree
    EXACTLY, everywhere, instead of agreeing to solver tolerance at one point.

    ``early_stop=0`` disables every stop and hands back a plain ascending grid, so a
    ``--no-early-stop`` run computes and logs its p's in the order it always did.

    ``traversal="ascending"`` is the pre-2026-09 shape: walk 0..n-1 and stop only on the
    high side, and -- the one exception to the endpoint rule above -- fill through
    ``p_max`` rather than anchoring it. It exists because the eta-pool cache holds 503
    sweeps produced that way and they back a published figure; anchoring the endpoint
    would make a fresh run disagree with those entries under their own key. It is a cache key (``sweep_cache._build_config``),
    so an ascending sweep and a middle-out one can never be pooled.
    """

    TRAVERSALS = ("middle_out", "ascending")

    def __init__(self, n_p: int, early_stop: int = 2,
                 high: float = HIGH_PERFECT, low: float = LOW_FLOOR,
                 traversal: str = "middle_out"):
        if n_p < 1:
            raise ValueError(f"n_p must be >= 1, got {n_p}")
        if early_stop < 0:
            raise ValueError(f"early_stop must be >= 0, got {early_stop}")
        if traversal not in self.TRAVERSALS:
            raise ValueError(f"traversal must be one of {self.TRAVERSALS}, "
                             f"got {traversal!r}")
        self.traversal = traversal
        self.n_p = int(n_p)
        self.early_stop = int(early_stop)
        self.high = float(high)
        self.low = float(low)
        # an ascending sweep has no low side and no middle seed: index 0 is its start,
        # and every point it visits belongs to the high run
        self.mid = self.n_p // 2 if traversal == "middle_out" else 0
        self.high_fill: List[int] = []
        self.low_fill: List[int] = []
        self.high_trigger: List[int] = []
        self.low_trigger: List[int] = []
        self._seen: Dict[int, float] = {}
        self._hi_run: List[int] = []     # newest run of high-side indices at/above `high`
        self._lo_run: List[int] = []     # newest run of low-side indices at/below `low`

    def observe(self, idx: int, nmi: float) -> None:
        """Record the NMI measured at ``idx``. A NaN breaks both runs, as it must:
        a p whose every replicate failed says nothing about its neighbours."""
        v = float(nmi)
        self._seen[idx] = v
        # visit order IS p order on each side (ascending high, descending low), so
        # appending keeps each run in the order the guardrail has to read it
        if idx >= self.mid:
            self._hi_run = self._hi_run + [idx] if v >= self.high else []
        if idx <= self.mid:
            self._lo_run = self._lo_run + [idx] if v <= self.low else []

    @property
    def measured(self) -> List[int]:
        return sorted(self._seen)

    def indices(self):
        """Yield the indices to compute, stopping each side on its own rule."""
        if self.early_stop <= 0:
            yield from range(self.n_p)
            return

        if self.traversal == "ascending":
            for idx in range(self.n_p):
                if len(self._hi_run) >= self.early_stop:
                    self.high_trigger = self._hi_run[-self.early_stop:]
                    self.high_fill = list(range(idx, self.n_p))
                    return
                yield idx
            return

        yield self.mid
        hi, lo = self.mid + 1, self.mid - 1
        hi_done = lo_done = False
        while not (hi_done and lo_done):
            if not hi_done:
                if hi >= self.n_p:
                    hi_done = True
                elif len(self._hi_run) >= self.early_stop:
                    hi_done = True
                    self.high_trigger = self._hi_run[-self.early_stop:]
                    self.high_fill = list(range(hi, self.n_p - 1))
                    if self.n_p - 1 not in self._seen:
                        yield self.n_p - 1       # the p_max anchor
                else:
                    yield hi
                    hi += 1
            if not lo_done:
                if lo < 0:
                    lo_done = True
                elif len(self._lo_run) >= self.early_stop:
                    lo_done = True
                    self.low_trigger = self._lo_run[-self.early_stop:]
                    self.low_fill = list(range(1, lo + 1))
                    if 0 not in self._seen:
                        yield 0                  # the p_min anchor
                else:
                    yield lo
                    lo -= 1


def _bipartition_agreement(partition_ref: np.ndarray, partition_avg: np.ndarray) -> float:
    matches_direct = int(np.sum(partition_ref == partition_avg))
    matches_flipped = int(np.sum(partition_ref != partition_avg))
    return 100.0 * max(matches_direct, matches_flipped) / len(partition_ref)


def _kmeans_bipartition(v: np.ndarray, min_split: int = 1) -> np.ndarray:
    """Data-driven 2-way split of a Fiedler vector via k-means (k=2).

    The Ng-Jordan-Weiss recipe: cluster the 1-D Fiedler entries into two
    groups, splitting at the variance-minimizing natural break rather than a
    fixed threshold (τ=0). Label assignment is arbitrary — downstream metrics
    (ARI / NMI / flip-max agreement) are label-invariant. ``min_split`` is a
    guard: if k-means peels off a cluster smaller than it, we log and keep the
    split (rare with the normalized Laplacian, which is the hypothesis).
    """
    from sklearn.cluster import KMeans

    vv = np.asarray(v, dtype=float).reshape(-1, 1)
    labels = KMeans(n_clusters=2, n_init=10, random_state=0).fit_predict(vv)
    partition = labels == 1
    smaller = min(int(partition.sum()), int((~partition).sum()))
    if smaller < min_split:
        log_warning('p_sweep_inner',
                    f"k-means split smaller side {smaller} < min_split {min_split}")
    return partition


def align_fiedler_by_dot_product(
    fiedler_vector: np.ndarray, reference_vector: np.ndarray
) -> np.ndarray:
    """Normalize and sign-align a Fiedler vector against a reference."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        try:
            v_normalized = _normalize_vector(fiedler_vector)
            u_normalized = _normalize_vector(reference_vector)
        except ValueError as e:
            log_warning('align', f"Normalization failed: {e}")
            return fiedler_vector
        dot_product = np.dot(v_normalized, u_normalized)
        if not np.isfinite(dot_product):
            log_warning('align', f"Non-finite dot product: {dot_product}")
            return fiedler_vector
        return -v_normalized if dot_product < 0 else v_normalized


def bootstrap_p_sweep_simple(
    M: np.ndarray,
    fiedler_ref: np.ndarray,
    p_values: List[float],
    bootstrap_reps: int = 50,
    seed: int = 0,
    num_gaps: int = 10,
    min_split: int = 1,
    partition_ref: Optional[np.ndarray] = None,
    early_stop: int = 0,
    low_floor: float = LOW_FLOOR,
    traversal: str = "middle_out",
    partition_method: str = "sigma2",
    laplacian: str = "unnormalized",
    extra_metrics: bool = False,
    progress_cb: Optional[Callable[[int, float, float], None]] = None,
    on_p_result: Optional[Callable[[int, float, Dict[str, Any]], None]] = None,
) -> Dict[str, Any]:
    """Bootstrap p-sweep with precomputed (M, fiedler_ref). Uniform sampling.

    Parameters
    ----------
    M : (n, n) similarity matrix.
    fiedler_ref : (n,) reference Fiedler vector of L(M).
    p_values : sampling rates to sweep.
    bootstrap_reps : bootstrap iterations per p (skipped when p ≈ 1.0).
    seed : base seed; rep i uses seed+i.
    num_gaps, min_split : passed to partition_taxa (sigma2 mode only).
    partition_ref : optional precomputed reference partition. If None,
        derived from fiedler_ref according to ``partition_method``.
    early_stop : if > 0, sweep the grid middle-out (see :class:`MiddleOut`) and stop
        each side after this many consecutive settled points -- NMI == 1.0 on the high
        side, NMI <= ``low_floor`` on the low side. The skipped p's are extrapolated
        under the fill contract below; ``p_values[0]`` is measured either way. 0 sweeps
        the whole grid ascending, as it always did.
    low_floor : the low side's "no recovery" level. Only read when ``early_stop`` > 0.
    traversal : ``"middle_out"`` (default) or ``"ascending"`` -- see :class:`MiddleOut`.
        ``"ascending"`` reproduces the pre-2026-09 shape the eta-pool cache holds.
    partition_method : how to derive partitions from Fiedler vectors.
        - ``"sigma2"`` (default): use ``partition_taxa`` (gap search by σ₂,
          honours ``num_gaps``/``min_split``). Can return degenerate splits
          when an outlier taxon dominates.
        - ``"sign"``: split by sign of the Fiedler vector (``v > 0`` vs
          ``v <= 0``). No σ₂ search; matches the eta-binning convention.
        - ``"kmeans"``: data-driven τ via k-means(k=2) on the Fiedler
          entries (Ng-Jordan-Weiss); splits at the variance-minimizing
          natural break, ``min_split`` guards tiny clusters.
    extra_metrics : also record the linear-algebra diagnostics the original
        pipeline-A runs kept -- ``sigma2_avg_M``/``_S``, mean/median/std of
        ``||S - M||_2`` and of the numerical ranks of ``S`` and ``L(S)``, plus the
        constant ranks of ``M`` and ``L(M)`` (see ``EXTRA_METRIC_KEYS``). Off by
        default, and when off this function's arithmetic is unchanged -- it backs a
        published figure. Costs one extra n x n accumulator and O(n^2) work per
        replicate.
    progress_cb : optional ``f(p_index, p, nmi)`` called after each p value, so a caller can
        drive a progress bar over a sweep that otherwise runs silently for minutes.
    on_p_result : optional ``f(p_index, p, row)`` called with that p's metrics the moment
        it finishes, so a caller can log them live. Separate from ``progress_cb`` because
        the bar only needs to tick, and because every existing caller passes a
        two-argument callback. Without it a tree that takes half an hour reports nothing
        between "loaded" and "done" -- the pre-refactor pipeline logged a line per p.
    laplacian : which Laplacian the per-bootstrap Fiedler is taken from.
        - ``"unnormalized"`` (default): ``L = Deg(S) - S``.
        - ``"normalized"``: ``L_sym = I - D^{-1/2} S D^{-1/2}``.
        The reference Fiedler is supplied by the caller (``fiedler_ref``),
        so pass a reference computed from the matching Laplacian.

    Returns
    -------
    dict with keys: ``p_values``, ``partition_agreement_M``,
    ``partition_ari_M``, ``partition_nmi_M``, ``sign_agreement``,
    ``dot_product``, ``partition_split_ref``, ``reference_partition_quality``,
    and ``measured`` (1.0 where the p was computed, 0.0 where it was extrapolated).
    ``partition_ari_M`` is the Adjusted Rand Index and ``partition_nmi_M``
    is the Normalized Mutual Information between the reference bipartition
    and the bootstrap-averaged bipartition (both label-permutation invariant).

    Fill contract
    -------------
    An extrapolated p carries ONLY what its side's rule supports, and NaN everywhere
    else. Filling a quantity the trigger does not imply is how the mechanisms this one
    replaces went wrong: a perfectly recovered partition says nothing about ``dot`` or
    about ``||S - M||_2``, both of which are continuous in the sub-sample.

    =====================  =========================  ==========================
    quantity               high tail (NMI == 1)       low tail (NMI <= floor)
    =====================  =========================  ==========================
    nmi, ari               1.0                        mean of the triggering p's
    partition_agreement_M  100.0                      mean of the triggering p's
    the partition itself   ``partition_ref``          None
    dot, sign_agreement    NaN                        NaN
    sigma2_avg_S, ranks    NaN                        NaN
      of S and L(S)
    sigma2_avg_M, ranks    their real value -- these are constants of M, not of p
      of M and L(M)
    =====================  =========================  ==========================

    ``on_p_result`` fires in VISIT order, which is middle-out, not ascending. Every
    caller logs ``p`` rather than the index, so the log stays readable; a caller that
    assumes ascending order must pass ``early_stop=0``.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    from spectraltree.spectral_tree_reconstruction import partition_taxa

    if partition_method not in ("sigma2", "sign", "kmeans"):
        raise ValueError(
            f"partition_method must be 'sigma2', 'sign' or 'kmeans', got {partition_method!r}")
    if laplacian not in ("unnormalized", "normalized"):
        raise ValueError(
            f"laplacian must be 'unnormalized' or 'normalized', got {laplacian!r}")

    if partition_ref is None:
        if partition_method == "sigma2":
            partition_ref, ref_quality, ref_split = compute_reference_partition_and_quality(
                fiedler_ref, M, num_gaps=num_gaps, min_split=min_split
            )
        else:
            partition_ref = (_kmeans_bipartition(fiedler_ref, min_split)
                             if partition_method == "kmeans" else fiedler_ref > 0)
            ref_quality = float('nan')
            n_true = int(np.sum(partition_ref))
            n_false = len(partition_ref) - n_true
            ref_split = (min(n_true, n_false), max(n_true, n_false))
    else:
        ref_quality = float('nan')
        n_true = int(np.sum(partition_ref))
        n_false = len(partition_ref) - n_true
        ref_split = (min(n_true, n_false), max(n_true, n_false))

    partition_ref_int = partition_ref.astype(int)

    # Every output is indexed by p, never appended to: the grid is visited middle-out,
    # so the slot a result belongs in is not the number of results so far.
    n_p = len(p_values)
    extra: Dict[str, List[float]] = ({k: [float('nan')] * n_p for k in EXTRA_METRIC_KEYS}
                                     if extra_metrics else {})
    rank_M = rank_L_M = float('nan')
    sigma2_ref_M = float('nan')
    if extra_metrics:
        L_M_full = (compute_normalized_laplacian(M) if laplacian == "normalized"
                    else compute_laplacian(M))
        rank_M = numerical_rank(M)
        rank_L_M = numerical_rank(L_M_full)
        sigma2_ref_M = sigma2_of_partition(M, partition_ref)
        del L_M_full                       # n x n; nothing below needs it

    def _set_extra(idx, sigma2_M, sigma2_S, op_errs, rank_S, rank_LS):
        """One row of diagnostics for one p. Lists may be empty (every rep failed)."""
        if not extra_metrics:
            return
        extra["sigma2_avg_M"][idx] = float(sigma2_M)
        extra["sigma2_avg_S"][idx] = float(sigma2_S)
        for name, vals in (("operator_norm_error", op_errs),
                           ("empirical_rank_S", rank_S),
                           ("empirical_rank_L_S", rank_LS)):
            arr = np.asarray(vals, dtype=float)
            for stat, fn in (("mean", np.mean), ("median", np.median), ("std", np.std)):
                extra[f"{stat}_{name}"][idx] = (float(fn(arr)) if arr.size
                                                else float('nan'))
        extra["empirical_rank_M"][idx] = rank_M
        extra["empirical_rank_L_M"][idx] = rank_L_M

    def _emit(idx: int, p: float) -> None:
        """Hand the caller the row just written for this p."""
        if on_p_result is None:
            return
        row = {"nmi": partition_nmi_M[idx], "ari": partition_ari_M[idx],
               "agreement": partition_agreement_M[idx],
               "sign_agreement": sign_agreement[idx], "dot": dot_product[idx],
               "partition": partitions[idx], "measured": measured[idx]}
        for k, vals in extra.items():
            row[k] = vals[idx]
        on_p_result(idx, float(p), row)

    partition_agreement_M: List[float] = [float('nan')] * n_p
    partition_ari_M: List[float] = [float('nan')] * n_p
    partition_nmi_M: List[float] = [float('nan')] * n_p
    # the bootstrap-averaged partition per p, so a caller can score it against another
    # reference (e.g. the true tree's split) without paying for the sweep twice
    partitions: List[Optional[np.ndarray]] = [None] * n_p
    sign_agreement: List[float] = [float('nan')] * n_p
    dot_product: List[float] = [float('nan')] * n_p
    measured: List[float] = [0.0] * n_p

    trav = MiddleOut(n_p, early_stop=early_stop, low=low_floor, traversal=traversal)
    for idx in trav.indices():
        p = float(p_values[idx])
        measured[idx] = 1.0

        if p >= 0.9999:
            # S is M exactly, so this is not an extrapolation: every quantity below is
            # the reference's own, including the ones a fill would have to leave NaN.
            partition_agreement_M[idx] = 100.0
            partition_ari_M[idx] = 1.0
            partition_nmi_M[idx] = 1.0
            sign_agreement[idx] = 100.0
            dot_product[idx] = 1.0
            partitions[idx] = partition_ref.copy()
            _set_extra(idx, sigma2_ref_M, sigma2_ref_M, [0.0], [rank_M], [rank_L_M])
            trav.observe(idx, 1.0)
            if progress_cb is not None:
                progress_cb(idx, p, 1.0)
            _emit(idx, p)
            continue

        aligned: List[np.ndarray] = []
        op_errs: List[float] = []
        rank_S_vals: List[float] = []
        rank_LS_vals: List[float] = []
        S_sum = None
        for i in range(bootstrap_reps):
            S = _subsample_matrix_entries(M, p, seed=seed + i, builder=None)
            try:
                L_S = (compute_normalized_laplacian(S) if laplacian == "normalized"
                       else compute_laplacian(S))
                f_est = compute_fiedler_from_laplacian(L_S)
            except Exception as e:
                log_warning('p_sweep_inner', f"Fiedler failed at p={p:.4g}, rep={i}: {e}")
                continue
            aligned.append(align_fiedler_by_dot_product(f_est, fiedler_ref))
            if extra_metrics:
                op_errs.append(operator_norm_error(S, M))
                rank_S_vals.append(numerical_rank(S))
                rank_LS_vals.append(numerical_rank(L_S))
                # sigma2_avg_S is read off the AVERAGED sub-sample, as in the original
                S_sum = S.astype(float).copy() if S_sum is None else S_sum + S

        if not aligned:
            # every slot already holds NaN; observe() sees the NaN and breaks both runs,
            # so a failed p never counts towards a stop
            _set_extra(idx, float('nan'), float('nan'), [], [], [])
            trav.observe(idx, float('nan'))
            _emit(idx, p)
            continue

        try:
            v_avg = _normalize_vector(np.mean(aligned, axis=0))
        except ValueError:
            v_avg = np.zeros_like(fiedler_ref)

        try:
            if partition_method == "sigma2":
                partition_avg = partition_taxa(v_avg, M, num_gaps, min_split)
            elif partition_method == "kmeans":
                partition_avg = _kmeans_bipartition(v_avg, min_split)
            else:
                partition_avg = v_avg > 0
            agr_M = _bipartition_agreement(partition_ref, partition_avg)
            ari_M = float(adjusted_rand_score(partition_ref_int, partition_avg.astype(int)))
            nmi_M = float(normalized_mutual_info_score(partition_ref_int, partition_avg.astype(int)))
        except Exception as e:
            log_warning('p_sweep_inner', f"partition_agreement failed at p={p:.4g}: {e}")
            agr_M = float('nan')
            ari_M = float('nan')
            nmi_M = float('nan')
            partition_avg = None
        partitions[idx] = (None if partition_avg is None
                           else np.asarray(partition_avg).astype(bool))
        partition_agreement_M[idx] = float(agr_M)
        partition_ari_M[idx] = ari_M
        partition_nmi_M[idx] = nmi_M
        sign_agreement[idx] = float(compute_sign_agreement(fiedler_ref, v_avg))
        try:
            dot_product[idx] = float(compute_fiedler_dot_product(fiedler_ref, v_avg))
        except ValueError as e:
            log_warning('p_sweep_inner', f"dot_product failed at p={p:.4g}: {e}")
        if extra_metrics:
            if partition_avg is None:
                _set_extra(idx, float('nan'), float('nan'), op_errs,
                           rank_S_vals, rank_LS_vals)
            else:
                part = np.asarray(partition_avg).astype(bool)
                s_avg = S_sum / max(len(aligned), 1) if S_sum is not None else None
                _set_extra(idx, sigma2_of_partition(M, part),
                           float('nan') if s_avg is None
                           else sigma2_of_partition(s_avg, part),
                           op_errs, rank_S_vals, rank_LS_vals)
        trav.observe(idx, nmi_M)
        if progress_cb is not None:
            progress_cb(idx, float(p), float(nmi_M))
        _emit(idx, p)

    # --- extrapolate the ends the traversal skipped ---------------------------------
    # Only nmi / ari / agreement / the partition are written. Everything else keeps the
    # NaN it was preallocated with, which is the fill contract: ``dot``, ``sign``, and
    # every S-side diagnostic are continuous in the sub-sample and are NOT implied by a
    # settled partition. ``run_export._stats`` drops NaN and reports ``n``, so the
    # extrapolated rows announce themselves in curves.csv.
    if trav.high_fill:
        # the high rule is NMI == 1.0 exactly, so the partition IS the reference and
        # 1.0 / 1.0 / 100.0 are deductions, not estimates. Guarded rather than assumed:
        # a relaxed threshold would make them estimates, and then nothing is filled.
        exact = all(partition_nmi_M[i] >= HIGH_PERFECT
                    and partition_agreement_M[i] == 100.0 for i in trav.high_trigger)
        for i in trav.high_fill:
            if exact:
                partition_nmi_M[i] = 1.0
                partition_ari_M[i] = 1.0
                partition_agreement_M[i] = 100.0
                partitions[i] = partition_ref.copy()
            if progress_cb is not None:
                progress_cb(i, float(p_values[i]),
                            float(partition_nmi_M[i]))
            _emit(i, float(p_values[i]))
    if trav.low_fill:
        # the low rule says only "still on the floor", so the floor's own measured level
        # is carried down rather than a baseline invented from the metric's definition
        held = {name: float(np.mean([vals[i] for i in trav.low_trigger]))
                for name, vals in (("nmi", partition_nmi_M),
                                   ("ari", partition_ari_M),
                                   ("agreement", partition_agreement_M))}
        for i in trav.low_fill:
            partition_nmi_M[i] = held["nmi"]
            partition_ari_M[i] = held["ari"]
            partition_agreement_M[i] = held["agreement"]
            if progress_cb is not None:
                progress_cb(i, float(p_values[i]),
                            float(partition_nmi_M[i]))
            _emit(i, float(p_values[i]))

    return {
        "measured": measured,
        "p_values": list(p_values),
        "partition_agreement_M": partition_agreement_M,
        "partition_ari_M": partition_ari_M,
        "partition_nmi_M": partition_nmi_M,
        "sign_agreement": sign_agreement,
        "dot_product": dot_product,
        "partition_split_ref": ref_split,
        "reference_partition_quality": ref_quality,
        "partitions": partitions,
        **extra,
    }
