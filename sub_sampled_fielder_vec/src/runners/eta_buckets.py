"""Which screened trees a sweep actually needs, bucketed by imbalance.

The paper reports at eta in {1, 5, 10, 15}. A sweep that takes every tree clearing the
validity gate pays, in full, for trees that never reach a figure -- the cost is linear in
the selected count (``experiment_run.plan``). This module is the selection that stops
that: put each screened tree in the bucket of the nearest target within ``eta_tol``, take
only as many per bucket as asked for, and drop everything that lands in no bucket.

The dropped trees ARE the saving, so the counts are not a detail to be averaged away.
Measured on the screen caches on disk, bucketing ``eta_S`` at tol=1:

    1000 taxa   32 valid for S  ->  eta01 0, eta05 1, eta10 0, eta15 3,  28 unbucketed
    6000 taxa   69 valid for S  ->  eta01 6, eta05 6, eta10 4, eta15 4,  49 unbucketed

``eta_S`` has median 24.6 at m=1000 and 17.8 at m=6000, i.e. the mass sits ABOVE the
largest target, so a tol=1 window keeps 12-29% of the valid trees. At m=6000 that is the
win (18 swept instead of 69). At m=1000 it is a four-tree run with two empty buckets.
Nothing here decides which of those is acceptable -- it reports occupancy so a caller can
show it and widen ``eta_tol`` before committing.

Membership is :func:`src.utils.eta_pool_cache.closest_target`, unchanged and reused: the
nearest integer target within the tolerance, else no bucket. The windows are therefore
DISJOINT WITH GAPS rather than a partition of the eta axis -- eta=3.4 belongs to nothing
when the targets are 1/5/10/15 -- and that is the intent, not an oversight. The eta pool's
rejection sampler accepts trees by the same rule, which is why a generated (pooled) source
fills its buckets by construction while a real one does not.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from ..utils.eta_pool_cache import bin_name, closest_target

# The paper's buckets and the per-bucket counts the eta-pool figures read with
# (``src.utils.eta_pool_sweep.ETA_TARGETS`` / ``ETA_SAMPLE_CAP``). Restated as parallel
# tuples because a ``RunSpec`` field has to survive ``dataclasses.asdict`` into
# config.json, and a dict keyed by int comes back from JSON keyed by str.
ETA_TARGETS: Tuple[int, ...] = (1, 5, 10, 15)
PER_BUCKET: Tuple[int, ...] = (10, 10, 3, 3)
ETA_TOL: float = 1.0

UNBUCKETED = None       # the dict key for "matched no target"


def window(target: int, eta_tol: float = ETA_TOL) -> Tuple[float, float]:
    """The closed interval a target accepts, clamped at 1.

    Purely for display -- the membership test is ``closest_target``, which also breaks
    ties between overlapping windows. The clamp is not a fudge: ``partition_metrics.eta``
    is larger-clan / smaller-clan and so cannot go below 1, which makes the lower half of
    the eta01 window unreachable at any tolerance. Printing "-1" there invites the reader
    to think trees are being missed below it.
    """
    return (max(1.0, float(target) - float(eta_tol)), float(target) + float(eta_tol))


def label(target: Optional[int]) -> str:
    """``5 -> 'eta05'``, ``None -> 'unbucketed'``. The bucket name is the pool's
    (:func:`eta_pool_cache.bin_name`) so a bucket, a pool directory and a generated tree
    id all spell it the same way and stay greppable together."""
    return "unbucketed" if target is UNBUCKETED else bin_name(target)


def tree_eta(row: dict, cap_ops: Sequence[str]) -> float:
    """The eta this selection weighs for one screen row.

    ``max`` over the operators the GATE names, which is the expression ``select_ids`` and
    ``run_menu._print_status`` already use. Sharing it is what stops a bucket panel
    promising a count the plan then contradicts.
    """
    return max((float(row.get(f"eta_{k}", 0.0)) for k in cap_ops), default=0.0)


def bucket_of(row: dict, cap_ops: Sequence[str],
              eta_targets: Sequence[int] = ETA_TARGETS,
              eta_tol: float = ETA_TOL) -> Optional[int]:
    """The target this tree belongs to, or ``None`` for no bucket."""
    return closest_target(tree_eta(row, cap_ops), [int(t) for t in eta_targets],
                          float(eta_tol))


def bucket_counts(ids: Sequence[str], rows: Dict[str, dict], cap_ops: Sequence[str],
                  eta_targets: Sequence[int] = ETA_TARGETS,
                  eta_tol: float = ETA_TOL) -> Dict[Optional[int], List[str]]:
    """``{target: [tree ids]}`` with an ``UNBUCKETED`` key, every target present.

    Ids are sorted inside each bucket, so which trees a cap keeps does not depend on the
    order the caller enumerated the directory in. Every requested target gets a key even
    when empty -- an empty bucket is a thing to display, not a thing to omit.
    """
    out: Dict[Optional[int], List[str]] = {int(t): [] for t in eta_targets}
    out[UNBUCKETED] = []
    for tree_id in ids:
        row = rows.get(tree_id)
        if row is None:
            continue              # unscreened: not this function's business to guess
        out[bucket_of(row, cap_ops, eta_targets, eta_tol)].append(tree_id)
    for v in out.values():
        v.sort()
    return out


def apply_caps(buckets: Dict[Optional[int], List[str]],
               eta_targets: Sequence[int] = ETA_TARGETS,
               per_bucket: Sequence[int] = PER_BUCKET,
               ) -> Tuple[List[str], Dict[int, Tuple[int, int]]]:
    """``(kept ids, {target: (got, want)} for every short bucket)``.

    The first ``want`` ids of each bucket, which is a plain prefix of a sorted list --
    the same idiom as ``benchmark_pipeline._balance`` and ``resolve_pool_ids``. Nothing
    is ever topped up from outside a window: an eta15 bucket holds trees inside eta15's
    window or it is short, and being short is reported rather than papered over.

    ``UNBUCKETED`` is never kept. That is the whole point -- those are the trees the run
    does not need.
    """
    caps = {int(t): int(c) for t, c in zip(eta_targets, per_bucket)}
    kept: List[str] = []
    shortfall: Dict[int, Tuple[int, int]] = {}
    for target in (int(t) for t in eta_targets):
        want = caps.get(target, 0)
        have = buckets.get(target, [])
        kept.extend(have[:want])
        if len(have) < want:
            shortfall[target] = (len(have), want)
    return sorted(kept), shortfall


def select(ids: Sequence[str], rows: Dict[str, dict], cap_ops: Sequence[str],
           eta_targets: Sequence[int] = ETA_TARGETS,
           per_bucket: Sequence[int] = PER_BUCKET,
           eta_tol: float = ETA_TOL,
           ) -> Tuple[List[str], Dict[Optional[int], List[str]], Dict[int, Tuple[int, int]]]:
    """Bucket then cap. Returns ``(kept, buckets, shortfall)``.

    ``buckets`` comes back alongside the selection so a caller can draw the occupancy it
    selected from without re-deriving it and risking a different answer.
    """
    buckets = bucket_counts(ids, rows, cap_ops, eta_targets, eta_tol)
    kept, shortfall = apply_caps(buckets, eta_targets, per_bucket)
    return kept, buckets, shortfall


def panel_lines(buckets: Dict[Optional[int], List[str]],
                eta_targets: Sequence[int] = ETA_TARGETS,
                per_bucket: Sequence[int] = PER_BUCKET,
                eta_tol: float = ETA_TOL, bar_width: int = 24,
                legend: bool = True) -> List[str]:
    """The occupancy table, as lines ready for ``run_menu._framed``.

    Rendered here rather than in the menu so the non-interactive ``--dry-run`` path shows
    the same panel as the prompt, from the same numbers.

    Columns::

        bucket      the target, named as the eta pool names its bins
        found       screened trees whose eta lands within the tolerance of it
        cap         the most this run will sweep from the bucket
        sweep       what actually runs: min(found, cap)

    The window itself is NOT a column: it is ``target +- eta_tol`` for every row, so a
    per-row copy of it spends twelve characters restating the tolerance once per bucket.
    The rule goes in the header line instead.

    The bar is scaled to the largest row INCLUDING "no bucket", because that row is
    usually the big one and a bar that silently excluded it would hide the saving.

    ``bar_width`` 24 keeps the widest row inside 72 characters, which still fits an
    80-column terminal once ``_framed`` adds its 4 of box. ``legend=False`` drops the
    footer for a redraw that has already shown it once.
    """
    caps = {int(t): int(c) for t, c in zip(eta_targets, per_bucket)}
    rows_out = [(int(t), len(buckets.get(int(t), []))) for t in eta_targets]
    n_unb = len(buckets.get(UNBUCKETED, []))
    peak = max([n for _, n in rows_out] + [n_unb, 1])

    def bar(n: int) -> str:
        filled = int(round(bar_width * n / peak))
        return "█" * filled + "·" * (bar_width - filled)

    # no trailing pad on the header: ``_framed`` sizes its box to the longest line, so a
    # run of spaces there would widen the box past what anything prints
    lines = [f"a tree joins its nearest target if within +-{eta_tol:g}",
             f"{'bucket':<11}{'found':>7}{'cap':>6}{'sweep':>6}"]
    total = 0
    short_any = False
    for target, n in rows_out:
        cap = caps.get(target, 0)
        take = min(n, cap)
        total += take
        # no count in the tag: with the columns named, "found 6 / cap 10" already says
        # it is short, and spelling it out again pushed the row past 80 columns
        mark = ""
        if n < cap:
            short_any = True
            mark = "  short"
        lines.append(f"{label(target):<11}{n:>7}{cap:>6}{take:>6}  {bar(n)}{mark}")
    lines.append(f"{'no bucket':<11}{n_unb:>7}{'-':>6}{0:>6}  {bar(n_unb)}  skipped")
    lines.append(f"{'TOTAL':<11}{'':>7}{sum(caps.values()):>6}{total:>6}"
                 f"  <- will sweep")
    if legend:
        lines += [
            "eta   = larger clan / smaller clan of the tree's reference split",
            "found = screened trees landing within the tolerance of that target",
            "cap   = most to sweep from that bucket;  sweep = min(found, cap)",
        ]
        if short_any:          # only explain the marker when one is on screen
            lines.append("short = fewer found than the cap asks for")
        lines += [
            "no bucket = too far from every target, so skipped. That is the time",
            "            saved; widen the tolerance to catch more of them.",
        ]
    return lines
