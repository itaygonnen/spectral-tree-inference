"""Data-source-agnostic core of the operator-comparison benchmark.

Screen -> cohort -> bootstrap p-sweep -> recovery curve + scale plot, exactly as
the ``distance_vs_similarity`` notebooks do it. The *only* thing that varies
between the real and generated twins is the loader callable, so both call
:func:`run_benchmark_pipeline` with a different one.

Every intermediate result is cached per tree as it is produced (see
:mod:`src.utils.benchmark_cache`), so killing the run and re-launching it against
the same output directory resumes rather than restarting.
"""
from __future__ import annotations

import os

# Pin BLAS to one thread per process before numpy is imported — N workers each
# running a dense eigh on an n×n matrix would otherwise oversubscribe the box.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
from dataclasses import asdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Callable, Dict, List, Optional, Sequence, Tuple  # noqa: E402

import numpy as np  # noqa: E402

from ..utils.benchmark_cache import (  # noqa: E402
    SCREEN_OF_METHOD, BenchmarkConfig, append_screen_row, check_meta,
    load_sweep, read_screen_table, save_sweep, sweep_path, update_meta,
    write_screen_npz,
)
from ..utils.benchmark_plots import (  # noqa: E402
    METHODS, plot_recovery_curve, plot_scale, summarize,
)
from ..utils.benchmark_results_json import write_results_json  # noqa: E402
from .benchmark_workers import run_pool, screen_one, sweep_one  # noqa: E402

__all__ = ["BenchmarkConfig", "run_benchmark_pipeline", "screen", "screen_until",
           "build_cohorts", "sweep", "aggregate_and_plot", "valid_by_cell"]


def _error_sink(errors: List[int], res: dict, verb: str) -> bool:
    """Log the first few failures; report whether ``res`` was one."""
    if "error" not in res:
        return False
    errors[0] += 1
    if errors[0] <= 3:
        print(f"  ! {res['tree']} failed to {verb}: {res['error']}", flush=True)
    return True


def screen(loader, tree_ids: Sequence[str], out_dir: Path, cfg: BenchmarkConfig,
           n_workers: int) -> List[dict]:
    """Stage 1 — per-operator validity gate on every tree's full matrix.

    A tree can be valid for one operator and not another, so validity and eta are
    recorded per operator rather than collapsed into a single verdict.
    """
    table = out_dir / "screen_table.csv"
    rows = read_screen_table(table)
    done = {r["tree"] for r in rows}
    todo = [t for t in tree_ids if t not in done]

    print(f"\nscreen: {len(done)} already done, {len(todo)} to go "
          f"({len(tree_ids)} total)")
    if todo:
        errors = [0]

        def _sink(row: dict) -> None:
            if _error_sink(errors, row, "screen"):
                return
            append_screen_row(table, row)
            rows.append(row)

        run_pool(screen_one, todo, n_workers, loader, cfg, _sink, "screened")
        if errors[0]:
            print(f"  ({errors[0]} trees failed to screen and were skipped)")

    order = {t: i for i, t in enumerate(tree_ids)}
    rows.sort(key=lambda r: order.get(r["tree"], len(order)))
    write_screen_npz(out_dir, rows, tree_ids)

    ops = cfg.screen_ops()
    if any(r.get(f"valid_{ops[0]}") is not None for r in rows):
        # Only the selected operators: an unscreened one would read as 0 valid.
        counts = "  ".join(
            f"{m}={sum(1 for r in rows if r.get(f'valid_{SCREEN_OF_METHOD[m]}'))}"
            for m in cfg.methods())
        print(f"  valid: {counts}  / {len(rows)} screened  "
              "(per operator — cohorts are independent)")
    else:
        print("  no newick trees — validity gate skipped for every tree")
    print(f"  table: {table}")
    return rows


def _cell_of(loader, tree_id: str) -> str:
    return loader.group_of(tree_id) or "all"


def valid_by_cell(loader, rows: List[dict],
                  cfg: BenchmarkConfig) -> Dict[str, Dict[str, int]]:
    """Per cell, how many trees came out valid for each selected operator."""
    counts: Dict[str, Dict[str, int]] = {}
    for r in rows:
        cell = counts.setdefault(_cell_of(loader, r["tree"]),
                                 {m: 0 for m in cfg.methods()})
        for m in cfg.methods():
            if r.get(f"valid_{SCREEN_OF_METHOD[m]}"):
                cell[m] += 1
    return counts


def screen_until(loader, tree_ids: Sequence[str], out_dir: Path,
                 cfg: BenchmarkConfig, n_workers: int,
                 target_per_cell: Optional[int] = None,
                 more_ids: Optional[Callable[[Dict[str, int]], List[str]]] = None,
                 max_rounds: int = 12) -> Tuple[List[dict], List[str]]:
    """Screen, and keep drawing replacements for trees that fail the gate.

    An invalid tree is not a result — it is a wasted draw, so when a cell ends up
    with fewer than ``target_per_cell`` valid trees for some operator, ``more_ids``
    is asked for that many more candidates and they are screened too. The loop
    ends when every cell is satisfied, when the supply is exhausted, or after
    ``max_rounds`` (a guard, not an expected exit).

    Returns ``(rows, ids)`` — every row on disk and every id drawn.
    """
    ids = list(tree_ids)
    rows = screen(loader, ids, out_dir, cfg, n_workers)
    if not target_per_cell or more_ids is None:
        return rows, ids

    seen = set(ids)
    for _ in range(max_rounds):
        counts = valid_by_cell(loader, rows, cfg)
        short: Dict[str, int] = {}
        for cell in {_cell_of(loader, t) for t in ids}:
            have = counts.get(cell)
            worst = min(have.values()) if have else 0
            if worst < target_per_cell:
                short[cell] = target_per_cell - worst
        if not short:
            break

        extra = [t for t in (more_ids(short) or []) if t not in seen]
        if not extra:
            print(f"  no replacements left for {len(short)} short cell(s): "
                  + ", ".join(f"{c} needs {n}" for c, n in sorted(short.items())))
            break
        print(f"\nreplacing invalid draws: {sum(short.values())} needed across "
              f"{len(short)} cell(s) -> screening {len(extra)} more")
        seen.update(extra)
        ids += extra
        rows = screen(loader, ids, out_dir, cfg, n_workers)
    return rows, ids


def _balance(loader, trees: Sequence[str], cap: int) -> List[str]:
    """Cap the list without letting one cell crowd out the others."""
    groups: Dict[Optional[str], List[str]] = {}
    for t in trees:
        groups.setdefault(loader.group_of(t), []).append(t)
    if len(groups) > 1 and None not in groups:
        # Plain sorting would put every 'bd_*' before every 'kingman_*' and the
        # cap would drop one cell entirely.
        per = max(1, cap // len(groups))
        return sorted(t for g in groups.values() for t in sorted(g)[:per])
    return sorted(trees)[:cap]


def build_cohorts(loader, rows: List[dict], tree_ids: Sequence[str],
                  cfg: BenchmarkConfig) -> Dict[str, List[str]]:
    """Stage 2 — one cohort per operator, each the trees *that operator* validates.

    Requiring a single tree to be valid under every operator throws away most of
    the draws (and at high eta, nearly all of them). Since each curve is an
    independent claim about its own operator, each gets its own cohort; the
    figures label the per-operator tree count so the difference stays visible.
    """
    ops = cfg.screen_ops()
    graded = [r for r in rows if any(r.get(f"valid_{op}") is not None
                                    for op in ops)]
    if not graded:
        print(f"cohort: no validity information — taking the first "
              f"{cfg.n_compare} trees for every operator")
        return {m: list(tree_ids[:cfg.n_compare]) for m in cfg.methods()}

    cohorts: Dict[str, List[str]] = {}
    for m in cfg.methods():
        ok = [r["tree"] for r in rows if r.get(f"valid_{SCREEN_OF_METHOD[m]}")]
        cohorts[m] = _balance(loader, ok, cfg.n_compare)
        print(f"cohort[{m}]: {len(ok)} valid -> {len(cohorts[m])} "
              f"after the n_compare={cfg.n_compare} cap")
    return cohorts


def sweep(loader, cohorts: Dict[str, Sequence[str]], out_dir: Path,
          cfg: BenchmarkConfig, n_workers: int) -> Dict[str, dict]:
    """Stage 3 — bootstrap p-sweep per tree, cached per tree.

    Each tree is swept only for the operators whose cohort it is in: a tree that
    only passed the L gate does not pay for a Griffing sweep whose curve nothing
    would plot.
    """
    (out_dir / "sweeps").mkdir(parents=True, exist_ok=True)
    union = sorted({t for ts in cohorts.values() for t in ts})
    members = {m: set(ts) for m, ts in cohorts.items()}
    ops_of = {t: [m for m in cfg.methods() if t in members.get(m, ())]
              for t in union}

    results: Dict[str, dict] = {}
    todo = []
    for ti, tid in enumerate(union):
        cached = load_sweep(sweep_path(out_dir, tid), cfg, keys=ops_of[tid])
        if cached is not None:
            results[tid] = cached
        else:
            # ti frozen by union order -> stable seeds across resumes
            todo.append((tid, ti, ops_of[tid]))

    print(f"\nsweep: {len(results)} already done, {len(todo)} to go "
          f"({len(union)} trees × {len(cfg.p_values)} p × "
          f"{cfg.bootstrap_reps} reps)")
    if not todo:
        return results

    errors = [0]

    def _sink(res: dict) -> None:
        if _error_sink(errors, res, "sweep"):
            return
        save_sweep(out_dir, res, cfg)
        results[res["tree"]] = res

    run_pool(sweep_one, todo, n_workers, loader, cfg, _sink, "swept")
    if errors[0]:
        print(f"  ({errors[0]} trees failed to sweep and were dropped)")
    return results


def aggregate_and_plot(cohorts: Dict[str, Sequence[str]],
                       results: Dict[str, dict], out_dir: Path,
                       cfg: BenchmarkConfig,
                       n_taxa: Optional[int] = None,
                       rows: Optional[List[dict]] = None,
                       loader=None, seq_len: Optional[int] = None) -> dict:
    """Stages 4+5 — stack each operator's own curves, save the figures and tables."""
    curves: Dict[str, np.ndarray] = {}
    rT: Dict[str, np.ndarray] = {}
    kept_by: Dict[str, List[str]] = {}
    for m in cfg.methods():
        kept = [t for t in cohorts.get(m, ()) if m in results.get(t, {})]
        if not kept:
            print(f"  ! {m}: no tree completed the sweep — dropped from the figures")
            continue
        missing = len(cohorts.get(m, ())) - len(kept)
        if missing:
            print(f"  aggregating {len(kept)}/{len(kept) + missing} trees for {m}")
        curves[m] = np.vstack([results[t][m] for t in kept])
        rT[m] = np.array([results[t]["rT"] for t in kept], float)
        kept_by[m] = kept
    if not curves:
        raise RuntimeError("no tree completed the sweep — nothing to plot")

    np.savez(out_dir / "compare_sweep.npz",
             p_values=np.asarray(cfg.p_values, float),
             **{f"nmi_{m}": c for m, c in curves.items()},
             **{f"rT_{m}": r for m, r in rT.items()},
             **{f"trees_{m}": np.array(kept_by[m], dtype=object)
                for m in curves})

    fig1 = plot_recovery_curve(cfg.p_values, curves,
                               out_dir / "recovery_curve.png", n_taxa=n_taxa)
    fig2 = plot_scale(cfg.p_values, curves, rT, out_dir / "scale_plot.png")

    # The per-p metric table: every metric at every p, per size and operator.
    tables: List[Path] = []
    if rows is not None and loader is not None:
        cell_of = {r["tree"]: loader.group_of(r["tree"]) for r in rows}
        n_taxa_of = {r["tree"]: r.get("n_taxa") or n_taxa for r in rows}
        tables = write_results_json(out_dir, cohorts, results, cfg, cell_of,
                                    n_taxa_of, seq_len=seq_len)

    summary = summarize(cfg.p_values, curves, rT)
    summary["tree_ids"] = kept_by
    summary["figures"] = [str(fig1), str(fig2)]
    summary["tables"] = [str(t) for t in tables]
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def run_benchmark_pipeline(loader, tree_ids: Sequence[str], out_dir: Path,
                           cfg: BenchmarkConfig,
                           num_workers: Optional[int] = None,
                           source_meta: Optional[dict] = None,
                           target_per_cell: Optional[int] = None,
                           more_ids: Optional[Callable[[Dict[str, int]],
                                                       List[str]]] = None) -> dict:
    """Screen -> cohorts -> sweep -> aggregate -> plot, resumable throughout.

    ``target_per_cell`` + ``more_ids`` turn the screen into a draw-until-satisfied
    loop: trees that fail an operator's validity gate are replaced instead of
    simply reducing the sample.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_ids = list(tree_ids)
    n_workers = num_workers or (os.cpu_count() or 1)

    check_meta(out_dir, cfg, {"source": source_meta or {},
                              "n_tree_ids": len(tree_ids),
                              "config": asdict(cfg)})
    print(f"output dir : {out_dir}")
    print(f"workers    : {n_workers}   trees: {len(tree_ids)}")

    rows, tree_ids = screen_until(loader, tree_ids, out_dir, cfg, n_workers,
                                  target_per_cell, more_ids)
    cohorts = build_cohorts(loader, rows, tree_ids, cfg)
    if not any(cohorts.values()):
        raise RuntimeError("every cohort is empty — no tree passed any validity "
                           "gate")
    # frozen: per-tree seed is 1000 * index in the sorted union
    update_meta(out_dir, cohorts=cohorts)

    results = sweep(loader, cohorts, out_dir, cfg, n_workers)
    n_taxa = next((r["n_taxa"] for r in rows if r.get("n_taxa")), None)
    seq_len = (source_meta or {}).get("seq_len")
    return aggregate_and_plot(cohorts, results, out_dir, cfg, n_taxa=n_taxa,
                              rows=rows, loader=loader, seq_len=seq_len)
