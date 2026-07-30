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
from typing import Dict, List, Optional, Sequence  # noqa: E402

import numpy as np  # noqa: E402

from ..utils.benchmark_cache import (  # noqa: E402
    OPERATORS, BenchmarkConfig, append_screen_row, check_meta, load_sweep,
    read_screen_table, save_sweep, sweep_path, update_meta, write_screen_npz,
)
from ..utils.benchmark_plots import (  # noqa: E402
    METHODS, plot_recovery_curve, plot_scale, summarize,
)
from .benchmark_workers import run_pool, screen_one, sweep_one  # noqa: E402

__all__ = ["BenchmarkConfig", "run_benchmark_pipeline", "screen", "build_cohort",
           "sweep", "aggregate_and_plot"]


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

    if any(r.get("valid_S") is not None for r in rows):
        counts = {op: sum(1 for r in rows if r.get(f"valid_{op}"))
                  for op in OPERATORS}
        both = sum(1 for r in rows if r.get("valid_S") and r.get("valid_D"))
        print(f"  valid: S={counts['S']}  D={counts['D']}  "
              f"L_sym={counts['Lsym']}  |  S∩D={both} / {len(rows)} screened")
    else:
        print("  no newick trees — validity gate skipped for every tree")
    print(f"  table: {table}")
    return rows


def build_cohort(loader, rows: List[dict], tree_ids: Sequence[str],
                 cfg: BenchmarkConfig) -> List[str]:
    """Stage 2 — trees valid under BOTH Fiedler-on-S and Griffing-on-D."""
    graded = [r for r in rows if r.get("valid_S") is not None]
    if not graded:
        print(f"cohort: no validity information — taking the first "
              f"{cfg.n_compare} trees")
        return list(tree_ids[:cfg.n_compare])

    both = [r["tree"] for r in graded if r.get("valid_S") and r.get("valid_D")]
    groups: Dict[Optional[str], List[str]] = {}
    for t in both:
        groups.setdefault(loader.group_of(t), []).append(t)

    if len(groups) > 1 and None not in groups:
        # Balance across generators — plain sorting would put every 'bd_*' before
        # every 'kingman_*' and the cap would drop one generator entirely.
        per = max(1, cfg.n_compare // len(groups))
        cohort = sorted(t for g in groups.values() for t in sorted(g)[:per])
        print(f"cohort: {len(both)} valid in both; balanced across "
              f"{len(groups)} generators -> {len(cohort)}")
    else:
        cohort = sorted(both)[:cfg.n_compare]
        print(f"cohort: {len(both)} valid in both -> {len(cohort)} after the "
              f"n_compare={cfg.n_compare} cap")
    return cohort


def sweep(loader, cohort: Sequence[str], out_dir: Path, cfg: BenchmarkConfig,
          n_workers: int) -> Dict[str, dict]:
    """Stage 3 — bootstrap p-sweep per tree, three operators, cached per tree."""
    (out_dir / "sweeps").mkdir(parents=True, exist_ok=True)
    results: Dict[str, dict] = {}
    todo = []
    for ti, tid in enumerate(cohort):
        cached = load_sweep(sweep_path(out_dir, tid), cfg)
        if cached is not None:
            results[tid] = cached
        else:
            todo.append((tid, ti))  # ti frozen by cohort order -> stable seeds

    print(f"\nsweep: {len(results)} already done, {len(todo)} to go "
          f"({len(cohort)} trees × {len(cfg.p_values)} p × "
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


def aggregate_and_plot(cohort: Sequence[str], results: Dict[str, dict],
                       out_dir: Path, cfg: BenchmarkConfig,
                       n_taxa: Optional[int] = None) -> dict:
    """Stages 4+5 — stack the per-tree curves, save both figures."""
    kept = [t for t in cohort if t in results]
    if not kept:
        raise RuntimeError("no tree completed the sweep — nothing to plot")
    if len(kept) < len(cohort):
        print(f"  aggregating {len(kept)}/{len(cohort)} trees "
              f"({len(cohort) - len(kept)} failed)")

    curves = {key: np.vstack([results[t][key] for t in kept])
              for key, _, _, _ in METHODS}
    rT = np.array([results[t]["rT"] for t in kept], float)

    np.savez(out_dir / "compare_sweep.npz",
             p_values=np.asarray(cfg.p_values, float),
             tree_ids=np.array(kept, dtype=object), rT=rT,
             nmi_G=curves["G"], nmi_L=curves["L"], nmi_Lsym=curves["Lsym"])

    fig1 = plot_recovery_curve(cfg.p_values, curves,
                               out_dir / "recovery_curve.png", n_taxa=n_taxa)
    fig2 = plot_scale(cfg.p_values, curves, rT, out_dir / "scale_plot.png")

    summary = summarize(cfg.p_values, curves, rT)
    summary["tree_ids"] = kept
    summary["figures"] = [str(fig1), str(fig2)]
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def run_benchmark_pipeline(loader, tree_ids: Sequence[str], out_dir: Path,
                           cfg: BenchmarkConfig,
                           num_workers: Optional[int] = None,
                           source_meta: Optional[dict] = None) -> dict:
    """Screen -> cohort -> sweep -> aggregate -> plot, resumable throughout."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_ids = list(tree_ids)
    n_workers = num_workers or (os.cpu_count() or 1)

    check_meta(out_dir, cfg, {"source": source_meta or {},
                              "n_tree_ids": len(tree_ids),
                              "config": asdict(cfg)})
    print(f"output dir : {out_dir}")
    print(f"workers    : {n_workers}   trees: {len(tree_ids)}")

    rows = screen(loader, tree_ids, out_dir, cfg, n_workers)
    cohort = build_cohort(loader, rows, tree_ids, cfg)
    if not cohort:
        raise RuntimeError("cohort is empty — no tree passed both validity gates")
    update_meta(out_dir, cohort=cohort)  # frozen: per-tree seed is 1000 * index

    results = sweep(loader, cohort, out_dir, cfg, n_workers)
    n_taxa = next((r["n_taxa"] for r in rows if r.get("n_taxa")), None)
    return aggregate_and_plot(cohort, results, out_dir, cfg, n_taxa=n_taxa)
