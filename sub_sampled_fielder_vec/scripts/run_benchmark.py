"""Interactive launcher for the operator-comparison recovery benchmark.

Runs the pipeline behind the ``distance_vs_similarity`` notebooks — validity
screen, both-valid cohort, bootstrap p-sweep, recovery curve, scale plot —
against either downloaded FASTA alignments or simulated trees, without opening a
notebook. Answer a handful of prompts (Enter takes the default everywhere) and it
runs to completion on every core of the machine.

Everything is cached per tree as it is produced, so an interrupted run can be
resumed by pointing it at the same output directory.

    python scripts/run_benchmark.py
"""
from __future__ import annotations

import os

# Must precede numpy, here and in the pipeline module: workers each run a dense
# eigendecomposition, and multi-threaded BLAS on top of that oversubscribes.
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import sys  # noqa: E402
from datetime import datetime  # noqa: E402
from pathlib import Path  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = PROJECT_ROOT.parent
for _p in (str(PROJECT_ROOT), str(REPO_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402

from src.utils.dataset_scan import (  # noqa: E402
    DatasetEntry, entry_from_paths, scan_roots,
)
from src.utils.interactive_ui import (  # noqa: E402
    confirm, get_input, get_menu_choice, print_divider, print_error,
    print_header, print_logo, print_option, print_success, print_warning,
)

DEFAULT_ROOTS = [REPO_ROOT / "data", PROJECT_ROOT / "data"]
# Informational only (printed as a hint). The former top-level synthesized/ dir was
# dissolved into per-section subdirs, so point at the notebook tree itself.
SYNTHESIZED_DIR = PROJECT_ROOT / "analysis"

SOURCES = [
    "real         — downloaded FASTA (+ optional Newick) trees",
    "generated    — simulated trees (choose model, sizes, eta bins)",
    "synthesized  — CBM-theory notebooks (not part of this tool)",
]


# ---------------------------------------------------------------------------
# Step 2a — real datasets
# ---------------------------------------------------------------------------

def choose_real_dataset() -> DatasetEntry | None:
    """Scan the known roots, print what was found, let the user pick or type a path."""
    print_header("Step 2 — pick a dataset")
    roots = list(DEFAULT_ROOTS)
    extra = get_input("Extra data root to scan (Enter to skip)", default="")
    if extra:
        roots.append(Path(extra).expanduser())

    print(f"\nscanning: {', '.join(str(r) for r in roots)}")
    entries = scan_roots(roots)

    if entries:
        print()
        for i, e in enumerate(entries, 1):
            print_option(str(i), e.describe(), highlight=e.runnable)
            print(f"      fasta : {e.fasta_dir or '—'}")
            print(f"      newick: {e.newick_dir or '— (validity screen skipped)'}")
    else:
        print_warning("no fasta/ directories found under those roots")

    print()
    choice = get_input("Choice — a number, or a path to a fasta directory",
                       default="1" if entries else "")
    if not choice:
        print_error("nothing selected")
        return None

    if choice.isdigit() and 1 <= int(choice) <= len(entries):
        entry = entries[int(choice) - 1]
    else:
        newick = get_input("Newick directory (Enter for none — skips the screen)",
                           default="")
        try:
            entry = entry_from_paths(Path(choice).expanduser(),
                                     Path(newick).expanduser() if newick else None)
        except NotADirectoryError as exc:
            print_error(str(exc))
            return None

    if not entry.runnable:
        print_error(f"'{entry.label}' has no usable trees ({entry.status})")
        return None
    if not entry.newick_stems:
        print_warning("no newick trees — the validity screen will be skipped and "
                      "the cohort will just be the first n_compare trees")
    return entry


# ---------------------------------------------------------------------------
# Step 2b — generated datasets
# ---------------------------------------------------------------------------

def prompt_generated():
    """Ask the simulation knobs; return (loader, tree_ids, source_meta, plan).

    ``plan`` is None unless the run draws from the eta pool, in which case the
    ids are provisional: the pool is topped up after the confirm, and the real id
    list is whatever ended up on disk.
    """
    from src.utils.generated_prompts import prompt_generated_plan
    from src.runners.benchmark_loaders import GeneratedLoader

    plan = prompt_generated_plan()
    loader = GeneratedLoader(
        seq_len=plan.seq_len,
        n_taxa=plan.n_values[0] if len(plan.n_values) == 1 else None,
        mutation_rate=plan.mutation_rate,
        params=plan.params,
    )
    # second plan is returned unconditionally: the replacement-draw supplier needs
    # it whether or not the eta pool is involved.
    return (loader, plan.tree_ids, plan.meta(),
            (plan if plan.pooled else None), plan)


# ---------------------------------------------------------------------------
# Step 3 — knobs
# ---------------------------------------------------------------------------

def _ask_operators(keys, labels):
    """Which operators to screen, sweep and plot — '1,3' style, all by default."""
    valid = {str(i): k for i, k in enumerate(keys, 1)}
    valid.update({k.lower(): k for k in keys})
    while True:
        raw = get_input("Operators — comma-separated numbers",
                        default=",".join(str(i) for i in range(1, len(keys) + 1)))
        picked = [valid.get(tok.strip().lower())
                  for tok in raw.split(",") if tok.strip()]
        if picked and all(picked):
            chosen = list(dict.fromkeys(picked))
            print(f"  → {', '.join(labels[k] for k in chosen)}")
            return chosen
        print_warning(f"unknown operator in {raw!r}")


def prompt_knobs(n_available: int, *, allow_subset: bool = True,
                 step: str = "3"):
    """Ask the few parameters that matter; every one is Enter-for-default.

    ``allow_subset`` is False for generated runs: there the tree count is derived
    from the model/size/eta grid, and truncating the derived list would drop
    whole cells rather than thin them.
    """
    from src.runners.benchmark_pipeline import BenchmarkConfig
    from src.utils.benchmark_cache import METHOD_KEYS, METHOD_LABELS

    print_header(f"Step {step} — benchmark parameters")

    for i, k in enumerate(METHOD_KEYS, 1):
        print_option(str(i), METHOD_LABELS[k])
    ops = _ask_operators(METHOD_KEYS, METHOD_LABELS)

    n_use = n_available
    if allow_subset:
        n_use = int(get_input(f"How many of the {n_available} trees to use",
                              default=str(n_available)))
    raw_p = get_input("p-values — comma-separated, or Enter for 20 log-spaced "
                      "0.01→1.0", default="")
    p_values = ([float(x) for x in raw_p.split(",") if x.strip()] if raw_p
                else list(np.logspace(-2, 0, 20)))
    reps = int(get_input("Bootstrap replicates per p", default="10"))
    n_compare = int(get_input("Cohort cap (trees carried into the p-sweep)",
                              default="100"))

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    default_out = PROJECT_ROOT / "results" / "runs" / f"{stamp}-benchmark"
    out = get_input("Output directory (an existing one resumes)",
                    default=str(default_out))

    cfg = BenchmarkConfig(p_values=p_values, bootstrap_reps=reps,
                          n_compare=n_compare, operators=ops)
    return cfg, Path(out).expanduser(), n_use


# ---------------------------------------------------------------------------

def show_synthesized() -> None:
    print_header("synthesized — CBM-theory notebooks")
    print(f"  {SYNTHESIZED_DIR}")
    print("\n  These verify the CBM/HBM spectral-gap theory directly "
          "(paper/fig04_hbm_spectral_gap.ipynb,\n  paper/fig02_pstar_synth_cbm.ipynb, …). "
          "They have no bootstrap p-sweep and no recovery\n  curve, so this "
          "benchmark tool does not apply to them — open them as notebooks.")
    print()


def summarize_run(cfg, out_dir: Path, n_trees: int, source_meta: dict,
                  *, provisional: bool = False) -> None:
    print_divider()
    print_header("Ready to run")
    print(f"  → Source     : {source_meta.get('kind')}  "
          f"{source_meta.get('label', '')}")
    print(f"  → Trees      : {n_trees}"
          + ("  (planned — the eta pool is topped up first)" if provisional
             else ""))
    print(f"  → p-values   : {len(cfg.p_values)} points "
          f"[{min(cfg.p_values):.1e}–{max(cfg.p_values):.1e}]")
    print(f"  → Bootstrap  : {cfg.bootstrap_reps} reps   "
          f"cohort cap {cfg.n_compare}")
    print(f"  → Workers    : {os.cpu_count()} (all cores)")
    print(f"  → Output     : {out_dir}")
    print()


def print_results(summary: dict, out_dir: Path) -> None:
    """The numbers, then where everything landed."""
    from src.utils.benchmark_cache import METHOD_LABELS

    print_success("done")
    print(f"\n  {'operator':26s} {'trees':>5s}  {'p*':>6s}  "
          f"{'NMI at smallest p':>17s}")
    for key, m in summary.get("methods", {}).items():
        ps = m.get("p_star_median")
        # METHOD_LABELS, not the figure label — the latter is LaTeX.
        print(f"  {METHOD_LABELS.get(key, key):26s} {m.get('n_trees', '?'):>5}  "
              f"{(f'{ps:.3f}' if ps else '  none'):>6s}  "
              f"{m.get('median_nmi_at_min_p', float('nan')):>17.3f}")
    print(f"\n  r(T) range   : {summary.get('rT_min', 0):.3f} – "
          f"{summary.get('rT_max', 0):.3f}")
    print(f"\n  results in   : {out_dir}")
    for name, what in (("summary.json", "these numbers, machine-readable"),
                       ("recovery_curve.png", "NMI vs p"),
                       ("scale_plot.png", "p* vs r(T)"),
                       ("compare_sweep.npz", "every per-tree curve"),
                       ("screen_table.csv", "per-tree validity + eta"),
                       ("sweeps/", "one npz per tree (the resume cache)")):
        print(f"    {name:20s} {what}")
    print()


def main() -> None:
    print_logo()
    print_header("Step 1 — data source")
    source = get_menu_choice("Which data source?", SOURCES, default_index=0)
    kind = source.split()[0]

    if kind == "synthesized":
        show_synthesized()
        return

    if kind == "real":
        from src.runners.benchmark_loaders import RealLoader
        entry = choose_real_dataset()
        if entry is None:
            return
        loader = RealLoader(fasta_dir=entry.fasta_dir, newick_dir=entry.newick_dir)
        tree_ids = entry.tree_ids
        source_meta = {"kind": "real", "label": entry.label,
                       "fasta_dir": str(entry.fasta_dir),
                       "newick_dir": str(entry.newick_dir or ""),
                       "data_key": {"kind": "real", "label": entry.label}}
        plan = plan_any = None
    else:
        loader, tree_ids, source_meta, plan, plan_any = prompt_generated()

    cfg, out_dir, n_use = prompt_knobs(
        len(tree_ids), allow_subset=(kind == "real"),
        step="3" if kind == "real" else "7")
    tree_ids = tree_ids[:n_use]
    summarize_run(cfg, out_dir, len(tree_ids), source_meta,
                  provisional=plan is not None)
    if not confirm("Confirm and start?"):
        print_warning("cancelled")
        return

    if plan is not None:
        # The expensive part of eta pooling — rejection sampling for the bins the
        # pool is short on — happens here, after the confirm. What comes back is
        # the id list that actually exists on disk.
        from src.utils.generated_prompts import resolve_pool_ids
        tree_ids = resolve_pool_ids(plan)
        if not tree_ids:
            print_error("the eta pool is empty for these parameters and nothing "
                        "could be built — nothing to run")
            return
        source_meta = {**source_meta, "n_trees": len(tree_ids)}
        print_success(f"eta pool ready — {len(tree_ids)} trees")

    target_per_cell = more_ids = None
    if plan_any is not None:
        from src.utils.generated_prompts import make_more_ids
        target_per_cell, more_ids = plan_any.per_cell, make_more_ids(plan_any)

    from src.runners.benchmark_pipeline import run_benchmark_pipeline
    print_divider()
    try:
        summary = run_benchmark_pipeline(loader, tree_ids, out_dir, cfg,
                                         source_meta=source_meta,
                                         target_per_cell=target_per_cell,
                                         more_ids=more_ids)
    except (ValueError, RuntimeError) as exc:
        print_error(str(exc))
        return
    except KeyboardInterrupt:
        print_warning("\ninterrupted — progress is on disk; rerun with the same "
                      f"output directory to resume:\n  {out_dir}")
        return

    print_divider()
    print_results(summary, out_dir)


if __name__ == "__main__":
    main()
