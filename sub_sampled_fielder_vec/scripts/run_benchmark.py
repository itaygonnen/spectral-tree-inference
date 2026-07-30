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
    Colors, confirm, get_input, get_menu_choice, print_divider, print_error,
    print_header, print_logo, print_option, print_success, print_warning,
)

DEFAULT_ROOTS = [REPO_ROOT / "data", PROJECT_ROOT / "data"]
SYNTHESIZED_DIR = PROJECT_ROOT / "analysis" / "theoretical_interpretation" / "synthesized"

SOURCES = [
    "real         — downloaded FASTA (+ optional Newick) trees",
    "generated    — simulated trees (Kingman / birth-death), no download needed",
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
    """Ask the simulation knobs; return (loader, tree_ids, source_meta)."""
    from analysis.theoretical_interpretation.utils.generated_data import (
        GENERATORS, build_ids,
    )
    from src.runners.benchmark_loaders import GeneratedLoader

    print_header("Step 2 — simulation settings")
    n_taxa = int(get_input("Number of taxa per tree", default="1000"))
    seq_len = int(get_input("Sequence length", default="1000"))
    n_per_gen = int(get_input(
        f"Trees per generator ({', '.join(GENERATORS)})", default="300"))

    tree_ids = build_ids(n_per_gen)
    print(f"  → {len(tree_ids)} trees "
          f"({n_per_gen} × {len(GENERATORS)} generators)")
    print(f"  {Colors.CYAN}note{Colors.RESET}: trees are deterministic in their id "
          "(seed = generator offset + index), so reruns rebuild the same data")
    return (GeneratedLoader(n_taxa=n_taxa, seq_len=seq_len), tree_ids,
            {"kind": "generated", "n_taxa": n_taxa, "seq_len": seq_len,
             "n_per_gen": n_per_gen})


# ---------------------------------------------------------------------------
# Step 3 — knobs
# ---------------------------------------------------------------------------

def prompt_knobs(n_available: int):
    """Ask the few parameters that matter; every one is Enter-for-default."""
    from src.runners.benchmark_pipeline import BenchmarkConfig

    print_header("Step 3 — benchmark parameters")

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
                          n_compare=n_compare)
    return cfg, Path(out).expanduser(), n_use


# ---------------------------------------------------------------------------

def show_synthesized() -> None:
    print_header("synthesized — CBM-theory notebooks")
    print(f"  {SYNTHESIZED_DIR}")
    print("\n  These verify the CBM/HBM spectral-gap theory directly "
          "(hbm_spectral_gap_verification.ipynb,\n  nonbalanced_flat_cbm.ipynb, …). "
          "They have no bootstrap p-sweep and no recovery\n  curve, so this "
          "benchmark tool does not apply to them — open them as notebooks.")
    print()


def summarize_run(cfg, out_dir: Path, n_trees: int, source_meta: dict) -> None:
    print_divider()
    print_header("Ready to run")
    print(f"  → Source     : {source_meta.get('kind')}  "
          f"{source_meta.get('label', '')}")
    print(f"  → Trees      : {n_trees}")
    print(f"  → p-values   : {len(cfg.p_values)} points "
          f"[{min(cfg.p_values):.1e}–{max(cfg.p_values):.1e}]")
    print(f"  → Bootstrap  : {cfg.bootstrap_reps} reps   "
          f"cohort cap {cfg.n_compare}")
    print(f"  → Workers    : {os.cpu_count()} (all cores)")
    print(f"  → Output     : {out_dir}")
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
                       "newick_dir": str(entry.newick_dir or "")}
    else:
        loader, tree_ids, source_meta = prompt_generated()

    cfg, out_dir, n_use = prompt_knobs(len(tree_ids))
    tree_ids = tree_ids[:n_use]
    summarize_run(cfg, out_dir, len(tree_ids), source_meta)
    if not confirm("Confirm and start?"):
        print_warning("cancelled")
        return

    from src.runners.benchmark_pipeline import run_benchmark_pipeline
    print_divider()
    try:
        summary = run_benchmark_pipeline(loader, tree_ids, out_dir, cfg,
                                         source_meta=source_meta)
    except (ValueError, RuntimeError) as exc:
        print_error(str(exc))
        return
    except KeyboardInterrupt:
        print_warning("\ninterrupted — progress is on disk; rerun with the same "
                      f"output directory to resume:\n  {out_dir}")
        return

    print_divider()
    print_success(f"done — {summary['n_trees']} trees aggregated")
    for name in ("recovery_curve.png", "scale_plot.png"):
        print(f"  {out_dir / name}")
    print()


if __name__ == "__main__":
    main()
