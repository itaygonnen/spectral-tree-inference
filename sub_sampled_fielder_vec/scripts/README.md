# scripts/

Entry points and the plotting libraries the notebooks import. Superseded one-shots live
in `legacy/`.

Not everything named `plot_*` is a script: seven of them are **libraries** that
notebooks `import`, so moving or renaming them breaks a notebook. The "Imported by"
column marks those.

## Paper pipeline — the eta-pool route (Figs 3, 5, 7, 8, 9)

| Script | Invocation | Writes | Consumed by |
|---|---|---|---|
| `build_eta_pool.py` | `python scripts/build_eta_pool.py --n <n> …` | `cache/pool_sample` | Figs 3, 5, 8, 9 |
| `build_eta_pool_parallel.py` | same + `mp.Pool` | `cache/pool_sample` | as above, faster |
| `build_sweeps.py` | `python scripts/build_sweeps.py --ns 3000 6000 --methods kmeans` | `cache/bootstrap_sweep` | Figs 3, 8, 9 |
| `run_experiment.py` | `python scripts/run_experiment.py` (edit `SWEEP_CONFIG` at the top; no CLI args) | `results/runs/…`, `cache/experiment_data` | **Fig 7** needs a `kingman_mean`+`uniform` run from this |
| `interactive_run.py` | `python scripts/interactive_run.py` | as above | menu-driven wrapper over the same pipeline |

`build_sweeps.py` defaults to `--methods kmeans`, the paper's operator. `sigma2` is the
expensive one — pass it only as cached-only so a new `n` cannot trigger a fresh run.

## Paper tooling

| Script | Invocation | Purpose |
|---|---|---|
| `sync_paper_figures.py` | `--check` / `--write-map` | Validates figure provenance against `analysis/paper_figures.py`; regenerates `PAPER_MAP.md`. `--check` is read-only and exits non-zero on problems. |
| `collate_open_items.py` | `python scripts/collate_open_items.py` | Builds `docs/overleafs/v9/OPEN_ITEMS.md` from `open-items/*.md`. Edit the fragments, never the collated file. |

## Per-method sweeps (supporting notebooks, no paper figure)

| Script | Invocation | Writes |
|---|---|---|
| `run_nj_sweep.py` | argparse | `results/runs/…`, `cache/distance_matrix` |
| `run_snj_sweep.py` | argparse | `results/runs/…` |
| `run_griffing_sweep.py` | argparse | `results/runs/…` |
| `run_real_data_benchmark.py` | argparse | `results/runs/real_data_benchmark/` |
| `run_benchmark.py` | interactive | `results/runs/…` + `screen_table.csv` |
| `nj_recompute_normalized_metrics.py` | argparse | rewrites NJ metrics in place |

## Plot libraries — imported, not run

| Module | Imported by |
|---|---|
| `plot_nj_p_star_vs_n.py` | `supporting/generated/nj_distance/nj_subsampling_balanced.ipynb` |
| `plot_nj_three_panel.py` | same |
| `plot_snj_p_star_vs_n.py` | `supporting/generated/tree_reconstruction/snj_subsampling.ipynb` |
| `plot_snj_three_panel.py` | same |
| `plot_bpart_eta_grid.py` | `supporting/generated/nj_distance/nj_subsampling_nonbalanced.ipynb`, `src/runners/bpart_synthetic.py` |
| `plot_bpart_overlay.py` | `nj_subsampling_nonbalanced.ipynb` |
| `merge_results.py` | `src/runners/experiment_runner_utils.py` (auto-plot step) |

Standalone plotters, run directly: `plot_bpart_eta_pool.py`,
`plot_distance_vs_similarity_grid.py`, `plot_fiedler_overlay.py`,
`plot_griffing_overlay.py`, `plot_real_data_benchmark.py`.

`plot_fiedler_overlay.py` (with `src/utils/fiedler_io.py`) reads output from
`src/runners/fiedler_sweep.py`. That runner has **no committed driver script** — the
`results/runs/*-fiedler_sweep_*` dirs exist but nothing in the repo reproduces them.
The runner was kept rather than deleted precisely because it is the only producer of
data those two modules read.

## `validation/`

`fiedler_plateau_validator.py`, `mutation_rate_validation.py` — figure-generating
sanity checks, not assertions. Not a test suite.

## `legacy/`

One-shot maintenance and superseded drivers, kept for reference:
`migrate_caches.py` (the cache migration that produced the current layout),
`cleanup_cache.py`, `run_experiment_distance.py`, `check_partition_quality.py`,
`run_n512_to_n4096_L10000.py`, `run_n8192_L10000.py` (hardcoded one-size sweeps
superseded by `build_sweeps.py --ns`), and `run_perfox_experiment.py` (was loose at the
package root).
