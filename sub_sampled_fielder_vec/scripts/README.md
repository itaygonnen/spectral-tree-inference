# scripts/

**Everything here is something you run.** Library code lives in `src/`.

That was not true before: six `plot_*` modules and `merge_results.py` were imported by
notebooks and by `src/runners/`, which meant `src/` depended on `scripts/` and this README
had to open with a warning. They now live in `src/plots/` and `src/utils/merge_results.py`.

Each script's **invocation contract differs**, and the filename does not tell you which —
so it is recorded here. Every file has a `__main__` guard, so "it has one" proves nothing.

| Script | How to invoke | Writes | Serves |
|---|---|---|---|
| `build_eta_pool.py` | `--n <n> …` | `cache/pool_sample` | Figs 3, 5, 8, 9 |
| `build_eta_pool_parallel.py` | `--n <n> …` (+ `mp.Pool`) | `cache/pool_sample` | as above, faster |
| `build_sweeps.py` | `--ns 3000 6000 --methods kmeans` | `cache/bootstrap_sweep` | Figs 3, 8, 9 |
| `sync_paper_figures.py` | `--check` / `--write-map` | `analysis/PAPER_MAP.md` | figure provenance |
| `collate_open_items.py` | no args (argparse, all optional) | `docs/overleafs/v9/OPEN_ITEMS.md` | the register |
| `run_real_data_benchmark.py` | `--out-dir …` etc. | `results/runs/real_data_benchmark/` | real-data notebook |
| `plot_real_data_benchmark.py` | `--metric …` | figures beside the results | real-data notebook |
| `nj_recompute_normalized_metrics.py` | **`<sweep_dir>`** positionally | rewrites NJ metrics in place | NJ notebooks |
| `plot_bpart_eta_pool.py` | **`<run_dir>`** positionally | PNGs in the run dir | bpart notebook |
| `plot_fiedler_overlay.py` | **`<run_dir>`** positionally | PNGs in the run dir | — (see note) |
| `plot_griffing_overlay.py` | **`<run_dir>`** positionally | `results/notebooks/05_nj_distance/…` | griffing |
| `plot_distance_vs_similarity_grid.py` | **no args** (paths hardcoded) | PNG under `results/runs/balanced_binary/` | exploratory |
| `interactive_run.py` | **interactive menu** (`d` = real FASTA cohorts) | `results/runs/…`, `cache/experiment_data`, real-cohort caches | pipeline A + real data |
| `run_benchmark.py` | **interactive menu** — see [docs/BENCHMARK_GUIDE.md](../docs/BENCHMARK_GUIDE.md) | `results/runs/…` + `results.json`, `screen_table.csv` | operator comparison |
| `run_experiment.py` | **edit `SWEEP_CONFIG` at the top, then run** — no CLI | `results/runs/…` | **Fig 7's prerequisite** |
| `run_nj_sweep.py` | **edit `SWEEP_CONFIG`, then run** — no CLI | `results/runs/…`, `cache/distance_matrix` | NJ notebooks |
| `run_snj_sweep.py` | **edit `SWEEP_CONFIG`, then run** — no CLI | `results/runs/…` | SNJ notebook |
| `run_griffing_sweep.py` | **edit `SWEEP_CONFIG`, then run** — no CLI | `results/runs/…` | griffing |
| `run_real_sweep.py` | `--cohort --stage --workers --p-min --p-points --reps --prefix` (`--list` to see cohorts) | `results/real_data/runs/<ts>-<name>/` + `_cache/` | real-cohort screening + sweep |

The four "edit-the-dict" scripts take **no arguments at all**. An earlier version of this
table wrongly listed three of them as `argparse`, and `nj_recompute_normalized_metrics.py`
too — it uses `sys.argv` directly.

`plot_fiedler_overlay.py` (with `src/utils/fiedler_io.py`) reads output from
`src/runners/fiedler_sweep.py`, whose **driver script was never committed** — the
`results/runs/*-fiedler_sweep_*` dirs exist but nothing in the repo reproduces them. The
runner was kept rather than deleted precisely because it is the only producer of data
those two modules read.

Validate the real-cohort path after touching it — it walks every call site and runs a
two-tree sweep in a scratch results root, in about a minute:

```bash
python -m analysis.utils.real_selftest
```

## `cluster/`
Not a script but a kit: `push.sh` (rsync the repo and a cohort to a Linux
box), `bootstrap.sh` (venv + `requirements-cluster.txt` + `pip install -e . --no-deps`)
and a README with the ssh workflow. `run_real_sweep.py` is the non-interactive twin of
`interactive_run.py`'s real-data branch and is what runs there under `nohup`.

## `validation/`

`fiedler_plateau_validator.py`, `mutation_rate_validation.py` — figure-generating sanity
checks, **not** assertions and not a test suite. Both take no arguments.
`mutation_rate_validation.py` currently fails to import: it wants `utils.utils` and
`utils.experiment_config`, a flat module layout that has never existed in this package
(the same is true at the `pre-cleanup-2026-07-30` tag, so this predates the cleanup).

## `legacy/`

One-shot maintenance and superseded drivers: `migrate_caches.py` (produced the current
cache layout), `cleanup_cache.py`, `run_experiment_distance.py`,
`check_partition_quality.py` (also broken on the same `utils.utils` import),
`run_n512_to_n4096_L10000.py` and `run_n8192_L10000.py` (hardcoded one-size sweeps
superseded by `build_sweeps.py --ns`), and `run_perfox_experiment.py`.
