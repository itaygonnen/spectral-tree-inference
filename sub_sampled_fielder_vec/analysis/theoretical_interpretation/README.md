# theoretical_interpretation

Paper-figure notebooks for the Fiedler-sub-sampling thesis (V.04). Each notebook is
self-contained: it builds (or loads from disk cache) a similarity / distance
matrix, runs a bootstrap p-sweep, and emits one or two figures plus optional
CSVs to `sub_sampled_fielder_vec/results/notebooks/<subdir>/<notebook>/`.

Re-running any cell is cheap — heavy computation is cached under
`sub_sampled_fielder_vec/cache/{full_matrix,sweep_trial,pool_sample,experiment_data,bootstrap_sweep}/`
through the unified `src.cache_io` API.

## Notebooks

| Subdir | Notebook | What it verifies |
|---|---|---|
| `01_cbm_theory/` | `balanced_binary_threshold` | p\* = C·log(n)/n on symmetric binary tree (α=0.9) |
| `01_cbm_theory/` | `decay_cbm_features` | HBM (decay-CBM) identities: μ(U), λ₂, Lemma 0.4 slack across η,m |
| `01_cbm_theory/` | `nonbalanced_flat_cbm` | p\* ∝ η(1+η)³ log(m) / [m(ρ − η·S_out)²] |
| `02_real_data_sweeps/` | `kingman_threshold_vs_theory` | Empirical p\* on JC69 + Kingman vs non-balanced CBM theory |
| `02_real_data_sweeps/` | `eta_pool_sweep_per_n` | Bootstrap p-sweeps on η-pooled Kingman samples, η ∈ {1,5,10,15} |
| `03_sampling_methods/` | `distance_vs_similarity` | Subsampling S directly vs subsampling D then S = exp(−α·D̂) |
| `03_sampling_methods/` | `nnm_vs_ipw` | Soft-Impute (NNM) vs IPW: p\* scaling and runtime |
| `04_tree_reconstruction/` | `snj_subsampling` | SNJ ‖R−R̂‖₂, σ₂ separation, RF distance vs p |
| `04_tree_reconstruction/` | `stdr_partition_recovery` | STDR recursive partition: first-layer Fiedler, bipartition Jaccard, ARI |
| `05_nj_distance/` | `nj_subsampling_balanced` | NJ ‖D−D̂‖₂, Q-criterion, RF distance vs p on balanced binary |
| `05_nj_distance/` | `nj_subsampling_nonbalanced` | NJ on unbalanced birth-death trees |
| `05_nj_distance/` | `griffing_distance_partition` | Griffing's $J D J$ leading-eigvec partition under sub-sampling (kingman: 50× lower p\* at large n) |
| `05_nj_distance/` | `fiedler_first_layer_recovery` | First-layer Fiedler sign-agreement on JC-similarity under uniform IPW; production sweep via `src/runners/fiedler_sweep.py` |

## Shared utilities

`utils/` is a sibling package imported by every notebook (and by
`scripts/build_eta_pool.py` / `scripts/build_eta_pool_parallel.py`). Do **not** move it.

- `balanced_binary`, `block_model` — closed-form S matrices and population Fiedler vectors
- `linalg_features` — μ(U), λ₂, spectral gap, Lemma 0.4 row
- `tree_features` — imbalance η, structural margin ρ, n_min, HBM helpers
- `spectral` — Fiedler of S, uniform mask, IPW, NNM (Soft-Impute)
- `recovery` — sign-agreement, p\* threshold detector
- `sweep` — bootstrap p-sweep loop with per-trial disk caching
- `cache` — two-level disk cache (full S + sub-sample) with `.complete` sentinels
- `plotting` — C-constant fit and two-panel recovery+scaling figure

## Outputs

Per-notebook CSV/PNG outputs land in `results/notebooks/<subdir>/<notebook>/`, e.g.
`results/notebooks/04_tree_reconstruction/stdr_partition_recovery/agg.csv`.
Each output dir also stores a `config.json` sidecar — re-running the notebook
with the same config loads the cached CSVs instead of recomputing.

The `05_nj_distance/` notebooks (and any production sweep launched via
`scripts/run_*_sweep.py`) write under `results/runs/<timestamp>-<run_name>/`
or `results/runs/<tree_model>/<sampling_method>/<timestamp>-<run_name>/` for the
`ExperimentRunner` pipeline. All run paths go through
`src.cache_io.run_dir(...)`.
