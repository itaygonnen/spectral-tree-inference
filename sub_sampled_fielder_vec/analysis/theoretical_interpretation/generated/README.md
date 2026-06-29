# generated

**Data source:** **dendropy** simulates a tree (Kingman coalescent or birth-death) and evolves
sequences along it (JC69); `S`/`D` are built from those sequences. This is most of the
experiments. The `eta_pool_sweeps/` notebooks read a pre-built pool of cached Kingman samples
binned by imbalance η.

| Topic | Notebook | What it produces |
|---|---|---|
| `eta_pool_sweeps/` | `eta_pool_sweep` | 3 operators (sign / σ₂-gap / k-means L_sym) on the η-pool. **Fig 1** p\* vs n; **Fig 2** recovery grid + p\* summary. `METRIC` and `AGG` toggles in §1 replace the old per_n / per_n_lsym / per_n_median notebooks |
| `eta_pool_sweeps/` | `fiedler_tree_partition_by_eta` | Fiedler vector & tree partition by η; coherence vs imbalance; validates the two η lower bounds |
| `eta_pool_sweeps/` | `kingman_threshold_vs_theory` | Empirical p\* on JC69 + Kingman vs non-balanced CBM theory |
| `sampling_methods/` | `distance_vs_similarity` | Sub-sampling `S` directly vs sub-sampling `D` then `S = exp(−α·D̂)` |
| `tree_reconstruction/` | `snj_subsampling` | SNJ ‖R−R̂‖₂, σ₂ separation, RF distance vs p |
| `nj_distance/` | `nj_subsampling_balanced` | NJ ‖D−D̂‖₂, Q-criterion, RF distance vs p on balanced binary |
| `nj_distance/` | `nj_subsampling_nonbalanced` | NJ on unbalanced birth-death trees |
| `distance_vs_similarity/` | `simulation_distance_vs_similarity` | 600-tree (Kingman + birth-death) recovery: Fiedler-on-S vs Griffing-on-D vs L_sym; screens, η histograms, recovery curves, p\* vs η / diameter. **Generated twin of the `real/` notebook.** Caches: `gen_*.npz` (stale scratch from prior runs is parked in `distance_vs_similarity/legacy/`) |
