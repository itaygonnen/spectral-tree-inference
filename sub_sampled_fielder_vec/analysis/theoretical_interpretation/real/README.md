# real

**Data source:** downloaded **FASTA + Newick** files — the 600-tree, 1000-taxon benchmark.
This is the only place real sequence data is used. FASTA lives under
`data/Datasets/Fasta 1000 taxa/extracted/fasta/`; Newick under
`data/real_datasets/Datasets/1000 taxa/newick/`.

| Topic | Notebook | What it produces |
|---|---|---|
| `distance_vs_similarity/` | `real_data_distance_vs_similarity` | Real 600-tree benchmark: Fiedler-on-S vs Griffing-on-D vs L_sym operator comparison, STDR validity gate, η histograms, recovery curves. The generated analogue is `generated/distance_vs_similarity/simulation_distance_vs_similarity`. Caches: `screen_*_600.npz`, `compare_sweep_*.npz`, `threshold_eta_600.npz`, `rT_cohort.npz` |
