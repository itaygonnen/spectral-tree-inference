# theoretical_interpretation

Paper-figure notebooks for the Fiedler-sub-sampling thesis. **Notebooks are organized by
data source**, because the source is the thing that was confusing when they were named by
topic. Three sources:

| Source | Meaning | Folder |
|---|---|---|
| **synthesized** | the similarity `S` / distance `D` is built **directly** from a closed-form block model (CBM/HBM, balanced binary) — no tree simulation | `synthesized/` |
| **generated** | **dendropy** simulates a tree (Kingman / birth-death) and evolves sequences (JC69); `S`/`D` come from those sequences | `generated/` |
| **real** | downloaded **FASTA + Newick** files (the 600-tree 1000-taxon benchmark) | `real/` |

Within each source, notebooks are grouped by topic. Because some topics (`sampling_methods`,
`tree_reconstruction`, `distance_vs_similarity`) have one notebook per source, those topic
folders appear under more than one source — that is intentional and honest.

```
theoretical_interpretation/
├── utils/        shared helper library (imported by every notebook — do NOT move)
├── figures/      shared figure output
│
├── synthesized/
│   ├── cbm_theory/            balanced_binary_threshold · decay_cbm_features ·
│   │                          flat_cbm_variance_recovery · nonbalanced_flat_cbm
│   ├── sampling_methods/      nnm_vs_ipw
│   └── tree_reconstruction/   stdr_partition_recovery
│
├── generated/
│   ├── eta_pool_sweeps/       eta_pool_sweep · fiedler_tree_partition_by_eta ·
│   │                          kingman_threshold_vs_theory
│   ├── sampling_methods/      distance_vs_similarity
│   ├── tree_reconstruction/   snj_subsampling
│   ├── nj_distance/           nj_subsampling_balanced · nj_subsampling_nonbalanced
│   └── distance_vs_similarity/ simulation_distance_vs_similarity
│
└── real/
    └── distance_vs_similarity/ real_data_distance_vs_similarity   ← the only real-data notebook
```

## Index

### synthesized — `S`/`D` from a block model
| Notebook | What it produces |
|---|---|
| `cbm_theory/balanced_binary_threshold` | p\* = C·log(n)/n on a symmetric binary tree (α=0.9) |
| `cbm_theory/decay_cbm_features` | HBM (decay-CBM) identities: μ(U), λ₂, Lemma 0.4 slack across η, m |
| `cbm_theory/flat_cbm_variance_recovery` | Flat-CBM cross-clan variance and sign recovery |
| `cbm_theory/nonbalanced_flat_cbm` | p\* ∝ η(1+η)³ log(m) / [m(ρ − η·S_out)²] |
| `sampling_methods/nnm_vs_ipw` | Soft-Impute (NNM) vs IPW: p\* scaling and runtime |
| `tree_reconstruction/stdr_partition_recovery` | STDR recursive partition: first-layer Fiedler, bipartition Jaccard, ARI |

### generated — dendropy tree + sequences
| Notebook | What it produces |
|---|---|
| `eta_pool_sweeps/eta_pool_sweep` | 3 operators (sign / σ₂-gap / k-means L_sym) on the η-pooled Kingman samples. **Fig 1** p\* vs n; **Fig 2** recovery grid + p\* summary. `METRIC` (nmi/ari/agreement) and `AGG` (mean/median) toggles replace the former per_n / per_n_lsym / per_n_median notebooks |
| `eta_pool_sweeps/fiedler_tree_partition_by_eta` | Fiedler vector & tree partition on cached Kingman samples by η; coherence vs imbalance; validates the two η lower bounds |
| `eta_pool_sweeps/kingman_threshold_vs_theory` | Empirical p\* on JC69 + Kingman vs non-balanced CBM theory |
| `sampling_methods/distance_vs_similarity` | Sub-sampling `S` directly vs sub-sampling `D` then `S = exp(−α·D̂)` |
| `tree_reconstruction/snj_subsampling` | SNJ ‖R−R̂‖₂, σ₂ separation, RF distance vs p |
| `nj_distance/nj_subsampling_balanced` | NJ ‖D−D̂‖₂, Q-criterion, RF distance vs p on balanced binary |
| `nj_distance/nj_subsampling_nonbalanced` | NJ on unbalanced birth-death trees |
| `distance_vs_similarity/simulation_distance_vs_similarity` | Generated 600-tree (Kingman + birth-death) recovery: Fiedler-on-S vs Griffing-on-D vs L_sym; screens, η histograms, recovery curves, p\* vs η / diameter |

### real — FASTA + Newick
| Notebook | What it produces |
|---|---|
| `distance_vs_similarity/real_data_distance_vs_similarity` | Real 600-tree benchmark: same operator comparison + validity gate + η histograms + recovery curves. The simulation twin above is the generated analogue |

## Shared utilities (`utils/`, plus one in `src/utils/`)
`utils/` is imported by every notebook via `from analysis.theoretical_interpretation.utils import …`
(notebooks add the project root to `sys.path` with a walk-up to `setup.py`, so this resolves from
any folder depth). **Do not move `utils/`.**

- `block_model`, `balanced_binary` — closed-form `S` and population Fiedler vectors
- `linalg_features` — μ(U), λ₂, spectral gap, Lemma 0.4 row
- `tree_features` — imbalance η, structural margin ρ, n_min, HBM helpers
- `distance_features` — Griffing operators, tree diameter
- `spectral`, `recovery`, `sweep`, `cache`, `plotting` — Fiedler/IPW/NNM, sign-agreement, p\* detector, disk cache, recovery figures
- `generated_data` — `make_generated`, `build_ids` (Kingman + birth-death loader for the generated notebooks)
- **`sweep_plots`** *(new)* — `plot_pstar_vs_n`, `plot_nmi_grid`, `plot_curves_per_n/per_eta`, `compute_pstar` for the eta-pool figures
- **`distance_similarity`** *(new)* — the shared screen/recovery/plot functions behind the two `distance_vs_similarity` notebooks (inject a loader + cache prefix)
- **`src/utils/eta_pool_sweep`** *(new)* — `discover_ns_and_samples`, `run_and_collect_sweeps`: the single source of truth for the "discover pool samples → bootstrap p-sweep → collect" loop (also used by `scripts/build_sweeps.py`)
