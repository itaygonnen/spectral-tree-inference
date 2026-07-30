# analysis

Notebooks behind the v9 manuscript (`docs/overleafs/v9/`). **The paper is the deliverable**;
these notebooks exist to produce and defend its figures.

Organized **by paper section, then by data source**. The section tells you what claim the
notebook serves; the source tells you where its numbers come from. Only 6 of the 16
notebooks produce a paper figure — the rest are labelled `supporting/` so nobody mistakes
scratch for a deliverable.

For the authoritative figure ↔ notebook ↔ claim ↔ cache ↔ command mapping, see
**[PAPER_MAP.md](PAPER_MAP.md)**. It is machine-checked by
`python scripts/sync_paper_figures.py --check`, which fails if a figure is missing, stale,
built from mixed runs, or absent from the map.

## Data sources

| Source | Meaning |
|---|---|
| **synthesized** | `S`/`D` built **directly** from a closed-form block model (CBM/HBM, balanced binary) — no tree simulation |
| **generated** | **dendropy** simulates a tree (Kingman / birth-death) and evolves sequences (JC69); `S`/`D` come from those sequences |
| **real** | downloaded **FASTA + Newick** (the 600-tree, 1000-taxon benchmark) |

## Layout

```
analysis/
├── PAPER_MAP.md   figure ↔ notebook ↔ claim ↔ cache ↔ rebuild command
├── utils/         shared helper library — imported by every notebook, do NOT move
├── figures/       scratch figure output (gitignored; NOT where paper figures live)
│
├── sec5_empirical/           §5 Empirical Results  (sections/empirical.tex)
│   ├── synthesized/  nonbalanced_flat_cbm            → Fig 2 (a–d)
│   └── generated/    eta_pool_sweep                  → Fig 3 (a–d), and Figs 8, 9
│
├── appD_hbm/                 App D: HBM extension  (sections/appendix-D.tex)
│   ├── synthesized/  hbm_spectral_gap_verification   → Fig 4
│   └── generated/    fiedler_tree_partition_by_eta   → Fig 5 (a–d)
│
├── appG_supplementary/       App G: Supplementary Figures  (sections/appendix-emp.tex)
│   ├── synthesized/  balanced_binary_threshold       → Fig 6
│   └── generated/    kingman_threshold_vs_theory     → Fig 7
│
└── supporting/               no paper figure — see PAPER_MAP.md for what each backs
    ├── synthesized/  cbm_theory/{flat_cbm_variance_recovery, decay_cbm_features}
    │                 sampling_methods/nnm_vs_ipw
    │                 tree_reconstruction/stdr_partition_recovery
    ├── generated/    distance_vs_similarity/simulation_distance_vs_similarity
    │                 nj_distance/{nj_subsampling_balanced, nj_subsampling_nonbalanced}
    │                 tree_reconstruction/snj_subsampling
    │                 sampling_methods/distance_vs_similarity
    └── real/         distance_vs_similarity/real_data_distance_vs_similarity
```

Two things the layout cannot express, so they are stated instead:

- **`eta_pool_sweep` spans two sections.** It produces Fig 3 (§5) *and* Figs 8 and 9
  (App G). It is filed under its primary figure; `PAPER_MAP.md` lists all three.
- **Directory letters follow the compiled PDF, not the filenames.** `appendix-emp.tex`
  compiles as Appendix **G** and `appendix-G.tex` as Appendix **F**. That mismatch is
  recorded for the author in `docs/overleafs/v9/open-items/19-cleanup.md`.

## Traps worth knowing before you re-run anything

- **A trial run can overwrite a paper asset.** The paper notebooks save straight into
  `docs/overleafs/v9/figures/`. When exploring, redirect `FIG_DIR` / `PAPER_FIG_DIR` to a
  scratch path first.
- **A cell that looks like setup may launch a sweep.** `nj_subsampling_nonbalanced` cell 15
  bootstraps `sys.path` *and* kicks off an n=4000 synthetic sweep. Do not assume the leading
  cells are cheap.
- **Some `notebook_dir(...)` arguments are historical cache keys**, not paths — e.g.
  `"01_cbm_theory/balanced_binary_threshold"` and `"03_sampling_methods/nnm_vs_ipw"`. The
  directories of those names are gone; the keys are kept so cached results still resolve.
  Renaming them orphans real data.
- **`kingman_threshold_vs_theory` (Fig 7) needs a `results/` run that is gitignored**, so it
  cannot run in a fresh clone until that sweep is regenerated. See `docs/RUNBOOK.md`.
- Editing a notebook cell clears its stored outputs. The 6 paper notebooks carry the
  published numbers in those outputs — prefer patching JSON `source` arrays over a
  round-trip that drops them.

## Shared utilities (`utils/`)

Imported as `from analysis.utils import …`. Notebooks put the
package root on `sys.path` via a depth-invariant walk up to `setup.py`, so this resolves
from any folder depth. **Do not move `utils/`.**

- `block_model`, `balanced_binary` — closed-form `S` and population Fiedler vectors
- `linalg_features` — μ(U), λ₂, spectral gap, Lemma 0.4 row
- `tree_features` — imbalance η, structural margin ρ, `n_min`, HBM helpers
- `spectral`, `recovery`, `sweep`, `cache`, `plotting` — Fiedler/IPW/NNM, sign-agreement,
  ARI/NMI, the p\* detector, disk cache, recovery figures
- `perturbation` — Davis–Kahan / sin-Θ machinery
- `generated_data` — `make_generated`, `build_ids` (Kingman + birth-death loader)
- `sweep_plots` — multi-figure η-pool plots (`plot_pstar_vs_n`, `plot_nmi_grid`, …)
- `sweep_plots_two_panel` — the per-η 2-panel figure tiled by LaTeX `subfigure`
- `src/utils/eta_pool_sweep` — single source of truth for the
  "discover pool samples → bootstrap p-sweep → collect" loop, shared with
  `scripts/build_sweeps.py`
- `src/utils/threshold_utils` — `find_discrete_threshold` (the p̂\* read-off),
  `fit_power_law`
