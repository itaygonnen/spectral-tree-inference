# analysis

Notebooks behind the v9 manuscript (`docs/overleafs/v9/`). **The paper is the deliverable**;
these notebooks exist to produce and defend its figures.

Four things, nothing nested inside them:

| | What |
|---|---|
| `paper/` | The 6 notebooks that produce a figure in the manuscript. Named `figNN_*` so `ls` answers "which notebook makes Figure 3?". |
| `supporting/` | The 12 that produce **no** paper figure. Flat, one file each. |
| `utils/` | Shared code, imported as `analysis.utils`. **Do not move.** |
| `notebooks_cache/` | Every `.npz` a notebook computes. Tracked, so re-plotting works in a fresh clone without re-running an expensive screen. |

Plus `legacy/` (superseded packages, unmaintained) and `figures/` (gitignored scratch —
paper figures go to `docs/overleafs/v9/figures/`, never here).

## paper/

| Notebook | Produces | Section |
|---|---|---|
| `fig02_pstar_synth_cbm.ipynb` | **Fig 2** (a–d) | §5, `empirical.tex` |
| `fig03_pstar_gen_kingman.ipynb` | **Fig 3** (a–d), **and Figs 8 + 9** | §5 + App G |
| `fig04_hbm_spectral_gap.ipynb` | **Fig 4** | App D |
| `fig05_fiedler_partitions.ipynb` | **Fig 5** (a–d) | App D |
| `fig06_pstar_balanced.ipynb` | **Fig 6** | App G |
| `fig07_identity_checks.ipynb` | **Fig 7** | App G |

`fig03` is the one notebook spanning two sections — Fig 3 in the main text plus Figs 8 and
9 in the appendix. The `figNN_` prefix names its primary figure only;
**[PAPER_MAP.md](PAPER_MAP.md)** is the complete mapping.

## supporting/

12 notebooks, no paper figure. Which claim each backs — or "exploratory, not cited" — is in
[PAPER_MAP.md](PAPER_MAP.md). Four of them (`distance_vs_similarity_*`, `pstar_*_distance`)
are the empirical evidence for the **distance route**, App F / `thm:main-dist`, which
carries no figure of its own.

## The authority on provenance

**[PAPER_MAP.md](PAPER_MAP.md)** maps every figure to its notebook, the claim it supports,
the cache it consumes, and the command that rebuilds it. It is **generated** from
`paper_figures.py`:

```bash
python scripts/sync_paper_figures.py --check      # validate; non-zero on problems
python scripts/sync_paper_figures.py --write-map  # regenerate PAPER_MAP.md
```

`--check` is not decorative. It fails if a figure is missing, if the paper references a
figure the manifest does not know about, if a mapped figure is unreferenced, or if a
notebook exists that is listed neither as a figure producer nor as supporting. Add a
notebook and the check tells you to declare it.

## Traps worth knowing before you re-run anything

- **A trial run can overwrite a paper asset.** The `paper/` notebooks save straight into
  `docs/overleafs/v9/figures/`. When exploring, point `FIG_DIR` / `PAPER_FIG_DIR` at a
  scratch path first.
- **A cell that looks like setup may launch a sweep.** `nj_subsampling_nonbalanced` cell 15
  bootstraps `sys.path` *and* kicks off an n=4000 synthetic sweep.
- **Some `notebook_dir(...)` arguments are historical cache keys, not paths** — e.g.
  `"01_cbm_theory/balanced_binary_threshold"` in `fig06`. Those directories are long gone;
  the keys are kept so cached results still resolve. Renaming one orphans real data.
- **`fig07` needs a `results/` run that is gitignored**, so it cannot run in a fresh clone
  until that sweep is regenerated. See `docs/RUNBOOK.md`.
- Editing a notebook cell clears its stored outputs, and the `paper/` notebooks carry the
  published numbers there. Patch the JSON `source` array rather than round-tripping.

## utils/

Imported as `from analysis.utils import …`. Notebooks put the package root on `sys.path` via
a depth-invariant walk up to `setup.py`, so this resolves from any depth.

- `block_model`, `balanced_binary` — closed-form `S` and population Fiedler vectors
- `linalg_features` — μ(U), λ₂, spectral gap, Lemma 0.4 row
- `tree_features` — imbalance η, structural margin ρ, `n_min`, HBM helpers
- `spectral`, `recovery`, `sweep`, `cache`, `plotting` — Fiedler/IPW/NNM, sign-agreement,
  ARI/NMI, the p\* detector, disk cache, recovery figures
- `perturbation` — Davis–Kahan / sin-Θ machinery
- `generated_data` — Kingman + birth-death loader
- `sweep_plots` — η-pool plots for the **appendix** figures (fits `C` by least squares)
- `sweep_plots_two_panel` — the per-η 2-panel figure for the **main** figures (fits `C` by
  median-of-ratios). The two are **not** interchangeable; see their headers.
- `src/utils/threshold_utils` — `find_discrete_threshold` (the p̂\* read-off), `fit_power_law`
