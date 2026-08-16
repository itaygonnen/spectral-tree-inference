# analysis

Notebooks behind the v9 manuscript (`docs/overleafs/v9/`). **The paper is the deliverable**;
these notebooks exist to produce and defend its figures.

Four things, nothing nested inside them:

| | What |
|---|---|
| `paper/` | The notebooks that produce a figure in the manuscript. Named `figNN_*` so `ls` answers "which notebook makes Figure 3?" — except the two distance-route producers, which still sit in `supporting/` (see below). |
| `supporting/` | 11 files: 8 that produce **no** paper figure, plus 3 figure producers that were promoted without moving (see below). Flat, one file each. |
| `utils/` | Shared code, imported as `analysis.utils`. **Do not move.** |
| `notebooks_cache/` | Every `.npz` a notebook computes. Tracked, so re-plotting works in a fresh clone without re-running an expensive screen. |

Plus `legacy/` (superseded packages, unmaintained) and `figures/` (gitignored scratch —
paper figures go to `docs/overleafs/v9/figures/`, never here).

## paper/

| Notebook | Produces | Section |
|---|---|---|
| `fig02_pstar_synth_cbm.ipynb` | **Fig 2** (a–d) | §5, `empirical.tex` |
| `fig03_pstar_gen_kingman.ipynb` | **Fig 3** (a–d) | §5, `empirical.tex` |
| `fig04_hbm_spectral_gap.ipynb` | **Fig 6** | App D |
| `fig05_fiedler_partitions.ipynb` | **Fig 7** (a–d) | App D |

**The `figNN_` prefix is the filename, not the printed number, and the two have drifted
apart.** Figures 4 and 5 are the distance-route pair in §5, so the App D figures shifted up
by two: `fig04_*` prints as Figure 6 and `fig05_*` as 7. Renaming the files to match would
only move the problem the next time a figure is inserted;
**[PAPER_MAP.md](PAPER_MAP.md)** is the authority, and `sync_paper_figures.py --check` is
what enforces it.

**2026-08-04:** `sections/appendix-emp.tex` (compiled App G, supplementary empirical
figures) was deleted, taking Figs 8–11 with it. Two notebooks left `paper/` as a result —
`fig06_pstar_balanced` → `supporting/pstar_balanced_nmi`, `fig07_identity_checks` →
`supporting/identity_checks_kingman` — and `fig03` no longer emits Figs 10/11 into the
paper (the cells still run; the PNGs are in `RETIRED`). `sections/appendix-H.tex` (the
S-vs-D comparison) was deleted the same day, taking `fig08_appH_estimator` and
`supporting/appendix_H_estimator_vs_operator` with it, so App **F** (`appendix-G.tex`) is
now the last appendix.

`fig03` is the one notebook spanning two sections — Fig 3 in the main text plus Figs 10 and
11 in the appendix.

## supporting/

8 notebooks with no paper figure. Which claim each backs — or "exploratory, not cited" — is
in [PAPER_MAP.md](PAPER_MAP.md). Two of them (`distance_vs_similarity_*`) remain the
non-figure evidence for the **distance route**, App F / `thm:main-dist`.

Those two run the **same pipeline on different data** (real 600-tree FASTA benchmark vs its
generated twin), so their shared spine lives in `src/`, not copied into both: screens in
`src/utils/screening.py`, the operator × threshold sweep and per-tree helpers in
`src/utils/operator_comparison.py`, the two comparison panels in
`src/plots/plot_operator_threshold.py`. Change a panel there and both notebooks stay
comparable — which is the only reason keeping both is worth anything. The operator and
threshold callables are passed **in** from the notebooks, because the σ₂ rule is
`spectraltree.partition_taxa` and `src/` must not import `spectraltree`.

Three files here **do** produce paper figures and are registered in `paper_figures.py`:
`eta_by_operator.ipynb` (Fig 8 a-d, `appendix-eta.tex`), and the two distance-route p* notebooks —
`pstar_flat_cbm_distance.ipynb` and `pstar_eta_pool_distance.ipynb` (Figs 4 and 5). They
still live here rather than in `paper/`: they
were written as supporting evidence and promoted afterwards. Move and rename them to
`paper/fig04_*` / `paper/fig05_*` only together with updating their `notebook=` fields in
`paper_figures.py`, or `--check` will fail with `NO-NOTEBOOK`.

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
- **A cell that looks like setup may launch a sweep.** `pstar_flat_cbm_distance` cell 2 and
  `pstar_eta_pool_distance` cell 2 bootstrap `sys.path` *and* kick off the grid.
- **Some `notebook_dir(...)` arguments are historical cache keys, not paths** — e.g.
  `"01_cbm_theory/balanced_binary_threshold"` in `pstar_balanced_nmi`. Those directories are
  long gone; the keys are kept so cached results still resolve. Renaming one orphans real data.
- **`identity_checks_kingman` needs a `results/` run that is gitignored**, so it cannot run in
  a fresh clone until that sweep is regenerated. See `docs/RUNBOOK.md`.
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
