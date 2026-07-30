# RUNBOOK — rebuilding the paper from nothing

Per-figure recovery instructions. The authoritative figure→notebook→cache mapping is
`analysis/PAPER_MAP.md`, generated from `paper_figures.py`;
this file adds the *operational* detail: cost, prerequisites, and traps.

## Rule zero: never let a trial run overwrite a paper asset

The paper notebooks save **straight into `docs/overleafs/v9/figures/`**. When exploring,
redirect the output first:

- most notebooks: set `FIG_DIR` to a scratch path
- `fiedler_tree_partition_by_eta`: set `PAPER_FIG_DIR` (its `FIG_DIR` is already scratch)

To preview a figure without republishing the paper asset, exec the cell headlessly with
`FIG_DIR` diverted, then write the PNG back into the `.ipynb` as a base64 `display_data`
output. Strip `%` magics before `exec`-ing notebook cells from a script.

Editing a notebook cell **clears its stored outputs**, and the 6 paper notebooks carry
the published numbers in those outputs. Prefer patching the JSON `source` arrays over a
round-trip that drops them.

## Check before you build

```bash
python scripts/sync_paper_figures.py --check    # is anything missing / unmapped / stale?
python scripts/sync_paper_figures.py            # md5 + mtime + producer per figure
```

`--check`'s STALE and MIXED-RUN lines are advisory mtime heuristics; a fresh clone gives
every file the same checkout time, so they say nothing there. The exact checks (MISSING,
UNMAPPED, UNUSED-MAP, ORPHAN, NO-NOTEBOOK) are the ones that gate.

## Cost model — extrapolate before launching

Cost per `(η, m)` cell is `N_TRIALS × |P_GRID|` sub-sampled Fiedler solves, each `O(m³)`
once `p > 0.1` (dense `scipy.linalg.eigh`; sparse `eigsh` below that).

| m | wall clock per cell |
|---|---|
| 1536 | ~39 s |
| 2816 | ~3 min |
| ~11000 | ~3 h |

Prefer fewer sizes or fewer trials over an overnight tail: the largest point buys one dot
on a log axis. Run long jobs in the background, smallest sizes first, with progress to a
log — the figure is assemblable before the tail finishes.

**Check what is already cached first.** A complete flat-CBM cell is 480 trial dirs
(24 four-decimal `p` keys × 20 seeds). See `docs/CACHE_AND_RESULTS.md` for current
contents.

## Build the paper PDF

```bash
cd docs/overleafs/v9
latexmk -pdf -interaction=nonstopmode thesis_v9.tex
grep -c "undefined" thesis_v9.log      # expected: 0
```

Expected output: **29 pages, 2,586,870 bytes, 0 undefined**. Run `latexmk` from *inside*
`v9/` — the shell cwd is not where you think after backgrounded commands. A truncated
`.aux`/`.out` from an interrupted run causes
`File ended while scanning use of \@newl@bel`; fix with `latexmk -C`, then rebuild.

`v9/` is self-contained: do **not** re-add `{../}` to `\graphicspath`. It previously let
Figure 5 resolve outside the tree, so uploading `v9/` alone to Overleaf failed.

---

## Figure 2 — synthesized flat CBM, `fig:pstar_synth`

Notebook: `paper/fig02_pstar_synth_cbm.ipynb`

No sweep prerequisite — `S` is closed-form, so the notebook is self-sufficient (it
caches trials under `cache/full_matrix` + `cache/sweep_trial` via `utils/cache.py`).

**Trap:** `m` must be divisible by `1+η` for `n₁ = m/(1+η)` to be exact, so **each η
needs its own m grid** (÷2, ÷6, ÷11, ÷16 for η = 1, 5, 10, 15). Reusing one grid
silently shifts the realised η — m=90 at η=15 gives 17.0.

## Figures 3, 8, 9 — generated Kingman, `fig:pstar_gen` / `recovery_grid` / `operator_sensitivity`

Notebook: `paper/fig03_pstar_gen_kingman.ipynb` (one notebook, three figures)

```bash
python scripts/build_eta_pool.py --n 500          # -> cache/pool_sample   (or build_eta_pool_parallel.py)
python scripts/build_sweeps.py --ns 500 --methods kmeans   # -> cache/bootstrap_sweep
```

Then run the notebook. For **Figure 9** you additionally need the other two operators:

```bash
python scripts/build_sweeps.py --ns 500 --methods sign sigma2 kmeans
```

`sigma2` is the expensive operator — pass it via `cached_only_methods` so a newly added
`n` cannot trigger a fresh expensive run.

## Figure 4 — HBM spectral gap, `fig:hbm_spectral_verification`

Notebook: `paper/fig04_hbm_spectral_gap.ipynb`

No prerequisite; HBM matrices are built in memory. **Both panels must come from one
execution** — a v8 version of this figure shipped with two panels from two different runs
(`open-items/00-R1.md [R1/23]`), which is why `--check` has a MIXED-RUN test.

## Figure 5 — Fiedler bipartitions by η, `fig:fiedler_partitions`

Notebook: `paper/fig05_fiedler_partitions.ipynb`

```bash
python scripts/build_eta_pool.py --n 500     # -> cache/pool_sample
```

Writes the four panels to `v9/figures/Fiedler_Bipartitions/` via `PAPER_FIG_DIR`. Its two
coherence panels go to the scratch `analysis/.../figures/` dir and are not used by the
paper.

## Figure 6 — balanced binary threshold, `fig:pstar_balanced`

Notebook: `paper/fig06_pstar_balanced.ipynb`

Reads `trials.csv` / `agg.csv` from
`results/notebooks/01_cbm_theory/balanced_binary_threshold/`. That `01_cbm_theory`
prefix is a **historical cache key, not a path** — the directory of that name no longer
exists in `analysis/`. Renaming the key orphans the 501 KB of cached trials behind this
figure.

The caption still prints `C=[TBD, C3]`; see `open-items/19-cleanup.md [CLEANUP/01]`.
When you do fill it in: `C` is the **median of the per-point ratios**
(`median_ratio_C`), never a linear-scale least-squares fit, and report the max/min spread
when it is large.

## Figure 7 — identity checks, `fig:identity_gen`

Notebook: `paper/fig07_identity_checks.ipynb`

**This is the one figure not reproducible from a clean checkout.** The notebook asserts,
in its *setup* cell, that this exists:

```
results/runs/kingman_mean/uniform/20260501-194101-kingman_mean_n500-8000_mu_0p1_uniform/
    n500_L10000/results.json
```

`results/` is gitignored, so in a fresh clone the notebook fails before any compute. To
regenerate, run **pipeline A** with `tree_model="kingman_mean"`,
`sampling_method="uniform"`, `taxa_values=[500…8000]`, `mutation_rate=0.1`,
`sequence_length=10000` — edit `SWEEP_CONFIG` at the top of
`scripts/run_experiment.py` (that pipeline takes no CLI arguments), or use
`scripts/interactive_run.py`. This is the most expensive item in the runbook.

Tracked as `open-items/19-cleanup.md [CLEANUP/08]`.

---

## Regenerating the open-items register

```bash
python scripts/collate_open_items.py     # open-items/*.md -> OPEN_ITEMS.md
```

Edit the fragments, never the collated file. Entry format is strict:
`### [ID/NN] title`, then `- **Anchor:**`, `- **Type:**`, `- **Item:**` (≤3 sentences),
`- **Options:**` (iff Type is `dilemma`), `- **Blocking:**` (must start yes/no). The
collator lints and will report entries it cannot parse as `0 entries` for that file.
