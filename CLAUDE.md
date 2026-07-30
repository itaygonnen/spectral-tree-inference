# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository layout (two distinct projects in one repo)

This repo contains two related but separate codebases. Don't mix imports between them.

1. **`spectraltree/`** — installable Python package (`pip install -e .`) implementing phylogenetic tree reconstruction methods (STDR, SNJ, NJ, RAxML wrapper, etc.). Stable, paper-backing code. Tests live in `tests/`.
2. **`sub_sampled_fielder_vec/`** — active research framework for analyzing Fiedler-vector quality under matrix sub-sampling. **Not pip-installed**; scripts add the package root to `sys.path` and import as `from src.<module> ...`. This is where most current development happens.

## Common commands

### `spectraltree` package
```bash
pip install -e .                              # install once
python -m unittest discover                   # run full test suite
python -m unittest tests.test_snj             # run a single module
python -m unittest tests.test_snj.TestSNJ.test_xxx   # run a single test
python setup.py build_sphinx                  # build docs (output in docs/_build/html)
```
Some third-party method tests (Forrest, RG, RAxML) are flaky — they should PASS or FAIL but never ERROR. RAxML and Octave-based tests need external binaries (`octave` via `conda install -c conda-forge octave`).

### `sub_sampled_fielder_vec` experiments
```bash
cd sub_sampled_fielder_vec
python scripts/interactive_run.py             # menu-driven launcher (recommended)
python scripts/run_experiment.py              # edit SWEEP_CONFIG dict at top, then run
```
Results land in `sub_sampled_fielder_vec/results/<timestamp>-<run_name>/n{taxa}_L{seq_len}/`. Matrix caches live in `sub_sampled_fielder_vec/src/cache/` (gitignored, reused across runs keyed by `(n_taxa, seq_len, mu, tree_model, seq_model)`).

There is no test suite for `sub_sampled_fielder_vec/` — see `tests/README.md` for why (a bare
`test*.py` in `.gitignore` made any test file uncommittable; now negated for that dir). Validate
changes with:

```bash
python -m src.cache_io                        # cache key/sentinel smoke checks
python scripts/sync_paper_figures.py --check  # figure provenance; non-zero on problems
python scripts/collate_open_items.py          # regenerates + lints the open-items register
```

`collate_open_items.py` currently **exits 1** on one pre-existing lint warning
(`13-B1.md [B1/29]`, Item longer than 3 sentences). It still writes `OPEN_ITEMS.md`
correctly — that non-zero exit is not a regression.

plus a small experiment (e.g. `taxa_values=[256], bootstrap_reps=2, p_values=[0.1, 0.5, 1.0]`).

## Architecture: `spectraltree`

Top-level `spectraltree/__init__.py` re-exports the public API — import via `import spectraltree` and use `spectraltree.balanced_binary(...)`, `spectraltree.Jukes_Cantor()`, etc. Inside the package, use **relative imports** (`from . import utils`).

Reconstruction pipeline:
1. `generation.py` simulates sequences on a tree (Dendropy-backed; `Jukes_Cantor`, `HKY`, `GTR`, ...).
2. `similarities.py` builds a similarity/distance matrix from observations (`paralinear_distance`, `JC_similarity_matrix`, `HKY_similarity_matrix(obs)(metaHKY)`).
3. A `ReconstructionMethod` subclass consumes that matrix:
   - `reconstruct_tree.py`: `NeighborJoining`, `TreeSVD`, base classes.
   - `snj.py`: `SpectralNeighborJoining` (uses 2nd singular value criterion).
   - `spectral_tree_reconstruction.py`: `STDR` (top-down spectral recovery).
   - `recursive_str.py`: `STR`.
   - `raxml_reconstruction.py`: external RAxML wrapper.
   - `choi_reconstruction.py`, `forrest_reconstruction.py`: legacy/optional, currently commented out of `__init__.py`.

Interface convention (post-2021 refactor): pass method-specific args to the constructor; call the instance with only `(observations, taxon_metadata)`. `spectral_tree_reconstruction` is the lone exception — leave its interface alone.

## Architecture: `sub_sampled_fielder_vec`

This framework asks: *if we replace the full similarity matrix `M` with a sub-sampled `S` (sampling rate `p`), how well does the Fiedler vector of `L(S)` recover the partition that `L(M)` would produce?*

### Data flow (one experiment)
`scripts/run_experiment.py` → `extract_config_values` → `ExperimentRunner.run()` → for each `(n_taxa, seq_len)` calls `sweep_for_params(cfg, n, L, run_dir)` in `src/runners/bootstrap_sweep.py`. That function:
1. Generates sequences on a tree (`src/models/`).
2. Builds the full similarity matrix `M` and computes the **reference Fiedler vector** at `p=1.0`.
3. For each `p` in `cfg.experiment.p_values`: runs a bootstrap loop of `K` reps, each time sub-sampling `M → S_k` (uniform / leveraged / LDS), maintaining a streaming average `S̄` via Welford, computing the Fiedler of `L(S_k)`, sign-aligning to the reference, and accumulating.
4. After the loop, computes partition-agreement metrics (`partition_agreement_M`, `partition_agreement_S`), `dot_product`, `sign_agreement` in `src/utils/metrics.py` and writes `results.json` + plots.

### Configuration system
Pydantic v2 models in `src/config/base_config.py`: `StructuredConfig` bundles `TreeConfig`, `SequenceConfig`, `ExperimentConfig`, `SamplingConfig`, `MetricsConfig`, `GuardrailsConfig`, `CacheConfig`, `OutputConfig`. Build via `src/config/presets.py:custom_config(...)`. The script-level `SWEEP_CONFIG` dict is flattened into these by `extract_config_values()`.

### Sampling methods (`src/core/sampling/`)
Pluggable via `SamplingConfig.method`:
- **`uniform`** — baseline; sample `O(p·n²)` upper-triangular entries uniformly.
- **`leveraged`** (IALM) — Phase 1 uniform → estimate leverage scores → Phase 2 leverage-biased sampling → IALM matrix completion. Highest accuracy, slow. `src/core/sampling/leveraged/`.
- **`lds`** — Leveraged Debiased Sampler: same two-phase sampling but uses an unbiased IPW estimator with effective inclusion probability `π_ij = p_0 + (1-p_0)·p_ij` instead of IALM completion. 10–100× faster than IALM. **Critical**: debias with `π_ij`, not `p_ij` alone — see `CHANGELOG.md` 2026-02-20.

`SamplingConfig` knobs that control behavior at low `p`:
- `allow_uniform_fallback=True` (safe mode) — falls back to uniform when leveraged budget is below theoretical minimum.
- `force_leveraged=True` (research mode) — proceeds with leveraged using a 90/10 phase split, no guardrails.
- `ialm_bypass_threshold` — skip IALM and return the sparse sample directly when `p` is high enough.

### Performance constraints (load-bearing)
The November-2025 optimizations in `src/core/` are why experiments at `n=8192` are tractable. Don't regress them:
- Fiedler computation uses partial eigendecomposition (`scipy.linalg.eigh` exploiting symmetry), not full SVD.
- Sub-sampling uses vectorized boolean masks, not `np.random.choice` on `n²` elements.
- Laplacians are computed once per `p` (when needed for metrics), not per bootstrap iteration.
- M-based metrics (operator norm, coherence) are cached before the `p`-loop.

### Cross-project boundary
`sub_sampled_fielder_vec` does **not** import from `spectraltree`. Tree generation lives in `src/models/tree_models.py`; sequence simulation in `src/models/sequence_models.py`. If you need something from `spectraltree`, port it — don't add a cross-import.

### Analysis notebooks

`analysis/` is flat and holds exactly four things:

- **`paper/`** — the 6 notebooks that produce a manuscript figure, named `figNN_*` so `ls`
  answers "which notebook makes Figure 3?". `fig03` also emits Figs 8 and 9.
- **`supporting/`** — the 12 that produce no paper figure. Flat, one file each.
- **`utils/`** — shared code, imported as `analysis.utils`. **Do not move it.**
- **`notebooks_cache/`** — every `.npz` a notebook computes, tracked so re-plotting works
  in a fresh clone. Data does NOT live beside notebooks.

There is no per-topic nesting: 18 notebooks previously sat in 16 directories, which
grouped nothing.

Its one sibling is **`analysis/legacy/`** — six superseded packages (`comparison`,
`generic_analysis`, `leveraged_sampling_analysis`, `notebooks`, `scripts`,
`spectral_analysis`). Unmaintained; don't add to them, and don't revive an import from
them.

**`PAPER_MAP.md` is the authority** for which notebook produces which figure, which
claim it supports, which cache it consumes, and how to rebuild it. It is *generated*
from `paper_figures.py` — edit that, then run
`python scripts/sync_paper_figures.py --write-map`. Validate with `--check`, which
cross-checks the manifest against every `\includegraphics` in `v9/sections/*.tex` and
exits non-zero on a mismatch.

Load-bearing invariants:

- **`src/` must never import from `analysis/`.** Two modules used to be pulled in by
  bare name through `sys.path` hacks (`partition_validity`, `tree_plots`), which meant
  `import src` only worked by accident. Both now live in `src/utils/`.
- Some `notebook_dir(...)` arguments are **historical cache keys, not paths** —
  `"01_cbm_theory/balanced_binary_threshold"` and `"03_sampling_methods/nnm_vs_ipw"`.
  The directories are gone; renaming the keys orphans real cached results.
- A cell that looks like setup may launch a sweep (`nj_subsampling_nonbalanced` cell 15
  starts an n=4000 run).
- Editing a notebook cell clears its stored outputs, and the 6 paper notebooks hold the
  published numbers there. Patch the JSON `source` arrays instead of round-tripping.

Key docs: `docs/RUNBOOK.md` (rebuild figure N, with costs), `docs/CACHE_AND_RESULTS.md`
(all 7 cache scopes, what's orphaned), `scripts/README.md` (which `plot_*` files are
libraries that notebooks import, so must not be moved).

## Thesis presentation (`sub_sampled_fielder_vec/docs/thesis_seminar.html`)

Single self-contained HTML file — all CSS, JS, and canvas drawing code inline. MathJax loaded from CDN (requires internet to render LaTeX). Open with `open sub_sampled_fielder_vec/docs/thesis_seminar.html`.

### Slide structure (21 slides)
| Part | Slides | Topic |
|------|--------|-------|
| I | 1–6 | Domain intro, STDR |
| II | 7–9 | The sub-sampling problem |
| III | 10–13 | Building S (CBM) |
| IV | 14–16 | Bridging tree topology to linear algebra |
| V | 17–21 | Empirical results |

### FIG registry pattern
All canvas figures live in `const FIG = { key: cv => drawFn(cv._x, cv._w, cv._h, opts) }`. To add a new figure: add an entry here, then reference it with `<canvas data-fig="key" width="W" height="H"></canvas>` in a slide. The deck engine calls `hidpi(cv)` and then `FIG[key](cv)` for every canvas on the current slide.

**Currently defined but not yet placed on any slide:** `phase`, `bern`, `compare`.

### `lineChart` options
```js
lineChart(ctx, W, H, {
  xData,          // array of p values (log x-axis)
  curves,         // [{label, color, y:[]}]
  yRange,         // default [0,1] — e.g. [0.4,1] to zoom
  threshold,      // red dashed horizontal line (e.g. 0.95)
  title,          // top label
  yLabel,         // rotated y-axis label, default 'NMI'
  compact,        // true → tighter padding (for small multi-panel layouts)
  noYAxis,        // true → hide y-axis ticks and label
  noLegend,       // true → hide in-chart legend
})
```

### CSS color semantics — don't cross these
- `--blue` / `--orange` — **clan A / clan B** exclusively. Don't use for neutral UI chrome.
- `--green` — recovery / success regions.
- `--red` — thresholds, split edges, error.
- `--line` / `--ink-dim` — neutral borders, secondary text.

### Python editing trap — LaTeX backslash corruption
When writing HTML content containing MathJax LaTeX from Python strings, `\r`, `\t`, `\n` are interpreted as control characters. `\rho` → CR + `ho`; `\text` → tab + `ext`. Always use **raw strings** (`r"..."`) or **double backslashes** (`\\rho`, `\\text`) when constructing HTML with Python.

## The v9 manuscript (`sub_sampled_fielder_vec/docs/overleafs/v9/`)

Thin master `thesis_v9.tex` + `sections/*.tex`; figures in `figures/`; open questions in
`open-items/` (collated to `OPEN_ITEMS.md`). Related Work is **not** at
`deferred/related.tex` — that dir is empty and the master's claim is stale; the outline is
`sections/related.tex`.

- **v9 is now TRACKED by git** — `.tex`, `figures/`, `open-items/`, and the PDF. It is also
  self-contained: every `\includegraphics` resolves inside `v9/figures/`, so do **not** re-add
  `{../}` to `\graphicspath` (that is what made Figure 5 depend on assets outside the tree).
  Build artifacts are excluded by `v9/.gitignore`.
- **`v7/`, `v8/` and `distance approach/` are still untracked AND gitignored** — for those,
  assume nothing is recoverable and snapshot before destructive edits.
- Appendix filenames do not match compiled letters: `appendix-G.tex` → App **F**,
  `appendix-emp.tex` → App **G**. The notebook dirs follow the *compiled* letters.
- ALWAYS build with `latexmk -pdf -interaction=nonstopmode thesis_v9.tex` **from inside the v9
  directory** (the shell cwd is not where you think after backgrounded commands), and verify with
  `grep -c "undefined" thesis_v9.log` against a baseline taken *before* editing. Zero is the
  expected value. A truncated `.aux`/`.out` from an interrupted run causes
  `File ended while scanning use of \@newl@bel`; fix with `latexmk -C` then rebuild.
- ALWAYS put restructure rationale in the `.tex` as `% [Reviewer Note: ...]` — invisible in the PDF.
- **NEVER attach a bare `\label{}` to unnumbered material** (after `\subsection*`, `\paragraph`, or
  mid-paragraph). It binds to the last stepped counter and `\Cref` then prints a figure or lemma
  number. When a heading or `definition` environment is dissolved, either keep a numbered anchor or
  retarget every reference — `\Cref{def:x}` → `\eqref{eq:x}`, `\Cref{sec:sub}` → the parent section.
- ALWAYS re-grep for orphaned targets after removing structure; appendices reference main-text
  labels heavily (`appendix-emp.tex` alone held 9 references to two dissolved subsections).
- ALWAYS preserve every `\cite{}` key when collapsing bullet lists into prose, and keep the
  claim discipline in `open-items/12-A3.md`: "sufficient", never "tight"/"optimal"; `cor:tolerance`
  is *necessary*; `cor:infeasible` means the guarantee goes silent, not a converse.

### Paper figure conventions (established in the v9 refine pass)

- **NEVER leave matplotlib titles or `suptitle` in a paper figure** — the LaTeX caption carries that
  text. Same for parameter labels: the η (or other sweep parameter) belongs in the `\subcaption`,
  not inside the image.
- ALWAYS emit **one image file per panel-group and tile in LaTeX**, not one giant multi-row PNG.
  The 2×2 pattern: `subfigure[t]{0.49\textwidth}` with
  `\captionsetup{singlelinecheck=false, justification=raggedright, skip=1pt}` and the `\subcaption`
  *before* `\includegraphics` (that is what puts the "(a) η=1" in the upper-left corner). Keep the
  outer `\label` so existing `\Cref`s resolve; add per-cell labels for individual citation.
- ALWAYS size type for the **printed** scale, not the PNG: a 7in-wide figure in a
  `0.49\textwidth` subfigure is scaled ~0.46×, so anything under ~13pt lands under 6pt on paper.
- At subfigure size, NEVER put a multi-entry legend or explanatory text inside the axes — state the
  colour order, marker roles and fitted constants once in the caption. Thin the ticks
  (`NMI_TICKS`, `N_TICKS` in `sweep_plots_two_panel.py`).
- ALWAYS reuse the same plotting function across figures that should look alike
  (`plot_pstar_pair`) and parameterize the differences (`n_ticks`, `size_label`) rather than
  copying the drawing code.
- ALWAYS plot **one** reference curve. Two renderings of the same rate (a fitted `C log n/n` beside
  the theorem's own constant) confuse rather than corroborate.
- Estimating `C` in `p* = C log n / n`: ALWAYS use the **median of the per-point ratios**
  (`median_ratio_C`), never a linear-scale least-squares fit — the latter is set almost entirely by
  the largest `p*`, i.e. by the smallest trees. Report the max/min spread of those per-point ratios
  when it is large; it is the honest statement of how well one constant fits.
- The theorem's own constant `8(1+η)³β₀²log m/[m·margin²]` (`cbm_sufficient_p`) is **not plottable
  on measured data**: rebuilt per size from pool medians it swings an order of magnitude between
  neighbouring `n`, and it is `nan` wherever the measured margin `ρ − η·S_out^max` goes
  non-positive (most sizes at η ≥ 10 in the Kingman pool). Keep it as a printed diagnostic.
- `p̂*` = smallest grid `p` with mean metric ≥ threshold. **0.90, not 0.95** — at 0.95 the read-off
  lands on the flat top of the transition where one noisy sample moves `p̂*` a whole grid step.
- **NEVER report ad-hoc diagnostics (fit exponents, RMSE) as if they were the paper's metrics.**
  The paper scores NMI; anything else is scaffolding and must be labelled as such.
- ALWAYS re-read the prose when a figure gains panels: adding η rows falsified the "increasing η
  shifts the curve upward" claim (the fitted constants are not monotone through η=15). Flag the
  contradiction instead of letting stale text stand.
- matplotlib mathtext: `\le`/`\ge` are **not** symbols — use `\leq`/`\geq` (this also bites in
  titles and annotations, not just labels).

### Notebook ↔ paper-figure workflow

- **NEVER let a preview run overwrite `docs/overleafs/v9/figures/`.** The figure cells save to
  `FIG_DIR`; when exploring, exec the cells headlessly with `FIG_DIR` redirected to scratch.
- Editing a notebook cell **clears its stored outputs**, so the notebook shows stale figures until
  it is re-run. To show a result without republishing the paper asset: exec the cell with `FIG_DIR`
  diverted, then write the PNG back into the `.ipynb` as a base64 `display_data` output.
- ALWAYS strip `%` magics before `exec`-ing notebook cells from a script.

### Sweep economics (`analysis/`)

- Cost per `(η, m)` cell is `N_TRIALS × |P_GRID|` sub-sampled Fiedler solves, and each solve is
  `O(m³)` once `p > 0.1` (dense `scipy.linalg.eigh`; sparse `eigsh` below that). Measured: m=1536
  → 39 s, m=2816 → 3 min, m≈11k → ~3 h. ALWAYS extrapolate before launching a full grid, and
  prefer fewer sizes / fewer trials over an overnight tail — the largest point buys one dot on a
  log axis.
- Every trial is cached individually, so runs are **resumable and safely interruptible**; a smaller
  `N_TRIALS` just reads the first seeds of an existing cache, keeping old and new η comparable.
- ALWAYS check what is already cached (`cache/full_matrix`, `cache/sweep_trial`; a complete flat-CBM
  cell has 480 trial dirs = 24 four-decimal `p` keys × 20 seeds) before recomputing.
- Flat-CBM grids: `m` must be divisible by `1+η` for `n1 = m/(1+η)` to be exact, so **each η needs
  its own m grid** (÷2, ÷6, ÷11, ÷16 for η = 1, 5, 10, 15). Reusing one grid silently shifts the
  realised η (m=90 at η=15 gives 17.0).
- Long runs go in the background with progress to a log, smallest sizes first, so the figure is
  assemblable before the tail finishes.

## Conventions worth knowing

- `gemini.md` (sister AI instructions) sets stylistic preferences: keep files ≲200 lines, one function = one logic, type-hint matrix inputs as `np.ndarray`, prefer `_v2`/`temp_` naming over overwriting working code during experiments. Apply these to new `sub_sampled_fielder_vec` code.
- LaTeX is fine in docstrings/comments when describing matrix math.
- The `legacy/` subdirs in both `spectraltree/` and `tests/` are historical; don't add to them.
