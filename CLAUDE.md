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

There is no formal test runner for `sub_sampled_fielder_vec/tests/` — the directory exists but is largely empty. Validate changes by running a small experiment (e.g. `taxa_values=[256], bootstrap_reps=2, p_values=[0.1, 0.5, 1.0]`).

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
`sub_sampled_fielder_vec/analysis/` holds Jupyter notebooks that consume `results/` JSON. The `theoretical_interpretation/` subdir is the current focus (see git status: figure 1/2/3 notebooks and `utils/{block_model,linalg_features,tree_features}.py`). These are research artifacts — treat them as scratch unless told otherwise.

## Conventions worth knowing

- `gemini.md` (sister AI instructions) sets stylistic preferences: keep files ≲200 lines, one function = one logic, type-hint matrix inputs as `np.ndarray`, prefer `_v2`/`temp_` naming over overwriting working code during experiments. Apply these to new `sub_sampled_fielder_vec` code.
- LaTeX is fine in docstrings/comments when describing matrix math.
- The `legacy/` subdirs in both `spectraltree/` and `tests/` are historical; don't add to them.
