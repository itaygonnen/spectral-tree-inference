# Changelog

All notable changes to the STDR framework are documented here.

---

## [2026-09-17] Middle-out early stop in the sweep, and L_sym off the sweep default

The p-grid is now swept from its middle outwards and each end stops once the curve has
settled, instead of measuring every p from `p_min` to 1.0. The strategy is the one
`middle_out_runner.py` has always had for the legacy single-tree flow; the code is new,
because that implementation has four defects (below).

### Added
- `p_sweep_inner.MiddleOut` — the traversal and both stop rules, shared by the similarity
  route (`bootstrap_p_sweep_simple`) and the distance route (`bpart_sweep_raw`). Serial,
  two cursors, one step per side. High side stops after N consecutive `NMI == 1.0`, low
  side after N consecutive `NMI <= 0.10`; `p_values[0]` is always measured so `p_min`
  anchors the floor with a real number. `traversal="ascending"` reproduces the pre-2026-09
  shape, which is what the eta-pool cache holds.
- `--early-stop N` / `--no-early-stop` on `scripts/run_sweep.py`, `RunSpec.early_stop`
  (default 2), and the matching menu question.
- `measured_{arm}`, a per-p column recording which points were computed. Additive — not in
  `operator_sweep.META_KEYS` — so no cached tree was invalidated.
- `scripts/check_early_stop.py` — replays the traversal over every cached curve and asserts
  both stop rules against data that predates the code, plus a structural check that the
  traversal partitions `range(n_p)` exactly for every `(n_p, early_stop, stop position)`.
  0 violations across the 8 cached trees.
- `scripts/migrate_sweep_cache_v4.py` — the one-off v3 -> v4 re-key.

### Changed
- **Sweeps read L(S) and B by default, not L_sym** (`experiment_run.DEFAULT_OPERATORS`).
  L_sym is a full Fiedler arm, roughly half a sweep's wall time where B is 3%
  (`GRIFFING_COST_SHARE`), and the recovery figure compares L against B. It stays in
  `ALL_OPERATORS`: the screen still records `valid_Lsym`/`eta_Lsym`, the `valid_Lsym` gate
  still works, `--operators S,Lsym,B` still sweeps it. This is the larger of the two
  savings here — early stop skips 21.7% of grid points across the cached trees but only
  0.4% of the DENSE ones (p > 0.1, where the O(m^3) solves are).
- Both inner sweeps index their outputs by p instead of appending, since visit order is no
  longer p order. Every callback fires in visit order; the log line already printed `p`.
- `sweep_cache.SCHEMA_VERSION` v3 -> v4: `early_stop_consecutive_100` renamed `early_stop`,
  `traversal` added. The 503 cached eta-pool sweeps were re-keyed and corrected in place
  rather than orphaned (~60 h of dense solves avoided); `partition_nmi_M` is byte-identical
  across the migration, and it is one of the only three fields `eta_pool_sweep` reads.

### Fixed
- The fill no longer fabricates what its trigger does not imply. An extrapolated p carries
  nmi/ari/agreement and, on the high side, the reference partition; `dot`,
  `sign_agreement`, `sigma2_avg_S` and every rank/operator-norm diagnostic are NaN. On one
  cached eta-pool entry the old fill had written `dot = 1.0` over 12 of 25 points, the
  measured value at the last real point being 0.9928.
- `load_tree_sweep` compares `early_stop` under a partial order, not for equality: a tree
  swept with a smaller (or no) stop is strictly better data and is reused, while a run that
  asked for measurement is never handed an extrapolated tail.

### Defects in `middle_out_runner.py` this does NOT inherit (that code is untouched)
1. **Index skipping.** When one side's guardrail fires it widens `workers_per_side` and
   then advances `offset` by it, so the other side's next block starts past where it
   stopped. At `n_p=20, workers=8`, firing the high guardrail in round 1 leaves indices
   2,3,4,5 never computed AND never filled. Two serial cursors cannot do that, and
   `check_early_stop.py` tests it.
2. `check_guardrails_high` fires at `>= 99.9` but `create_guardrail_high_result` writes
   `100.0` — a fabrication at the trigger, not only in the tail.
3. That fill asserts `S = M`: `operator_norm_error = 0`, every `_S` metric equal to its
   `_M` counterpart. Untrue at p = 0.3 with a perfectly recovered partition.
4. `create_below_threshold_result` writes `partition_agreement = 50.0` as a "random
   baseline". Flip-max agreement's random baseline at the 20:1 splits `MAX_ETA` allows is
   ~95%, and this pipeline scores NMI, whose measured floor is 0.000-0.048.

---

## [2026-09-17] Removed the leveraged (IALM) and LDS samplers

Uniform sampling is now the only sub-sampling method. Nothing in the live pipeline
(`experiment_run.py` → `operator_screen` / `griffing_sweep` / `p_sweep_inner`) had
selected `leveraged` or `lds` since the screen→gate→p-sweep flow landed, and every
result the manuscript rests on is uniform.

### Removed
- `src/core/sampling/leveraged/` (14 files, ~1,400 lines): IALM solver, singular-value
  thresholding, leverage-score estimation, non-uniform sampling, the LDS/IPW sampler.
- `src/utils/sampling_logger.py` — existed only to dump leverage scores per `p`.
- `SamplingConfig` fields `theta`, `target_rank`, `tau_floor_multiplier`, `ialm_max_iter`,
  `ialm_tol`, `ialm_bypass_threshold`, `force_leveraged`, `allow_uniform_fallback`,
  `log_sampling_diagnostics`, `prob_formula`, and their `custom_config` / `SWEEP_CONFIG`
  pass-throughs.
- 15 always-NaN diagnostic columns from `results.json` (`phase1_s1..s3`, `leverage_max`,
  `leverage_std`, `leverage_sum`, `leverage_symmetry_error`, `ialm_iterations`,
  `phase1_sufficiency`, `spectral_gap`, `phase1_actual`, `phase2_actual`,
  `fallback_to_uniform`, `tau_floor`, `phase1_budget_fraction`). 133 columns → 88.
- ~225 lines from `bootstrap_sweep.py`: the sampler-kwargs branching, the per-replicate
  leverage-metric collection, the diagnostics-save block, the leveraged/LDS run summary.

### Changed
- `SamplingConfig.method` is `Literal["uniform"]`. It is kept because it names the
  results directory (`results/<tree_model>/<method>/`); no cache key encodes it, so
  every existing cache stays valid.
- `get_sampler` raises `ValueError` on anything but `"uniform"`; `SamplingConfig`
  rejects `"leveraged"`/`"lds"` at validation time.
- `results.json` config columns are now `sampling_method` + `sampling_matrix_kind`
  (was method/theta/target_rank/tau_floor_multiplier/prob_formula).

### Removed — the analysis and docs for it
- `analysis/legacy/leveraged_sampling_analysis/` (19 files, 2,679 lines): the IO + metrics +
  visualization stack for leveraged/LDS runs, including the 4th copy of
  `plot_phase_transition_scaling`.
- `analysis/legacy/scripts/analyze_lds_experiment.py` (374 lines) — the directory is now gone.
- `docs/deprecated/` in full: `LEVERAGED_SAMPLING.md`, `LDS_SAMPLING.md`, `LDS_MIGRATION.md`
  (1,500+ lines documenting a deleted API). `docs/INDEX.md`, `docs/ARCHITECTURE.md`,
  `docs/CONFIGURATION.md`, `analysis/legacy/README.md` and `CLAUDE.md` updated to match.

`analysis/legacy/comparison/notebook_utils.py` keeps two `sys.path.insert` calls naming
`leveraged_sampling_analysis/io`. They were already broken before this pass (they omit the
`legacy/` path segment); recorded in `analysis/legacy/README.md`, not fixed — legacy is
quarantined.

### Not touched
- `analysis/supporting/nnm_vs_ipw.ipynb` — its NNM arm uses Soft-Impute in
  `analysis/utils/spectral.py`, never the deleted IALM solver.

### Verified
- All touched modules import; `python -m src.cache_io` passes; `run_sweep.py --dry-run`
  and `sync_paper_figures.py --check` unchanged (13 figures, 4 floats, exit 0).
- End-to-end `sweep_for_params` at n=128, L=400, p∈{0.5, 1.0}: `partition_agreement_M`
  = 100% at both, incremental save writes `results.json` with 88 columns.

---

## [2026-07-30] Repository restructure around the v9 paper

Two cleanup passes. **No mathematics, figure, or plotted curve changed** — all 23 paper
figure md5s are unchanged and the PDF still builds to 29 pages / 2,586,870 bytes.

### Fixed — the paper was not in version control

- `.gitignore` had three over-broad patterns. A bare `figures` excluded
  `docs/overleafs/v9/figures/`, and the `.tex` files were never added, so the manuscript
  lived entirely outside git. A bare `test*.py` matched `tests/test_*.py`, which is why that
  directory has never held a single test.
- Committed modules that tracked notebooks already imported but that were never added:
  `analysis/utils/sweep_plots_two_panel.py` (Figs 2, 3), `src/utils/screening.py`,
  `src/utils/eta_pool_griffing.py`, and `paper/fig04_hbm_spectral_gap.ipynb` — the sole
  producer of Figure 4. A fresh clone previously raised `ImportError` on five notebooks.
- v9 is now self-contained: Figure 5 resolved through `\graphicspath{{../}}` to assets
  outside the tree, so uploading `v9/` alone to Overleaf failed.
- `src/` no longer imports from `analysis/`. Two modules (`partition_validity`, `tree_plots`)
  were pulled in by bare name via `sys.path` hacks, so `import src` only ever worked by
  accident; both now live under `src/`.
- Untracked 10.9 MB of third-party PDFs, 35 machine-specific run logs, and `last_run.json`
  (mutable state) — all committed by mistake or by habit.

### Added

- `analysis/PAPER_MAP.md` — figure ↔ notebook ↔ claim ↔ cache ↔ rebuild command, *generated*
  from `analysis/paper_figures.py` and validated by
  `scripts/sync_paper_figures.py --check`, which cross-checks the manifest against every
  `\includegraphics` in the paper. The previous mapping was stale and misattributed Figure 7.
- `docs/RUNBOOK.md`, `docs/CACHE_AND_RESULTS.md`, `scripts/README.md`,
  `analysis/legacy/README.md`, `tests/README.md`.
- `docs/overleafs/v9/open-items/19-cleanup.md` — 10 findings referred to the author.

### Changed

- `analysis/` flattened to `paper/` (6 notebooks named `figNN_*`), `supporting/` (12),
  `utils/` and `notebooks_cache/`. The `theoretical_interpretation/` layer and all per-topic
  subdirectories are gone — 18 notebooks had been spread over 16 directories, and `.npz`
  data sat beside 2 notebooks but not the other 16. Imports are now `analysis.utils`.
- `scripts/` now holds only things you run; the six imported `plot_*` libraries moved to
  `src/plots/` and `merge_results.py` to `src/utils/`.
- Six superseded analysis packages archived under `analysis/legacy/`.
- Deleted as dead: `configs/` (zero code references), `examples/` (zero importers),
  `analysis/utils/distance_similarity.py` (843 lines, documented as shared but never
  imported), `docs/ANALYSIS_GUIDES.md` (documented notebooks that no longer exist).
- `docs/METRICS.md` now documents **NMI** as the metric the paper reports; it previously
  mentioned neither NMI nor ARI.

---

## [2026-02-20] Critical Bug Fix: LDS Debiasing Probability

### Fixed

**CRITICAL: Incorrect debiasing probability in LDS** (`src/core/sampling/leveraged/lds_sampler.py:315-325, 339`)

**The Bug:**
- LDS uses two phases: Phase 1 (uniform, probability $p_0$) and Phase 2 (leveraged, probability $p_{ij}$)
- Algorithm requires debiasing with **effective inclusion probability**: $\pi_{ij} = p_0 + (1-p_0) \cdot p_{ij}$
- Code was debiasing **all entries** (both phases) with only Phase 2 probability $p_{ij}$
- **Impact**: Biased estimator $\mathbb{E}[\hat{M}] \neq M$ → violates unbiasedness guarantee
- **Symptom**: Worse performance at low p-values where Phase 1 dominates the sample

**The Fix:**
```python
# Before (WRONG):
X_hat_sparse = compute_debiased_estimator(matrix, Omega, p_matrix)  # Uses p_{ij} only

# After (CORRECT):
p_0 = phase1_actual / n_upper
pi_matrix = p_0 + (1.0 - p_0) * p_matrix  # Effective inclusion probability
X_hat_sparse = compute_debiased_estimator(matrix, Omega, pi_matrix)
```

**Why It Matters:**
- **Unbiasedness**: The theoretical guarantee $\mathbb{E}[\hat{M}] = M$ is critical for Davis-Kahan eigenvector stability
- **Inclusion-Exclusion**: Entry (i,j) is sampled if: (1) sampled in Phase 1 OR (2) missed in Phase 1 AND sampled in Phase 2
- **Low p impact**: At low p-values, Phase 1 budget dominates (30% of total) → most entries debiased with wrong probability
- **Bias direction**: Phase 1 entries debiased with $p_{ij}$ (high-leverage) instead of $p_0$ (uniform) → under-amplification → signal loss

**Validation:**
- Expected to improve partition agreement at low p-values (p < 0.01)
- Should see more stable performance across p-value range
- Theoretical properties now correctly implemented

---

## [2026-02-18] Critical Bug Fix: Leverage Score Computation

### Fixed

**CRITICAL: sklearn TruncatedSVD normalization bug** (`src/core/sampling/leveraged/compute_leverage_scores.py:78-87`)

**The Bug:**
- `sklearn.TruncatedSVD.fit_transform()` returns **U*Σ** (coordinates in reduced space), NOT orthonormal **U**
- Leverage scores were computed from scaled vectors: `μ_i = (n/r) * ||U*Σ[i,:]||²`
- This squared the singular values into the leverage scores: `||U*Σ||² = ||U||² * σ²`
- **Impact**: Leverage scores were astronomically wrong (billions instead of ~1.0)
- **Symptom**: Negative correlation with ground truth, RMSE in billions (10⁹)

**The Fix:**
```python
# Before (WRONG):
U = svd_model.fit_transform(X_sparse)  # Returns U*Σ, not U!

# After (CORRECT):
U_sigma = svd_model.fit_transform(X_sparse)
s = svd_model.singular_values_
s_safe = s.copy()
s_safe[s_safe < 1e-12] = 1.0  # Numerical stability
U = U_sigma / s_safe[None, :]  # Normalize to get orthonormal U
```

**Why It Matters:**
- Leverage scores should average ~1.0 (since they sum to n)
- The bug caused a scaling chain: IPW amplification (X/p) → huge singular values → squared in leverage calculation
- Phase 2 sampling probabilities were completely wrong, making LDS effectively random sampling

**Validation:**
- Ground truth comparison now shows **positive correlation** (r > 0.9)
- RMSE reduced from 10⁹ to ~1.0
- Leverage score distributions now match theoretical expectations

**Discovered through:**
- Comparative analysis in `analysis/notebooks/lds_experiment_comparison.ipynb`
- Ground truth validation against full SVD leverage scores

---

## [2026-02-14] HLDT → LDS Rename

### Changed
- **Renamed HLDT → LDS (Leveraged Debiased Sampler)**
  - More descriptive name capturing the algorithm's two phases (leveraged sampling → debiased estimator)
  - Updated all code, documentation, and configuration files
  - **Breaking Changes**:
    - Configuration: `sampling_method="hldt"` → `sampling_method="lds"`
    - Python imports: `from .leveraged.hldt_sampler import HLDTSampler` → `from .leveraged.lds_sampler import LDSSampler`
    - Class name: `HLDTSampler` → `LDSSampler`
    - Parameter: `force_hldt` → `force_lds`
    - File renames:
      - `hldt_sampler.py` → `lds_sampler.py`
      - `HLDT_SAMPLING.md` → `LDS_SAMPLING.md`
      - `HLDT_MIGRATION.md` → `LDS_MIGRATION.md`
  - **Why**: "Leveraged Debiased Sampler" describes **what** the algorithm does, while "HLDT" only indicates **who** invented it
  - **Historical Note**: Below entries use "HLDT" to reflect terminology at the time of implementation. The underlying algorithm (from Huang, Liu, Du, and Tao, 2018) remains unchanged.

---

## [2026-02-14] Bug Fixes

### Fixed
1. **HLDT parameter passing bug** (`src/runners/bootstrap_sweep.py:247-272`)
   - Fixed `TypeError: HLDTSampler.__init__() got an unexpected keyword argument 'ialm_max_iter'`
   - Issue: IALM-specific parameters were being passed to HLDTSampler
   - Solution: Conditionally build method_kwargs based on sampling method
   - Now correctly passes:
     - Common params to all samplers: theta, target_rank, allow_uniform_fallback
     - IALM-only params to leveraged: ialm_max_iter, ialm_tol, ialm_bypass_threshold, force_leveraged
     - HLDT-only params to hldt: tau_floor_multiplier, force_hldt

2. **Sparse matrix support for metrics** (`src/core/utils.py:36-63`, `src/core/metric_computer.py:29-91`)
   - Fixed `ValueError: Sparse arrays/matrices are not supported by this function`
   - Issue: HLDT returns sparse matrices but metric computation used dense-only scipy.linalg functions
   - Solution:
     - Updated `compute_laplacian()` to handle both sparse and dense matrices
     - Added `_ensure_dense()` method to MetricComputer to convert sparse to dense for metrics
     - FiedlerVectorComputer already handled sparse correctly
   - Trade-off: Metrics conversion to dense acceptable (computed once per p-value, not in tight loop)

3. **HLDT diagnostic logging support** (`src/runners/bootstrap_sweep.py:718-721`)
   - Fixed diagnostic logging not saving HLDT leverage scores
   - Issue: Diagnostic logging looked for 'leverage_scores' but HLDT uses 'leverage_scores_regularized'
   - Solution: Check for both key names with fallback: `metrics.get('leverage_scores') or metrics.get('leverage_scores_regularized')`

4. **Misleading sampling summary logs** (`src/runners/bootstrap_sweep.py:963-979`)
   - Fixed log always showing "Method: leveraged" regardless of actual sampling method
   - Issue: Hardcoded "leveraged" string in summary, IALM statistics shown even for HLDT
   - Solution:
     - Use `cfg.sampling.method` to display correct method name
     - Conditionally show IALM statistics only for `method="leveraged"`
     - Show HLDT-specific message for `method="hldt"`
   - Impact: Logs now correctly show "Method: hldt" when HLDT is used

---

## [2026-02-13] HLDT Sampling Implementation

**Major Feature**: Added HLDT (Huang-Liu-Du-Tao) sampling method for fast spectral preservation

### Added
- **HLDTSampler** class (`src/core/sampling/leveraged/hldt_sampler.py`)
  - Single-shot debiased estimator (10-100x faster than IALM)
  - Three-phase pipeline: uniform → leveraged → debiasing
  - Sparse matrix output for efficiency
  - Comprehensive diagnostics and fallback support

- **Core modules**:
  - `compute_debiased_estimator.py` - Debiased estimator X̂ = X/p
  - Regularization floor for leverage scores (τ_floor = mean(μ))

- **Configuration**:
  - `sampling_method="hldt"` option
  - `sampling_tau_floor_multiplier` parameter
  - `sampling_allow_uniform_fallback` support for HLDT

- **Documentation**:
  - `docs/HLDT_SAMPLING.md` - Complete user guide (469 lines)
  - `docs/HLDT_MIGRATION.md` - Technical migration details

- **Testing**:
  - `scripts/test_hldt_quick.py` - Validation test suite
  - Phase transition validation (50% → 97% agreement)

### Performance
- **Speed**: 57.5x faster than IALM at p=0.05
- **Accuracy**: 95%+ agreement at moderate p-values (p>0.2)
- **Scalability**: Suitable for large trees (n>500 taxa)

### Modified
- `src/core/sampling/leveraged/compute_leverage_scores.py` - Added regularization support
- `src/core/sampling/__init__.py` - Added HLDT to sampler factory
- `src/config/base_config.py` - Added HLDT parameters
- `src/runners/bootstrap_sweep.py` - Added HLDT diagnostic logging
- `scripts/interactive_run.py` - Added HLDT to both UI paths

---

## [2026-02-13] Leveraged Sampling Configuration Flow

**Focus**: Enhanced configuration control for leveraged sampling research mode

### Added
- **`allow_uniform_fallback`** parameter
  - Safe Mode (`True`, default): Falls back to uniform at low p
  - Research Mode (`False`): Proceeds with leveraged sampling using 90/10 split

- **`InsufficientBudgetError`** custom exception
  - Prevents accidental catching by generic handlers
  - More precise error handling

### Modified
- **Configuration threading** through 8 layers:
  1. `scripts/interactive_run.py` - UI collection
  2. `src/runners/experiment_runner_utils.py` - Extraction
  3. `src/config/presets.py` - `custom_config()` signature
  4. `src/config/base_config.py` - `SamplingConfig` dataclass
  5. `src/runners/bootstrap_sweep.py` - Sampler instantiation
  6. `src/core/sampling/leveraged/sampler.py` - Final usage

### Fixed
- **Cache selection bug**: Menu sorted by n_taxa but selection used unsorted list
  - Added sorting in selection handler (`scripts/interactive_run.py:342`)

### Documentation
- Updated `docs/LEVERAGED_SAMPLING.md` with Research Mode section
- Updated `docs/ARCHITECTURE.md` with SamplingConfig details
- Updated `docs/CONFIGURATION.md` with mode examples

---

## [Previous Update] Interactive Launcher & Diagnostic Logging

**Focus**: User experience improvements and diagnostic validation system

### Added

#### Interactive Launcher
- **`scripts/interactive_run.py`** - Menu-driven experiment configuration
  - Gemini-style colored UI
  - Last run re-execution (`[r]` option)
  - Cached matrix selection for instant loading
  - Interactive parameter entry with defaults

#### Persistent Caching System
- **`src/utils/persistent_cache.py`** - Disk-based matrix storage
  - Location: `src/cache/{experiment_params}/`
  - Cached data: Similarity matrix, observations, tree, Fiedler reference
  - Instant loading for re-runs (no recomputation)

#### Diagnostic Logging System
- **`src/utils/sampling_logger.py`** - Diagnostic data saver
  - Saves leverage scores and Phase 2 probabilities
  - Sparse `.npz` format per p-value
  - Location: `{run_dir}/sampling_data/p_{p:.4f}.npz`

- **`analysis/notebooks/leverage_sampling_explorer.ipynb`**
  - Validation tool for leverage score quality
  - Three-panel visualization (heatmap, scatter, probs)
  - Correlation analysis validates Phase 1 quality

#### Configuration
- **`log_sampling_diagnostics`** flag
  - Enable diagnostic data logging
  - Works with both interactive and script-based workflows

### Modified
- `src/core/sampling/leveraged/sampler.py` - Extended metrics tracking
- `src/runners/bootstrap_sweep.py` - Diagnostic logging integration
- `src/utils/interactive_ui.py` - UI components (logo, colors, inputs)

### Documentation
- Updated `README.md` with interactive launcher quick start
- Updated `docs/ARCHITECTURE.md` with caching and diagnostic system
- Updated `docs/LEVERAGED_SAMPLING.md` with diagnostic logging section
- Updated `docs/ANALYSIS_GUIDES.md` with explorer notebook guide

### Research Findings
From n=1024, p=0.1438 experiment:
- **Correlation**: -0.0728 (nearly zero!)
- **Sampling rate**: 3.54%
- **Insight**: Phase 1 uniform sampling at low p doesn't reliably estimate leverage scores
- **Validates**: Theoretical minimum budget requirement (4·n·r·log(n))

---

## Earlier History

For earlier changes, see commit history:
```bash
git log --oneline --decorate --graph
```

Key milestones:
- Initial STDR framework implementation
- Partition-based metrics (gap-based thresholding)
- Bootstrap averaging with Welford's algorithm
- Leveraged sampling with IALM
- Parallel processing support
- Comprehensive metrics and diagnostics

---

## Documentation Improvements (Ongoing)

### 2026-02-13
- **Created `docs/INDEX.md`** - Central navigation hub for all 19 documentation files
- **Consolidated changelogs** - This file replaces scattered update summaries
- **Deprecated `docs/LEVERAGED_SAMPLING_USAGE.md`** - Content moved to main guide

### Quality Improvements
- Consistent cross-referencing between docs
- Clear categorization (Getting Started / User Guides / Advanced / Meta)
- Quick reference by task
- AI agent navigation guide

---

## Backward Compatibility

All changes maintain backward compatibility:
- ✅ Default behaviors unchanged
- ✅ Existing configurations continue to work
- ✅ New parameters are optional with sensible defaults
- ✅ Old script-based workflow still supported

---

**For detailed technical information, see individual documentation files in `docs/`**

**For complete commit history**: `git log --all --decorate --oneline --graph`
