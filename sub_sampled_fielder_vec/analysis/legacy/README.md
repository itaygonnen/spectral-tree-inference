# analysis/legacy — superseded, unmaintained, kept for reference

Everything here predates the v9 manuscript and is **superseded by
`analysis/`**, which holds the paper's notebooks. Nothing in
this directory is imported by `src/`, `scripts/`, or any paper notebook — that was
verified before the move, and the one live dependency was extracted first (see below).

Do not add to these. If you need something here, port it into
`analysis/utils/` or `src/utils/` rather than reviving the import.

Last substantive activity: **Dec 2025 – Apr 2025**, versus July 2026 for the paper
notebooks.

## What was extracted before archiving

`comparison/phase_transition_utils.py` was **not** dead. Four live importers pulled
from it, two of them paper-critical (`utils/sweep_plots.py` → Fig 3,
`utils/recovery.py` → every notebook). Its three exported functions —
`find_discrete_threshold`, `fit_power_law`, `evaluate_power_law` — now live in
`src/utils/threshold_utils.py`, differential-tested against the originals over ~600
random cases including the NaN and insufficient-data paths. The copies here are the
originals, retained only as history; the canonical versions are in `src/`.

## Contents

| Path | What it does | Condition |
|---|---|---|
| `comparison/` | uniform vs leveraged head-to-head runner, sigmoid/power-law phase-transition fits, notebook glue | `notebook_utils.py` loads plot modules **by file path**, which is brittle. Its 3 `.md` files are redirect stubs to `docs/ANALYSIS_GUIDES.md`, **which has been deleted** -- it documented notebooks that no longer exist. Same for the two `spectral_analysis/*/README.md` stubs. |
| `generic_analysis/` | 15-module numerical-linear-algebra diagnostics library (Davis–Kahan, eigen spectrum, IPR, rank/coherence, operator norm) | **Never imported in this checkout** — no importer anywhere and no `__pycache__` was ever generated. 1,630 lines. |
| `leveraged_sampling_analysis/` | IO + metrics + visualization stack for the LDS / leveraged-sampling experiments | Was in use (has `__pycache__` from May 2026). `metrics/diagnostics.py` and `metrics/leverage.py` both define `compute_leverage_concentration`. Its docstring example path `results/kingman_mean/lds/…` no longer exists; the live layout nests under `results/runs/`. |
| `notebooks/` | `matrix_utils.py` only, 728 lines, "for target matrix analysis notebooks" | **Contains no notebooks.** The helper library has no consumer. |
| `scripts/` | `analyze_lds_experiment.py`, cross-`n` LDS analysis | Overlaps `leveraged_sampling_analysis/` and `plot_phase_transition.py`. |
| `spectral_analysis/` | two pipelines: `sweep_params_analysis/` (p_crit detection, BBP crossing, scaling laws, 10 plot modules) and `target_quality_anlysis/` (pre-flight diagnostics on full similarity matrices, plus a stability K-trial variant) | **86 MB — 78% of all of `analysis/`.** Holds 12 committed Dec-2025 run directories under `target_quality_anlysis/analysis_results/` (231 PDF, 197 TXT, 84 SVG, 84 NWK). Directory name has a typo (`anlysis`), kept because renaming it breaks those run dirs. |
| `plot_phase_transition.py` | single-run phase-transition plot script | Was at `analysis/` root. |

## Known duplication (recorded, not fixed)

- **Four** implementations of phase-transition plotting/fitting coexist:
  `generic_analysis/phase_transition_plotting.py`,
  `leveraged_sampling_analysis/visualization/scaling.py` (both define a
  `plot_phase_transition_scaling`), `spectral_analysis/sweep_params_analysis/phase_transition/`,
  and `plot_phase_transition.py`.
- `spectral_analysis/target_quality_anlysis/output/partition_validity.py::check_partition_valid_in_tree`
  duplicates `src/utils/partition_validity.py`, which is the version the paper
  notebooks import.
- `target_quality_anlysis/cli/target_analysis_main.py` and `stability_analysis_main.py`
  are 294 lines each and near-identical; likewise `output/tables.py` vs
  `output/stability_tables.py`.

The 86 MB of committed run output is inventoried in `docs/CACHE_AND_RESULTS.md` and was
deliberately **not** deleted in this pass.
