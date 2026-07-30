# Documentation index

Start at the package [README](../README.md). This page says which document answers which
question — and, honestly, which ones are stale.

## Current — written or verified 2026-07-30

| Document | Answers |
|---|---|
| [../README.md](../README.md) | What this is, what the three pipelines are, how to verify the repo |
| [RUNBOOK.md](RUNBOOK.md) | How to rebuild Figure N from nothing, and what it costs |
| [CACHE_AND_RESULTS.md](CACHE_AND_RESULTS.md) | What is in `cache/` and `results/`, what is orphaned, what is safe to prune |
| [../analysis/PAPER_MAP.md](../analysis/PAPER_MAP.md) | Which notebook produces which figure, and which claim it supports |
| [../analysis/README.md](../analysis/README.md) | Notebook layout (section → source) and its traps |
| [../scripts/README.md](../scripts/README.md) | What each script does; which `plot_*` files are libraries, not scripts |
| [../analysis/legacy/README.md](../analysis/legacy/README.md) | What the superseded analysis packages were |
| [../tests/README.md](../tests/README.md) | Why there is no test suite, and what verification does exist |
| [overleafs/v9/OPEN_ITEMS.md](overleafs/v9/OPEN_ITEMS.md) | Open questions and known defects in the paper (116 entries, generated) |

## Reference — accurate on the engine, predates the paper restructure

These describe `src/` and pipeline A correctly, but their paths and examples predate the
July-2026 reorganization. Trust the mechanism, verify the paths.

| Document | Caveat |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | System architecture. Covers pipeline A only. |
| [CONFIGURATION.md](CONFIGURATION.md) | `StructuredConfig` / Pydantic reference. Still current. |
| [METRICS.md](METRICS.md) | Metric definitions. Predates NMI becoming the paper's score. |
| [INTERACTIVE_GUIDE.md](INTERACTIVE_GUIDE.md) | `scripts/interactive_run.py` walkthrough. Still current. |
| [../CHANGELOG.md](../CHANGELOG.md) | Notably the 2026-02-20 LDS debiasing fix (`π_ij`, not `p_ij`). |

## Sampling methods

| Method | Speed | Accuracy | Use when | Docs |
|---|---|---|---|---|
| `uniform` | fastest | baseline | default, quick tests | `src/core/sampling/uniform/` |
| `lds` | fast | high | large trees, production | [deprecated/LDS_SAMPLING.md](deprecated/LDS_SAMPLING.md) |
| `leveraged` (IALM) | slow | highest | small trees, max accuracy | [deprecated/LEVERAGED_SAMPLING.md](deprecated/LEVERAGED_SAMPLING.md) |

Those two documents were moved to `deprecated/` — the prose was superseded, but the
samplers themselves are live in `src/core/sampling/`. `deprecated/` also holds
`LDS_MIGRATION.md`.

## Other artifacts here

- `thesis_seminar.html` — 21-slide self-contained deck (inline CSS/JS/canvas, MathJax
  from CDN). Documented in the project `CLAUDE.md`.
- `papers/` — third-party reference PDFs. Only two are tracked: `README.md` and
  `Spectral_top-down_recovery_of_latent_tree_models.pdf` (the STDR paper this work
  builds on). The other five are gitignored — they were briefly tracked by mistake
  when a `git add -A` swept 10.9 MB of them into commit `b53e2f5`.
- `overleafs/v7/`, `v8/`, `distance approach/` — superseded manuscript versions,
  gitignored. **Only `v9/` is tracked.**

## Previously listed here and now gone

The earlier index linked 11 paths that no longer resolve: the three `LDS_*` /
`LEVERAGED_SAMPLING` docs (moved to `deprecated/`), five READMEs inside the analysis
packages (now under `analysis/legacy/`), and `DOCUMENTATION_UPDATE_SUMMARY.md` /
`RECENT_CHANGES.md` (deleted). Recorded so a reader who remembers them does not go
hunting.
