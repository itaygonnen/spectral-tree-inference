# Documentation index

Start at the package [README](../README.md). This page says which document answers which
question — and, honestly, which ones are stale.

## Current — written or verified 2026-07-30

| Document | Answers |
|---|---|
| [../README.md](../README.md) | What this is, what the three pipelines are, how to verify the repo |
| [RUNBOOK.md](RUNBOOK.md) | How to rebuild Figure N from nothing, and what it costs |
| [BENCHMARK_GUIDE.md](BENCHMARK_GUIDE.md) | `scripts/run_benchmark.py` — every prompt, eta pooling, outputs, resume rules |
| [CACHE_AND_RESULTS.md](CACHE_AND_RESULTS.md) | What is in `cache/` and `results/`, what is orphaned, what is safe to prune |
| [../analysis/PAPER_MAP.md](../analysis/PAPER_MAP.md) | Which notebook produces which figure, and which claim it supports |
| [../analysis/README.md](../analysis/README.md) | Notebook layout (section → source) and its traps |
| [../scripts/README.md](../scripts/README.md) | What each script does; which `plot_*` files are libraries, not scripts |
| [../analysis/legacy/README.md](../analysis/legacy/README.md) | What the superseded analysis packages were |
| [../tests/README.md](../tests/README.md) | Why there is no test suite, and what verification does exist |

## Reference — accurate on the engine, predates the paper restructure

These describe `src/` and pipeline A correctly, but their paths and examples predate the
July-2026 reorganization. Trust the mechanism, verify the paths.

| Document | Caveat |
|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | System architecture. Covers pipeline A only. |
| [CONFIGURATION.md](CONFIGURATION.md) | `StructuredConfig` / Pydantic reference. Still current. |
| [METRICS.md](METRICS.md) | Metric definitions. Predates NMI becoming the paper's score. |
| [INTERACTIVE_GUIDE.md](INTERACTIVE_GUIDE.md) | `scripts/interactive_run.py` walkthrough. **Rewritten 2026-09-16**: both menu branches now run one experiment (screen → gate → sweep → one run directory). The single-tree flow it used to describe is `interactive_run.single_tree_menu()` / `scripts/run_experiment.py`. |
| [../CHANGELOG.md](../CHANGELOG.md) | Notably the 2026-02-20 LDS debiasing fix (`π_ij`, not `p_ij`). |

## Sampling

`uniform` is the only method: sample `O(p·n²)` upper-triangular entries uniformly and
debias by `1/p`. Code in `src/core/sampling/uniform/`.

The `leveraged` (IALM) and `lds` samplers, their `analysis/legacy/` diagnostic stack,
and the whole `docs/deprecated/` directory (`LEVERAGED_SAMPLING.md`, `LDS_SAMPLING.md`,
`LDS_MIGRATION.md`) were **removed on 2026-09-17**. History is in git.

## Other artifacts here

- `thesis_seminar.html` — 21-slide self-contained deck (inline CSS/JS/canvas, MathJax
  from CDN). Documented in the project `CLAUDE.md`.
- `papers/` — third-party reference PDFs. Only two are tracked: `README.md` and
  `Spectral_top-down_recovery_of_latent_tree_models.pdf` (the STDR paper this work
  builds on). The other five are gitignored — they were briefly tracked by mistake
  when a `git add -A` swept 10.9 MB of them into commit `b53e2f5`.
- `overleafs/v7/`, `v8/`, `distance approach/` — superseded manuscript versions,
  gitignored. **`v9/` and `v10/` are tracked; `v10/` is the live manuscript.**

## Previously listed here and now gone

The earlier index linked 11 paths that no longer resolve: the three `LDS_*` /
`LEVERAGED_SAMPLING` docs (deleted 2026-09-17 with the samplers they documented), five
READMEs inside the analysis packages (now under `analysis/legacy/`), and
`DOCUMENTATION_UPDATE_SUMMARY.md` / `RECENT_CHANGES.md` (deleted). Recorded so a reader
who remembers them does not go hunting.
