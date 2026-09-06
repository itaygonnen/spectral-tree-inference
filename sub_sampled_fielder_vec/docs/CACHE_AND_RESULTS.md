# Cache and results inventory

Measured 2026-07-30. **Nothing here was deleted** — this is the record that makes a
later pruning decision safe, not the pruning itself.

Totals: `cache/` **3.2 GB**, `results/` **4.4 GB**, `data/` **594 MB** — the repo-level `data/cohorts/`
now holds **4.0 GB**: `data/cohorts/6000 taxa/` (added 2026-08-31) holds 100 alignments
of 6000 taxa x 5000 sites (2.8 GB) beside all 3000 true trees (586 MB), with `1000 taxa` alongside it, extracted from
`sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz`; the remaining 2900 alignments stay in
the archive (~90 GB extracted). All three are gitignored. `logs/` (35 files) **is** committed.

Everything is regenerable in principle, but not cheaply: rebuilding `pool_sample` alone
means re-running the η-pool builder at every size. Cost model (from `CLAUDE.md`):
m=1536 → 39 s, m=2816 → 3 min, m≈11k → ~3 h, and each solve is `O(m³)` once `p > 0.1`.

## `cache/` — seven scopes

Roots are declared in `src/cache_io.py`; `CACHE_ROOT = <package>/cache`. Every entry is
marked complete by a `.complete` sentinel, so runs are **resumable and safely
interruptible**.

| Scope | Size | Keys | Complete entries | Written by | Read by |
|---|---|---|---|---|---|
| `pool_sample` | **2.2 G** | 7 | 296 | `scripts/build_eta_pool*.py` → `src/utils/eta_pool_cache.py` | Figs 3, 5, 8, 9 |
| `sweep_trial` | 449 M | 91 | **36,857** | *nothing in `src/`* — **orphaned** | nothing |
| `experiment_data` | 441 M | 16 | 16 | `src/utils/persistent_cache.py` ← `bootstrap_sweep.py` | pipeline A |
| `distance_matrix` | 52 M | 16 | 16 | `src/runners/nj_sweep.py` | NJ sweeps |
| `full_matrix` | 42 M | 291 | 291 | *nothing in `src/`* — **orphaned** | nothing |
| `bpart_sweep` | 3.3 M | 40 | 315 | `src/utils/bpart_sweep_cache.py` | bpart notebooks |
| `bootstrap_sweep` | 2.7 M | 503 | 503 | `src/utils/sweep_cache.py` ← `build_sweeps.py` | Figs 3, 8, 9 |

**491 MB across two orphaned scopes.** `sweep_trial` (36,857 entries — 96% of all cache
entries on disk) and `full_matrix` have no reader *or* writer left in `src/`; both were
produced by the pre-`edb397b` layout and migrated by `scripts/legacy/migrate_caches.py`.
They are the obvious pruning candidates, but confirm first that no notebook reaches them
through `analysis/utils/cache.py`, which has its own key
scheme (see below).

### Three coexisting key schemes

Not unified, because unifying them would move every key and orphan the lot:

1. **`cache_io.make_key`** — sorted kwargs, e.g.
   `L10000_mu0p1000_n0500_num_classes0004_pop_size1p0000_seqJC69_treekingman`.
   Used by `experiment_data`, `full_matrix`, `distance_matrix`.
2. **`eta_pool_cache.param_key`** — hand-rolled, with its own private `_fmt_value`
   duplicating `cache_io._fmt_value`:
   `n0500_L10000_mu0p1000_kingman_pop1p0000_JC69`, then
   `eta{TT:02d}/sample_{NNNN:04d}/`.
3. **`sweep_cache.compute_sweep_key`** — flat, `<param_key>__eta{TT}__sample_{NNNN}__<sha1(config)[:12]>`,
   over `{p_values(9dp), bootstrap_reps, seed, num_gaps, min_split,
   early_stop_consecutive_100, partition_method, matrix_kind, distance_alpha,
   schema_version="v3"}`. Canonical keys at the current defaults:
   `sign=acb9d64195e6`, `sigma2=6f5687277cea`, `kmeans=234a0bd360f7`.

`laplacian` is added to the sweep key **only when non-default**, so adding the
normalized-Laplacian option did not move a single existing key.

Plus two non-disk layers: an in-memory similarity cache keyed on a hash of the
observations (`src/utils/random_entries.py`), and a resumable `screen_table.csv`
(`src/utils/benchmark_cache.py`).

### The cache-root argument is a no-op

Every shim (`persistent_cache`, `sweep_cache`, `eta_pool_cache`, `bpart_sweep_cache`)
**ignores its `cache_root` parameter** and uses `cache_io.CACHE_ROOT`. Until this pass,
`src/utils/eta_pool_sweep.py` passed `<package>/src/cache` — a path that does not exist.
Harmless only because the argument is discarded; now it imports `CACHE_ROOT` instead.
If you ever make a shim honour that argument, audit all call sites first.

## `results/` — 4.4 GB

Layout is **inconsistent — six conventions coexist**, measured:

| Convention | Count | Example |
|---|---|---|
| flat `runs/<ts>-<name>/n<n>_L<L>/` | **60** | `runs/20260609-151136-bpart_kingman_L10000/n512_L10000/` |
| nested `runs/<tree_model>/<sampling_method>/<ts>-<name>/` | 5 runs, 3 trees | `runs/kingman_mean/uniform/20260501-194101-…/` |
| **4-level** `runs/<tree_model>/<matrix_kind_alpha>/<sampling_method>/<ts>-<name>/` | 7 runs | `runs/balanced_binary/distance_a1p000/uniform/20260512-163359-…/` |
| named, **no timestamp** | 1 | `runs/real_data_benchmark/n1000/` |
| `notebooks/<NN_topic>/<notebook>/` | 5 leaves | `notebooks/01_cbm_theory/balanced_binary_threshold/` |
| stray run dir at `results/` top level | 1 | `20260522-114045-fiedler_sweep_birth_death_L10000` (partial; predates `e0bb68f`) |

Nesting is chosen by `src/cache_io.py:238-262` `run_dir()` — which has **zero callers**.
The convention is actually implemented three other places:
`src/runners/experiment_runner.py:514-545` (duplicates the logic inline; always nests,
since both kwargs are always truthy), `src/runners/experiment_runner_utils.py:52-99`
(writes to `results/<tree_model>/…` with **no `runs/` segment**, and adds the
`distance_a{alpha}` level), and the per-method sweep drivers plus three notebooks, which
hand-roll the flat form.

Nested paths **restate information the leaf name already carries** — e.g.
`runs/kingman_mean/uniform/20260501-194101-kingman_mean_n500-8000_mu_0p1_uniform/` names
the tree model and the sampling method twice each, because
`experiment_runner_utils.py:35-49 generate_run_prefix` bakes them into the prefix that is
then nested underneath them.

### Paper-backing

| Path | Why it matters |
|---|---|
| `results/runs/kingman_mean/uniform/20260501-194101-kingman_mean_n500-8000_mu_0p1_uniform/` | **Prerequisite of the identity-check figure** (App G's `fig:identity_gen` until that appendix was deleted on 2026-08-04). `supporting/identity_checks_kingman.ipynb` asserts this exists in its setup cell, so it cannot run in a fresh clone without re-running this sweep. |
| `results/notebooks/01_cbm_theory/balanced_binary_threshold/` | **Figure 6's inputs** — `trials.csv` (501 KB), `agg.csv`, `config.json`. The `01_cbm_theory` prefix is a historical cache key; renaming it orphans this. |
| `results/notebooks/03_sampling_methods/nnm_vs_ipw/` | outputs of the `nnm_vs_ipw` supporting notebook |
| `results/runs/20260609-*-bpart_*` | distance-route (bpart) evidence, June 2026 — the newest generation |

### Superseded and scratch — where the space is

**Four dirs hold 3.1 GB, 71% of `results/`:**
`20260512-193843-snj_sweep` (1.0 G), `20260514-085600-nj_sweep` (784 M),
`20260514-174813-nj_sweep` (732 M), `20260513-001823-snj_notebook` (613 M). The two
`nj_sweep` dirs are superseded by the May-19 `*_L10000` runs.

Also: named smoke tests (`*-griffing_smoke`, `*-fiedler_smoke`, `*-nj_smoke`,
`*-nj_meanimp_smoke`, `*-nj_sweep_birth_death_smoke`), and **22 same-minute duplicate
reruns** — `fiedler_sweep_birth_death_L10000` ×5, `fiedler_sweep_kingman_L10000` ×4,
`griffing_sweep_*` ×3 each, `bpart_*` ×6/×6/×5, plus `20260527-104008` and
`20260527-104013` five seconds apart (the first evidently aborted).

The `*-fiedler_sweep_*` dirs are the only output of `src/runners/fiedler_sweep.py`, whose
driver script was never committed — `scripts/plot_fiedler_overlay.py` and
`src/utils/fiedler_io.py` still read them.

## `analysis/legacy/spectral_analysis/` — 86 MB in the source tree

`target_quality_anlysis/analysis_results/` holds **12 committed Dec-2025 run
directories** (231 PDF, 197 TXT, 84 SVG, 84 NWK) — 78% of all of `analysis/`. Unlike
`results/`, this is **tracked**, so it is in every clone's history. Each run dir has its
own `config.json` and is self-contained, so nothing dangles; it is simply stale output
checked into source. Removing it from the working tree would not shrink the clone unless
history is rewritten.

## If you do prune later

1. Lowest risk: `cache/sweep_trial` + `cache/full_matrix` (**491 MB**, no reader or
   writer in `src/`) — but grep `analysis/utils/cache.py`
   first, it addresses caches by its own scheme.
2. Next: the four superseded `results/runs` dirs (**3.1 GB**) and the named smoke tests.
3. Keep unconditionally: `cache/pool_sample`, `cache/bootstrap_sweep`,
   the `kingman_mean/uniform/20260501-…` run, and `results/notebooks/`.
