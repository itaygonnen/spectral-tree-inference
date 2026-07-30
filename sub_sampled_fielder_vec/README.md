# sub_sampled_fielder_vec

Research code behind the manuscript in **`docs/overleafs/v9/`**. The paper is the
deliverable; everything here exists to produce, defend, or reproduce its claims.

**Core question.** STDR partitions a phylogeny by the sign of the Fiedler vector of
`L(S)`, where `S` is a similarity matrix built from sequence data. If we only observe a
random `p`-fraction of `S`'s entries, how large must `p` be for the recovered
bipartition to still be correct — and how does that threshold scale with the number of
taxa `m` and the top-split imbalance `η = max(n₁,n₂)/min(n₁,n₂)`?

The paper's answer (Theorem 1): `p = Θ(log m / m)` suffices at fixed imbalance, with `η`
entering as `(1+η)³` through coherence inflation and spectral-gap erosion.

## Start here

| I want to… | Read |
|---|---|
| know where a paper figure comes from | `analysis/PAPER_MAP.md` |
| rebuild a figure from scratch | `docs/RUNBOOK.md` |
| understand the notebook layout | `analysis/README.md` |
| know what a script does | `scripts/README.md` |
| know what's in cache/ and results/ | `docs/CACHE_AND_RESULTS.md` |
| see open questions on the paper | `docs/overleafs/v9/OPEN_ITEMS.md` |

Verify the repo is consistent:

```bash
python scripts/sync_paper_figures.py --check     # figure provenance; exits non-zero on problems
python -m src.cache_io                           # cache key/sentinel smoke checks
```

## Repository shape

```
sub_sampled_fielder_vec/
├── docs/overleafs/v9/        THE PAPER — tex, figures, open-items register (tracked)
├── analysis/                the paper's notebooks, by section then data source
│   └── legacy/              superseded, unmaintained (see its README)
├── src/                     the engine: config, sampling, models, runners, caching
├── scripts/                 entry points + the plot libraries notebooks import
├── cache/  results/  data/  gitignored, large, regenerable-but-expensive
└── logs/                    run logs (gitignored)
```

`src/` must never import from `analysis/`. The reverse is fine and normal.

## Three pipelines

There is no single "the experiment runner" — three distinct pipelines coexist, writing
different cache scopes. Knowing which one you are in is most of the battle.

| | Pipeline | Entry point | Config lives in | Cache scopes | Serves |
|---|---|---|---|---|---|
| **A** | legacy bootstrap sweep | `scripts/run_experiment.py`, `scripts/interactive_run.py` | a `SWEEP_CONFIG` dict at the top of the script, or the interactive menu — **no CLI args** | `experiment_data`, `sweep_trial` | historical results, and **Figure 7's prerequisite run** |
| **B** | η-pool + sweeps | `scripts/build_eta_pool.py` → `scripts/build_sweeps.py` | module constants in `src/utils/eta_pool_sweep.py` (`METHOD_SPECS`, `ETA_TARGETS`, `P_VALUES`, `BOOTSTRAP_REPS`) | `pool_sample`, `bootstrap_sweep` | **Figures 3, 5, 8, 9** |
| **C** | per-method sweeps | `scripts/run_{nj,snj,griffing}_sweep.py`, `run_benchmark.py` | `src/config/presets.custom_config(...)` inline | `distance_matrix`, `bpart_sweep` | supporting notebooks (distance route) |

Pipeline B is the one the paper's main figures depend on. Its call chain:

```
scripts/build_sweeps.py:main
  -> src/utils/eta_pool_sweep.py:run_and_collect_sweeps
       -> src/utils/eta_pool_cache.py:list_completed_samples / load_pool_entry
       -> src/utils/sweep_cache.py:compute_or_load_sweep
            -> src/runners/p_sweep_inner.py:bootstrap_p_sweep_simple   <- the actual compute
```

Every trial is cached individually behind a `.complete` sentinel, so **runs are
resumable and safely interruptible**, and lowering `N_TRIALS` just reads the first seeds
of an existing cache.

## Method: one p-sweep

1. Simulate a tree and evolve sequences (`src/models/`), or build `S` in closed form
   from a block model (`analysis/utils/block_model.py`).
2. Build the full similarity matrix `M` and its reference Fiedler vector at `p = 1`.
3. For each `p`: sub-sample `M → S_k` (`K` bootstrap reps), take the Fiedler vector of
   `L(S_k)`, sign-align it to the reference, and average.
4. Round the averaged vector to a bipartition and score it against the truth.

The **rounding rule** matters and is a paper choice: `kmeans` (k=2 on the Fiedler entries
of the normalized Laplacian `L_sym`) is `PAPER_METHOD`. `sign` (split at τ=0) and
`sigma2` (σ₂-gap search) remain available — `sigma2` is the expensive one and Figure 9
compares all three. Scoring is **NMI**, with `p̂*` the smallest grid `p` reaching
NMI ≥ **0.90** (not 0.95: at 0.95 the read-off lands on the flat top of the transition
where one noisy sample moves `p̂*` a whole grid step).

## Sampling methods (`src/core/sampling/`)

- **`uniform`** — baseline; `O(p·m²)` upper-triangular entries, sampled uniformly.
- **`lds`** — Leveraged Debiased Sampler: two-phase leverage-biased sampling with an
  unbiased IPW estimator. 10–100× faster than IALM. **Debias with the effective
  inclusion probability `π_ij = p₀ + (1−p₀)·p_ij`, not `p_ij`** — see `CHANGELOG.md`
  2026-02-20.
- **`leveraged`** — the same two-phase sampling followed by IALM matrix completion.
  Highest accuracy, slow.

## Performance constraints — do not regress

The November-2025 optimizations are why `m = 8192` is tractable:
partial eigendecomposition via `scipy.linalg.eigh` (not full SVD); vectorized boolean
masks for sub-sampling (not `np.random.choice` over `m²`); Laplacians computed once per
`p`, not per bootstrap rep; M-based metrics cached before the `p`-loop.

## Tests

There are none. `.gitignore` carried a bare `test*.py` pattern, which silently ignored
anything added under `tests/` — fixed, but no suite has been written. The only
self-verifying code is `python -m src.cache_io`. See `tests/README.md`.
