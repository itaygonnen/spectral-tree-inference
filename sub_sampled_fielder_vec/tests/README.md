# tests/

**There is no test suite.** This is a statement of fact, not a placeholder.

`tests/__init__.py` describes `unit/`, `integration/` and `performance/` subdirectories
that have never existed.

## Why nothing is here

The repo-root `.gitignore` carried a bare `test*.py` pattern, which matched
`sub_sampled_fielder_vec/tests/test_*.py` at any depth. Any test file added here was
silently ignored by git and could never be committed — so the directory could not fill up
even if someone tried.

That pattern is now negated for this directory (`!sub_sampled_fielder_vec/tests/test_*.py`),
verified with `git check-ignore`. Adding `tests/test_foo.py` will work.

## What verification does exist

- `python -m src.cache_io` — the one piece of self-verifying code in ~21,000 lines.
  Asserts key determinism, `.complete` sentinel discipline, save/load round-trip,
  `get_or_compute` cache hits, and `extend` merge behaviour.
- `python scripts/sync_paper_figures.py --check` — validates figure provenance against
  the manifest and cross-checks it against the paper's `\includegraphics` calls.
- `python scripts/collate_open_items.py` — lints the open-items register format.
- `scripts/validation/{fiedler_plateau_validator,mutation_rate_validation}.py` —
  figure-generating sanity checks, **not** assertions.
- `examples/*.py` — demo scripts, no assertions.

There is no `pytest.ini`, `pyproject.toml`, `conftest.py`, or CI config anywhere in the
package.

## If you write tests

Start where a silent failure would be most expensive:

1. **Cache key stability.** Three independent key schemes coexist
   (`cache_io.make_key`, `eta_pool_cache.param_key`, `sweep_cache.compute_sweep_key`) and
   a change to any of them orphans gigabytes without erroring. Pin the canonical keys:
   `sign=acb9d64195e6`, `sigma2=6f5687277cea`, `kmeans=234a0bd360f7`.
2. **`bootstrap_p_sweep_simple` rounding rules.** `sign` / `sigma2` / `kmeans` on a
   hand-built matrix with a known bipartition.
3. **Metric invariances.** NMI and ARI must be invariant to a global Fiedler sign flip;
   `find_discrete_threshold` must return a `p` from the grid, never an interpolated one.

Do not write tests that need `cache/`, `results/` or `data/` — all three are gitignored,
so such a test passes only on this machine.
