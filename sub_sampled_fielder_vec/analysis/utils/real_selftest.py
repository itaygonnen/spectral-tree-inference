"""Smoke check for the real-cohort pipeline: every call site, no heavy compute.

    python -m analysis.utils.real_selftest

Exists because a refactor dropped a parameter from ``sweep_cache_dir`` and the stale
call site in the interactive menu only failed *after* the user had answered eight
prompts and started a run. This walks the same code paths in a second: the cohort
loader, both cache paths, the config/CSV/JSON/plot writers, and a two-tree, two-p sweep
of the smallest cohort in a scratch results root. Nothing it does touches
``results/real_data``.

Exits non-zero on the first failure, so it is usable in a pre-run hook or CI.
"""
from __future__ import annotations

import contextlib
import io
import os
import sys
import tempfile
import traceback
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


@check
def cohorts_and_paths() -> str:
    from analysis.utils.real_cohorts import (cache_dir, get_cohort, list_cohorts,
                                             new_run_dir, screen_cache_path,
                                             sweep_cache_dir)
    cohorts = list_cohorts()
    if not cohorts:
        return "no cohorts on disk (data/cohorts/ empty) -- path checks only"
    c = cohorts[0]
    get_cohort(c.name)
    for fn in (cache_dir, screen_cache_path, sweep_cache_dir):
        fn(c.name)                            # one positional arg, and only one
    new_run_dir("selftest")
    m, seq_len = c.shape()
    return f"{len(cohorts)} cohort(s); {c.name}: {len(c.ids())} trees, m={m}, L={seq_len}"


@check
def exports_on_an_empty_run() -> str:
    """Every writer must cope with a run that produced nothing."""
    from analysis.utils.real_cohorts import list_cohorts, new_run_dir
    from analysis.utils.real_results import export_run, write_config
    names = [c.name for c in list_cohorts()][:1]
    run_dir = new_run_dir("selftest empty")
    write_config(run_dir, dict(reps=1, p_points=2), names, "sweep")
    files = export_run(run_dir, {n: [] for n in names})
    return f"{len(files) + 1} file(s) written for an empty run"


@check
def screen_and_sweep_two_trees() -> str:
    """The real thing, at the smallest size the data allows."""
    import numpy as np
    from analysis.utils.real_cohorts import (list_cohorts, screen_cache_path,
                                             sweep_cache_dir)
    from analysis.utils.real_eta_screen import run_real_eta_screen
    from analysis.utils.real_recovery_sweep import run_sweep, sweep_meta
    from analysis.utils.real_results import export_run

    cohorts = sorted(list_cohorts(), key=lambda c: c.shape()[0])
    if not cohorts:
        return "skipped: no cohorts on disk"
    c = cohorts[0]
    ids = c.ids(2)
    m, _ = c.shape()

    rows = run_real_eta_screen(ids, screen_cache_path(c.name),
                               cohort_name=c.name, workers=1)
    p_values = np.array([0.05, 1.0])
    run_sweep(ids, sweep_cache_dir(c.name), p_values, reps=1,
              cohort_name=c.name, m=m)
    sweep_meta(1, 10, 5, m)                   # keyword-free signature still valid

    from analysis.utils.real_cohorts import new_run_dir
    run_dir = new_run_dir("selftest sweep")
    files = export_run(run_dir, {c.name: ids})
    names = {f.name for f in files}
    for required in ("screening.csv", "curves.csv", "per_tree.csv", "summary.json"):
        if required not in names:
            raise AssertionError(f"{required} not written (got {sorted(names)})")
    return f"{c.name}: {len(rows)} screened, {len(ids)} swept, {len(files)} files exported"


@check
def menu_config_helpers() -> str:
    """The interactive layer's pure helpers, which the prompts sit on top of."""
    from analysis.utils.real_cohorts import list_cohorts
    from analysis.utils.real_interactive import (RULES, _default_rule, _plan,
                                                 _select_ids, _verdicts)
    cohorts = list_cohorts()[:2]
    if not cohorts:
        return "skipped: no cohorts on disk"
    verdicts_by = {c.name: _verdicts(c) for c in cohorts}
    rule = _default_rule(cohorts, verdicts_by)
    # the same keys _ask_config produces, so a missing one shows up here not mid-run
    from analysis.utils.real_interactive import _ask_config
    import inspect
    cfg = dict(max_trees=1, workers=1, rule=rule, p_min=0.01, p_points=2, reps=1,
               prefix="selftest", display_mode="progress", max_eta=20.0)
    src = inspect.getsource(_ask_config)
    missing = [k for k in ("max_trees", "workers", "rule", "p_min", "p_points", "reps",
                           "prefix", "display_mode", "max_eta")
               if f'"{k}"' not in src and f"'{k}'" not in src]
    if missing:
        raise AssertionError(f"_ask_config no longer sets {missing}")
    for is_screen in (True, False):
        rows = _plan(cohorts, cfg, is_screen, verdicts_by)
        if len(rows) != len(cohorts):
            raise AssertionError("plan lost a cohort")
    for name in RULES:
        _select_ids(cohorts[0].ids(3), verdicts_by[cohorts[0].name], name)
    return f"default gate: {rule!r}"


def main() -> None:
    # a scratch results root: the smoke check must never write into the real one
    tmp = tempfile.mkdtemp(prefix="real_selftest_")
    os.environ["STR_RESULTS_DIR"] = tmp
    print(f"scratch results root: {tmp}\n")

    failed = 0
    for fn in CHECKS:
        # the checks run real screens and sweeps; their progress output is noise here,
        # so it is captured and only replayed when something fails
        buf = io.StringIO()
        try:
            # stderr too: that is where tqdm draws, and a half-drawn bar between two
            # "ok" lines is exactly the noise this check should not add
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                note = fn()
            print(f"  ok    {fn.__name__}: {note}")
        except Exception:
            failed += 1
            print(f"  FAIL  {fn.__name__}")
            out = buf.getvalue().strip()
            if out:
                print("        --- output before the failure ---")
                for line in out.splitlines()[-8:]:
                    print(f"        {line}")
            traceback.print_exc(file=sys.stdout)
    print()
    if failed:
        print(f"{failed} of {len(CHECKS)} checks failed")
        sys.exit(1)
    print(f"all {len(CHECKS)} checks passed")


if __name__ == "__main__":
    main()
