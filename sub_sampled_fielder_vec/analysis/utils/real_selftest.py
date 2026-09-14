"""Smoke check for the screen/sweep pipeline: every call site, no heavy compute.

    python -m analysis.utils.real_selftest

Exists because a refactor dropped a parameter from ``sweep_cache_dir`` and the stale
call site in the interactive menu only failed *after* the user had answered eight
prompts and started a run. This walks the same code paths in a second: the dataset
loader, both cache paths, the config/CSV/JSON/plot writers, and a two-tree, two-p sweep
of the smallest dataset in a scratch results root. Nothing it does touches
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
def datasets_and_paths() -> str:
    from analysis.utils.real_datasets import get_dataset, list_datasets
    from src.utils.run_paths import (cache_dir, new_run_dir, screen_cache_path,
                                     sweep_cache_dir)
    datasets = list_datasets()
    if not datasets:
        return "no datasets on disk (data/datasets/ empty) -- path checks only"
    c = datasets[0]
    get_dataset(c.name)
    for fn in (cache_dir, screen_cache_path, sweep_cache_dir):
        fn(c.name)                            # one positional arg, and only one
    new_run_dir("selftest")
    m, seq_len = c.shape()
    return f"{len(datasets)} dataset(s); {c.name}: {len(c.ids())} trees, m={m}, L={seq_len}"


@check
def exports_on_an_empty_run() -> str:
    """Every writer must cope with a run that produced nothing."""
    from analysis.utils.real_datasets import list_datasets
    from src.utils.run_export import export_run, write_config
    from src.utils.run_paths import new_run_dir
    names = [c.name for c in list_datasets()][:1]
    run_dir = new_run_dir("selftest empty")
    write_config(run_dir, dict(reps=1, p_points=2), names, "sweep")
    files = export_run(run_dir, {n: [] for n in names})
    return f"{len(files) + 1} file(s) written for an empty run"


@check
def screen_and_sweep_two_trees() -> str:
    """The real thing, at the smallest size the data allows."""
    import numpy as np
    from analysis.utils.real_datasets import list_datasets
    from src.runners.operator_screen import run_screen
    from src.runners.operator_sweep import run_sweep, sweep_meta
    from src.utils.run_export import export_run
    from src.utils.run_paths import screen_cache_path, sweep_cache_dir

    datasets = sorted(list_datasets(), key=lambda c: c.shape()[0])
    if not datasets:
        return "skipped: no datasets on disk"
    c = datasets[0]
    ids = c.ids(2)
    m, _ = c.shape()

    loader = c.loader()
    rows = run_screen(ids, screen_cache_path(c.name), loader,
                      source=c.name, workers=1)
    p_values = np.array([0.05, 1.0])
    run_sweep(ids, sweep_cache_dir(c.name), loader, p_values, reps=1,
              source=c.name, m=m)
    sweep_meta(1, 10, 5, m)                   # keyword-free signature still valid

    from src.utils.run_paths import new_run_dir
    run_dir = new_run_dir("selftest sweep")
    files = export_run(run_dir, {c.name: ids})
    names = {f.name for f in files}
    for required in ("screening.csv", "curves.csv", "per_tree.csv", "summary.json"):
        if required not in names:
            raise AssertionError(f"{required} not written (got {sorted(names)})")
    return f"{c.name}: {len(rows)} screened, {len(ids)} swept, {len(files)} files exported"


@check
def menu_config_helpers() -> str:
    """The shared runner's pure helpers, which both front-ends sit on top of."""
    import dataclasses
    import inspect

    from analysis.utils.real_datasets import list_datasets
    from analysis.utils.real_interactive import _ask_config
    from src.runners.experiment_run import (GATES, RunSpec, default_gate, gate_key,
                                            plan, select_ids, verdicts)
    datasets = list_datasets()[:2]
    if not datasets:
        return "skipped: no datasets on disk"
    names = [c.name for c in datasets]
    verdicts_by = {n: verdicts(n) for n in names}
    gate = default_gate(names, verdicts_by)

    # every RunSpec field the prompts are supposed to set must appear in _ask_config,
    # so a field added to the spec and forgotten in the menu shows up here, not mid-run
    src = inspect.getsource(_ask_config)
    asked = {f.name for f in dataclasses.fields(RunSpec)} - {"stage", "num_gaps",
                                                             "min_split"}
    missing = [k for k in sorted(asked) if k not in src]
    if missing:
        raise AssertionError(f"_ask_config no longer sets {missing}")

    sources = [c.source(1) for c in datasets]
    for stage in ("screen", "sweep"):
        spec = RunSpec(stage=stage, max_trees=1, workers=1, gate=gate,
                       p_min=0.01, p_points=2, reps=1, prefix="selftest")
        rows = plan(spec, sources, verdicts_by)
        if len(rows) != len(datasets):
            raise AssertionError("plan lost a dataset")
    for name in GATES:
        select_ids(datasets[0].ids(3), verdicts_by[names[0]], name)
    if gate_key("valid_S") != "valid_L":
        raise AssertionError("the deprecated --dataset-rule valid_S alias is broken")
    return f"default gate: {gate!r}"


@check
def both_front_ends_agree() -> str:
    """The menu and the command line must build the SAME run out of the same answers.

    This is the check the merge exists for: the two used to carry a selection
    implementation each, and wrote different vocabularies into the same summary.json
    field. Anything that drifts again fails here rather than on the professor's cluster.
    """
    import dataclasses
    import importlib.util

    from analysis.utils.real_datasets import list_datasets
    from src.runners.experiment_run import RunSpec, plan, verdicts
    datasets = list_datasets()[:1]
    if not datasets:
        return "skipped: no datasets on disk"
    name = datasets[0].name

    path = _ROOT / "scripts" / "run_sweep.py"
    loader = importlib.util.spec_from_file_location("_run_real_sweep", path)
    cli = importlib.util.module_from_spec(loader)
    loader.loader.exec_module(cli)

    # the command line the professor would type (with the deprecated gate alias), and
    # the menu's equivalent answers
    args = cli.build_parser().parse_args(
        ["--dataset", name, "--stage", "sweep", "--limit", "3",
         "--dataset-rule", "valid_S", "--max-eta", "20", "--p-points", "2",
         "--reps", "1"])
    from_cli = cli.spec_from_args(args)
    from_menu = RunSpec(stage="sweep", max_trees=3, gate="valid_L", max_eta=20.0,
                        p_points=2, reps=1, workers=from_cli.workers,
                        p_min=from_cli.p_min, operators=from_cli.operators)
    if dataclasses.asdict(from_cli) != dataclasses.asdict(from_menu):
        diff = {k: (v, dataclasses.asdict(from_menu)[k])
                for k, v in dataclasses.asdict(from_cli).items()
                if v != dataclasses.asdict(from_menu)[k]}
        raise AssertionError(f"the two front-ends build different runs: {diff}")
    v = {name: verdicts(name)}
    sources = cli.sources_from_args(args)
    a, b = plan(from_cli, sources, v), plan(from_menu, sources, v)
    if [r.selected for r in a] != [r.selected for r in b]:
        raise AssertionError("the two front-ends select different trees")
    if list(from_cli.p_values()) != list(from_menu.p_values()):
        raise AssertionError("the two front-ends build different p-grids")
    return f"{name}: same gate, same {len(a[0].selected)} trees, same grid"


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
