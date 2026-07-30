#!/usr/bin/env python3
"""Migrate legacy caches and notebook-output dirs into the unified cache_io layout.

Source → Destination mapping
============================

    cache/theoretical_interpretation/full/<key>/
        → cache/full_matrix/<key>/
        (rename fiedler_full.npz → fiedler.npz, inner array key too)

    cache/theoretical_interpretation/subsampled/<fk>/p..seed../
        → cache/sweep_trial/<fk>/p..seed../
        (rewrite agreement.json → metrics.json={"sign_agreement": v})

    src/cache/<old_key>/                      (persistent_cache)
        → cache/experiment_data/<new_key>/
        (new_key is computed from metadata.json via cache_io.make_key)

    src/cache/eta_pool/<pkey>/eta<TT>/sample_NNNN/
        → cache/pool_sample/<pkey>/eta<TT>/sample_NNNN/

    src/cache/eta_pool/<pkey>/eta<TT>/sample_NNNN/sweeps/<skey>/
        → cache/bootstrap_sweep/<pkey>__<eta..>__<sample..>__<skey>/

    analysis/<sub>/<nb>_outputs/
        → results/notebooks/<sub>/<nb>/

Idempotent: re-running skips already-complete destinations. Flags:
    --dry-run       print the plan, do nothing
    --verify        re-read every migrated entry to confirm sentinels
    --prune-empty   remove empty legacy directories afterwards
    --quiet         suppress per-entry logs
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve()
PROJECT_ROOT = HERE.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.cache_io import (  # noqa: E402
    CACHE_ROOT,
    RESULTS_ROOT,
    SENTINEL,
    make_key,
    full_matrix,
    experiment_data,
    pool_sample,
    sweep_trial,
    bootstrap_sweep,
)


# ----- migration helpers ---------------------------------------------------


def _make_logger(verbose: bool) -> Callable[[str], None]:
    return print if verbose else (lambda _msg: None)


def _safe_move(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))


def _new_persistent_key(meta: dict[str, Any]) -> str:
    """Compute the new experiment_data key from a legacy metadata.json."""
    kwargs: dict[str, Any] = {
        "n":    int(meta["n_taxa"]),
        "L":    int(meta["seq_len"]),
        "mu":   float(meta["mutation_rate"]),
        "tree": str(meta["tree_model"]),
        "seq":  str(meta["seq_model"]),
    }
    for k, v in (meta.get("tree_params") or {}).items():
        if k == "num_taxa":
            continue
        kwargs[k] = v
    for k, v in (meta.get("seq_params") or {}).items():
        if k == "mutation_rate":
            continue
        kwargs[k] = v
    return make_key(**kwargs)


# ----- migrators -----------------------------------------------------------


def migrate_theoretical_full(*, dry_run: bool, log: Callable[[str], None]) -> int:
    src_root = PROJECT_ROOT / "cache" / "theoretical_interpretation" / "full"
    if not src_root.exists():
        return 0
    log(f"[theoretical_full] {src_root}")
    moved = 0
    for src in sorted(src_root.iterdir()):
        if not src.is_dir() or not (src / SENTINEL).exists():
            continue
        dst = CACHE_ROOT / "full_matrix" / src.name
        if (dst / SENTINEL).exists():
            log(f"  skip (dst complete): full_matrix/{src.name}")
            continue
        if dry_run:
            log(f"  [DRY] move full/{src.name} -> full_matrix/{src.name} "
                f"+ rename fiedler_full -> fiedler")
            moved += 1
            continue
        _safe_move(src, dst)
        old_npz = dst / "fiedler_full.npz"
        new_npz = dst / "fiedler.npz"
        if old_npz.exists():
            arr = np.load(old_npz)["fiedler_full"]
            np.savez_compressed(new_npz, fiedler=arr)
            old_npz.unlink()
            (dst / SENTINEL).touch()
        moved += 1
        log(f"  moved full_matrix/{src.name}")
    log(f"  -> {moved} entries")
    return moved


def migrate_theoretical_subsampled(*, dry_run: bool, log: Callable[[str], None]) -> int:
    src_root = PROJECT_ROOT / "cache" / "theoretical_interpretation" / "subsampled"
    if not src_root.exists():
        return 0
    log(f"[theoretical_subsampled] {src_root}")
    moved = 0
    for src_fk in sorted(src_root.iterdir()):
        if not src_fk.is_dir():
            continue
        for src_trial in sorted(src_fk.iterdir()):
            if not src_trial.is_dir() or not (src_trial / SENTINEL).exists():
                continue
            dst_trial = CACHE_ROOT / "sweep_trial" / src_fk.name / src_trial.name
            if (dst_trial / SENTINEL).exists():
                continue
            if dry_run:
                log(f"  [DRY] move subsampled/{src_fk.name}/{src_trial.name} "
                    f"-> sweep_trial; rewrite agreement.json -> metrics.json")
                moved += 1
                continue
            _safe_move(src_trial, dst_trial)
            old = dst_trial / "agreement.json"
            if old.exists():
                payload = json.loads(old.read_text())
                metrics = {"sign_agreement": float(payload["agreement"])}
                (dst_trial / "metrics.json").write_text(json.dumps(metrics, indent=2))
                old.unlink()
                (dst_trial / SENTINEL).touch()
            moved += 1
    log(f"  -> {moved} trials")
    return moved


def migrate_persistent_cache(*, dry_run: bool, log: Callable[[str], None]) -> int:
    src_root = PROJECT_ROOT / "src" / "cache"
    if not src_root.exists():
        return 0
    log(f"[persistent_cache] {src_root}")
    moved = 0
    skipped = 0
    for src in sorted(src_root.iterdir()):
        if not src.is_dir() or src.name == "eta_pool":
            continue
        if not (src / SENTINEL).exists():
            continue
        meta_path = src / "metadata.json"
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text())
                new_key = _new_persistent_key(meta)
            except Exception as e:
                log(f"  warn: metadata parse failed ({e}); preserving key {src.name}")
                new_key = src.name
        else:
            log(f"  warn: no metadata.json; preserving key {src.name}")
            new_key = src.name
        dst = CACHE_ROOT / "experiment_data" / new_key
        if (dst / SENTINEL).exists():
            log(f"  skip (dst complete): experiment_data/{new_key}")
            skipped += 1
            continue
        if dry_run:
            log(f"  [DRY] move src/cache/{src.name} -> experiment_data/{new_key}")
            moved += 1
            continue
        _safe_move(src, dst)
        moved += 1
        log(f"  moved experiment_data/{new_key}")
    log(f"  -> {moved} entries (skipped {skipped})")
    return moved


def migrate_eta_pool(*, dry_run: bool, log: Callable[[str], None]) -> tuple[int, int]:
    src_root = PROJECT_ROOT / "src" / "cache" / "eta_pool"
    if not src_root.exists():
        return 0, 0
    log(f"[eta_pool] {src_root}")
    pool_moved = 0
    sweeps_moved = 0
    for src_pk in sorted(src_root.iterdir()):
        if not src_pk.is_dir():
            continue
        for src_eta in sorted(src_pk.iterdir()):
            if not src_eta.is_dir() or not src_eta.name.startswith("eta"):
                continue
            for src_sample in sorted(src_eta.iterdir()):
                if not src_sample.is_dir() or not src_sample.name.startswith("sample_"):
                    continue

                # 1) Extract sweeps/<skey>/ children into bootstrap_sweep.
                src_sweeps = src_sample / "sweeps"
                if src_sweeps.exists():
                    for src_sw in sorted(src_sweeps.iterdir()):
                        if not src_sw.is_dir() or not (src_sw / SENTINEL).exists():
                            continue
                        sweep_id = (
                            f"{src_pk.name}__{src_eta.name}"
                            f"__{src_sample.name}__{src_sw.name}"
                        )
                        dst_sw = CACHE_ROOT / "bootstrap_sweep" / sweep_id
                        if (dst_sw / SENTINEL).exists():
                            continue
                        if dry_run:
                            log(f"  [DRY] move sweep {src_sw.relative_to(PROJECT_ROOT)} "
                                f"-> bootstrap_sweep/{sweep_id}")
                        else:
                            _safe_move(src_sw, dst_sw)
                        sweeps_moved += 1
                    if not dry_run and src_sweeps.exists() and not any(src_sweeps.iterdir()):
                        src_sweeps.rmdir()

                # 2) Move the sample dir itself.
                if not (src_sample / SENTINEL).exists():
                    continue
                dst_sample = (
                    CACHE_ROOT / "pool_sample"
                    / src_pk.name / src_eta.name / src_sample.name
                )
                if (dst_sample / SENTINEL).exists():
                    continue
                if dry_run:
                    log(f"  [DRY] move pool {src_sample.relative_to(PROJECT_ROOT)} "
                        f"-> pool_sample/{src_pk.name}/{src_eta.name}/{src_sample.name}")
                else:
                    _safe_move(src_sample, dst_sample)
                pool_moved += 1

        # 3) Preserve the per-param manifest.json next to the bins.
        old_manifest = src_pk / "manifest.json"
        if old_manifest.exists():
            new_manifest = CACHE_ROOT / "pool_sample" / src_pk.name / "manifest.json"
            if not new_manifest.exists():
                if dry_run:
                    log(f"  [DRY] move manifest {old_manifest.relative_to(PROJECT_ROOT)} "
                        f"-> pool_sample/{src_pk.name}/manifest.json")
                else:
                    _safe_move(old_manifest, new_manifest)

    log(f"  -> {pool_moved} pool samples + {sweeps_moved} sweeps")
    return pool_moved, sweeps_moved


def migrate_results_to_runs(*, dry_run: bool, log: Callable[[str], None]) -> int:
    """Move every production run dir at ``results/`` top-level under ``results/runs/``.

    Two source layouts are handled:
      * flat: ``results/<ts>-<run_name>/``                   → ``results/runs/<ts>-<run_name>/``
      * nested: ``results/<tree_model>/<sampling_method>/<ts>-<run_name>/`` →
                ``results/runs/<tree_model>/<sampling_method>/<ts>-<run_name>/``

    The whole tree-model subtree (including its non-timestamp children) is
    moved in one go for the nested case, so loose files like
    ``balanced_binary/comparison_distance_vs_similarity_n256.png`` land
    under ``runs/balanced_binary/...`` too. ``results/notebooks/`` and
    ``results/runs/`` are left alone.
    """
    src_root = PROJECT_ROOT / "results"
    if not src_root.exists():
        return 0
    log(f"[results_runs] {src_root}")
    moved = 0
    runs_root = src_root / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)

    # Anything at top level that isn't "notebooks" or "runs" is a production run.
    for entry in sorted(src_root.iterdir()):
        if entry.name in ("notebooks", "runs"):
            continue
        if not entry.is_dir() and not entry.is_file():
            continue
        dst = runs_root / entry.name
        if dst.exists():
            log(f"  skip (dst exists): runs/{entry.name}")
            continue
        if dry_run:
            log(f"  [DRY] move results/{entry.name} -> results/runs/{entry.name}")
        else:
            _safe_move(entry, dst)
            log(f"  moved runs/{entry.name}")
        moved += 1
    log(f"  -> {moved} entries")
    return moved


def migrate_notebook_outputs(*, dry_run: bool, log: Callable[[str], None]) -> int:
    src_root = PROJECT_ROOT / "analysis" / "theoretical_interpretation"
    if not src_root.exists():
        return 0
    log(f"[notebook_outputs] {src_root}")
    moved = 0
    for sub in sorted(src_root.iterdir()):
        if not sub.is_dir() or sub.name in ("utils", "__pycache__"):
            continue
        for child in sorted(sub.iterdir()):
            if not child.is_dir() or not child.name.endswith("_outputs"):
                continue
            notebook_stem = child.name[: -len("_outputs")]
            dst = RESULTS_ROOT / "notebooks" / sub.name / notebook_stem
            if dst.exists():
                continue
            if dry_run:
                log(f"  [DRY] move {child.relative_to(PROJECT_ROOT)} "
                    f"-> results/notebooks/{sub.name}/{notebook_stem}")
            else:
                _safe_move(child, dst)
            moved += 1
    log(f"  -> {moved} notebook output dirs")
    return moved


# ----- verify + prune ------------------------------------------------------


def verify(*, log: Callable[[str], None]) -> int:
    log("Verifying migrated entries...")
    scopes = {
        "full_matrix":     full_matrix,
        "experiment_data": experiment_data,
        "pool_sample":     pool_sample,
        "sweep_trial":     sweep_trial,
        "bootstrap_sweep": bootstrap_sweep,
    }
    errors = 0
    for name, scope in scopes.items():
        if not scope.root.exists():
            log(f"  [{name}] empty (ok)")
            continue
        count = 0
        for sentinel in scope.root.rglob(SENTINEL):
            d = sentinel.parent
            try:
                _ = scope.load(*d.relative_to(scope.root).parts)
                count += 1
            except Exception as e:
                log(f"  ERROR [{name}] {d}: {e}")
                errors += 1
        log(f"  [{name}] {count} entries OK")
    return errors


def prune_empty(*, dry_run: bool, log: Callable[[str], None]) -> int:
    candidates = [
        PROJECT_ROOT / "cache" / "theoretical_interpretation" / "full",
        PROJECT_ROOT / "cache" / "theoretical_interpretation" / "subsampled",
        PROJECT_ROOT / "cache" / "theoretical_interpretation",
        PROJECT_ROOT / "src" / "cache" / "eta_pool",
    ]
    removed = 0
    for c in candidates:
        if not c.exists():
            continue
        has_file = any(p.is_file() for p in c.rglob("*"))
        if has_file:
            log(f"  skip non-empty: {c.relative_to(PROJECT_ROOT)}")
            continue
        if dry_run:
            log(f"  [DRY] remove empty: {c.relative_to(PROJECT_ROOT)}")
        else:
            shutil.rmtree(c)
            log(f"  removed empty: {c.relative_to(PROJECT_ROOT)}")
        removed += 1
    return removed


# ----- entry point ---------------------------------------------------------


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--prune-empty", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    log = _make_logger(verbose=not args.quiet)
    total = 0
    total += migrate_theoretical_full(dry_run=args.dry_run, log=log)
    total += migrate_theoretical_subsampled(dry_run=args.dry_run, log=log)
    total += migrate_persistent_cache(dry_run=args.dry_run, log=log)
    pool_n, sweeps_n = migrate_eta_pool(dry_run=args.dry_run, log=log)
    total += pool_n + sweeps_n
    total += migrate_notebook_outputs(dry_run=args.dry_run, log=log)
    total += migrate_results_to_runs(dry_run=args.dry_run, log=log)

    suffix = " (DRY RUN)" if args.dry_run else ""
    print(f"\nMigration: {total} entries moved{suffix}")

    if args.verify:
        if args.dry_run:
            print("\n(verify skipped: --dry-run)")
        else:
            errs = verify(log=log)
            if errs:
                print(f"\nVerify: {errs} ERRORS")
                return 1
            print("Verify: ok")

    if args.prune_empty:
        print("\nPruning empty legacy dirs...")
        prune_empty(dry_run=args.dry_run, log=log)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
