"""Compute normalised diagnostic metrics for an existing NJ sweep — no rerun.

Motivation: the raw ``||D - D̂||_2`` metric in the Figure-8 overlay grows
linearly with n at fixed p, which is a dimensional artifact of the spectral
norm (and of mean-imputation's all-ones offset), not a sub-sampling failure.
The user's intuition — "at fixed p, larger n samples more entries and should
give a more accurate estimate" — is correct for *per-entry* error. To
visualise that, we add three scale-invariant variants of the same metric:

    linf          = ||D - D̂||_∞                  (max absolute entry error)
    rel_specnorm  = ||D - D̂||_2 / ||D||_2        (relative spec-norm error)
    rel_frob      = ||D - D̂||_F / ||D||_F        (relative Frobenius error)

This script regenerates D and D̂ from the cfg seeds recorded in the existing
``nj_meta.json``; it does NOT rerun NJ. Output lands in a sibling
``nj_meta_extra.json`` per (n, L) directory. The loader at
``src/utils/nj_io.py`` opportunistically merges the extras.

Usage::

    python scripts/nj_recompute_normalized_metrics.py <sweep_dir>
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Dict, List

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from src.config.presets import custom_config
from src.core.similarity_builder import SimilarityMatrixBuilder
from src.runners.nj_sweep import _build_truth, _impute_mean
from src.utils.metrics import estimate_operator_norm_diff


def _rebuild_cfg(meta: Dict):
    return custom_config(
        num_taxa=int(meta["n_taxa"]),
        sequence_length=int(meta["seq_len"]),
        mutation_rate=float(meta.get("mutation_rate", 0.1)),
        tree_model=str(meta.get("tree_model", "balanced_binary")),
        seq_model=str(meta.get("seq_model", "JC69")),
        p_values=[float(p) for p in meta["p_values"]],
        bootstrap_reps=int(meta["bootstrap_reps"]),
        sampling_method="uniform",
        matrix_kind="distance",
    )


def _recompute_for_dir(sub: Path) -> Path:
    meta_path = sub / "nj_meta.json"
    with open(meta_path) as f:
        meta = json.load(f)

    n_taxa = int(meta["n_taxa"])
    seq_len = int(meta["seq_len"])
    cfg = _rebuild_cfg(meta)
    imputation = str(meta.get("imputation", "mean"))

    # Best-effort replay: seed numpy so tree + sequence simulation are
    # reproducible run-to-run (the original sweep didn't seed explicitly,
    # so this won't bit-match the on-disk run, but the statistical claim
    # is unchanged).
    np.random.seed(int(cfg.experiment.seed))

    print(f"[recompute] n={n_taxa} L={seq_len} (rebuilding ground truth)...", flush=True)
    t0 = time.time()
    tree, observations, D, taxa_metadata = _build_truth(cfg, n_taxa, seq_len)
    # ``estimate_operator_norm_diff`` accepts (S, M); pass (D, zeros) by using
    # the dedicated np.linalg.norm path for ||D||_2 (one call, large n OK
    # because we do it once per directory).
    print(f"  built D in {time.time()-t0:.1f}s; computing ||D||...", flush=True)
    t0 = time.time()
    D_norm_2 = float(np.linalg.norm(D, ord=2))
    D_norm_F = float(np.linalg.norm(D, ord="fro"))
    D_norm_inf = float(np.max(np.abs(D)))
    print(f"  ||D||_2 = {D_norm_2:.4g}   ||D||_F = {D_norm_F:.4g}   "
          f"||D||_inf = {D_norm_inf:.4g}  ({time.time()-t0:.1f}s)", flush=True)

    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="distance")
    sampler = builder.sampler

    p_values: List[float] = [float(p) for p in meta["p_values"]]
    bootstrap_reps: int = int(meta["bootstrap_reps"])
    seed_base: int = int(cfg.experiment.seed)

    per_p_extra: List[Dict] = []
    for p_idx, p in enumerate(p_values):
        t0 = time.time()
        linf_vals: List[float] = []
        rel_spec_vals: List[float] = []
        rel_frob_vals: List[float] = []
        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            D_hat = sampler.sample(D, p, seed=seed)
            if imputation == "mean":
                D_hat = _impute_mean(D_hat)
            diff = D - D_hat
            linf_vals.append(float(np.max(np.abs(diff))))
            spec_diff = float(estimate_operator_norm_diff(D_hat, D, n_iter=20))
            rel_spec_vals.append(spec_diff / max(D_norm_2, 1e-30))
            rel_frob_vals.append(float(np.linalg.norm(diff, ord="fro"))
                                  / max(D_norm_F, 1e-30))
        per_p_extra.append({
            "p": float(p),
            "linf_mean": float(np.mean(linf_vals)),
            "linf_std": float(np.std(linf_vals)),
            "linf_per_rep": linf_vals,
            "rel_specnorm_mean": float(np.mean(rel_spec_vals)),
            "rel_specnorm_std": float(np.std(rel_spec_vals)),
            "rel_specnorm_per_rep": rel_spec_vals,
            "rel_frob_mean": float(np.mean(rel_frob_vals)),
            "rel_frob_std": float(np.std(rel_frob_vals)),
            "rel_frob_per_rep": rel_frob_vals,
        })
        print(f"  p={p:.4g}: linf={per_p_extra[-1]['linf_mean']:.4g}, "
              f"rel_spec={per_p_extra[-1]['rel_specnorm_mean']:.4g}, "
              f"rel_frob={per_p_extra[-1]['rel_frob_mean']:.4g}  "
              f"({time.time()-t0:.1f}s)", flush=True)

    extra = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": p_values,
        "bootstrap_reps": bootstrap_reps,
        "imputation": imputation,
        "D_norm_2": D_norm_2,
        "D_norm_F": D_norm_F,
        "D_norm_inf": D_norm_inf,
        "seed_base": seed_base,
        "note": ("Generated by nj_recompute_normalized_metrics.py. "
                 "Ground truth was rebuilt via best-effort seed replay; "
                 "values are statistically equivalent to the original sweep "
                 "but not bit-identical."),
        "per_p": per_p_extra,
    }
    out_path = sub / "nj_meta_extra.json"
    with open(out_path, "w") as f:
        json.dump(extra, f, indent=2)
    print(f"  wrote {out_path}", flush=True)
    return out_path


def main():
    if len(sys.argv) < 2:
        print("usage: nj_recompute_normalized_metrics.py <sweep_dir>",
              file=sys.stderr)
        sys.exit(1)
    sweep_dir = Path(sys.argv[1]).resolve()
    if not sweep_dir.is_dir():
        print(f"not a directory: {sweep_dir}", file=sys.stderr)
        sys.exit(2)
    subs = [d for d in sorted(sweep_dir.iterdir())
            if d.is_dir() and (d / "nj_meta.json").exists()]
    if not subs:
        print(f"no nj_meta.json subdirs under {sweep_dir}", file=sys.stderr)
        sys.exit(3)
    # Process smallest n first so failures surface fast.
    subs = sorted(subs, key=lambda d: int(json.load(open(d / "nj_meta.json"))
                                          ["n_taxa"]))
    written = []
    for sub in subs:
        try:
            written.append(_recompute_for_dir(sub))
        except Exception as e:
            print(f"  ERROR in {sub}: {e!r}", file=sys.stderr, flush=True)
    print("\nDone. Wrote:")
    for p in written:
        print(f"  {p}")


if __name__ == "__main__":
    main()
