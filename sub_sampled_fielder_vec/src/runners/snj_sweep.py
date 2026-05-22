"""SNJ-with-subsampling sweep over ``p`` for a fixed (n_taxa, seq_len).

Mirrors the shape of ``bootstrap_sweep.sweep_for_params`` but runs the
SNJ-paper metrics from Jaffe & Kluger:

  1. ``||R - R̂||_2``  (Thm 4.2 quantity)
  2. σ₂ separation     (Fig. 3 analog, computed at the initial SNJ step)
  3. Robinson-Foulds    (paper's empirical metric)

Sampling is plain uniform with no completion: unsampled entries are 0,
sampled entries are IPW-scaled by ``1/p`` (see
``src/core/sampling/uniform/sampler.py``).
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import dendropy
from dendropy.calculate import treecompare

import spectraltree

from ..config import StructuredConfig
from ..core.similarity_builder import SimilarityMatrixBuilder
from ..models import get_tree_factory, get_sequence_factory
from ..utils.logging import log_info
from ..utils.metrics import estimate_operator_norm_diff
from ..utils.snj_sigma2_separation import compute_sigma2_separation
from ..utils.summaries import save_json


def _build_truth(cfg: StructuredConfig, n_taxa: int, seq_len: int
                 ) -> Tuple[object, np.ndarray, np.ndarray, object]:
    """Generate tree + observations + full similarity R + TaxaMetadata."""
    tree_factory = get_tree_factory(cfg.tree.model, cfg.tree.params)
    tree = tree_factory()
    seq_factory = get_sequence_factory(cfg.sequence.model, cfg.sequence.params)
    seq_model = seq_factory()

    char_matrix, taxa_metadata = spectraltree.simulate_sequences(
        seq_len=seq_len,
        tree_model=tree,
        seq_model=seq_model,
        mutation_rate=cfg.sequence.params["mutation_rate"],
    )

    R = spectraltree.JC_similarity_matrix(char_matrix)
    # Numerical safety: clip to [0, 1] (paralinear formula can dip negative
    # for very short sequences).
    R = np.clip(R, 0.0, 1.0)
    np.fill_diagonal(R, 1.0)
    return tree, char_matrix, R, taxa_metadata


def _run_snj(R_hat: np.ndarray, taxa_metadata, bifurcating: bool, alpha: float):
    """Run SNJ on the (possibly subsampled+IPW-rescaled) similarity matrix."""
    snj = spectraltree.SpectralNeighborJoining(
        similarity_metric=spectraltree.JC_similarity_matrix,
        scorer=spectraltree.sv2,
        bifurcating=bifurcating,
        alpha=alpha,
    )
    return snj.reconstruct_from_similarity(R_hat, taxa_metadata)


def _rf_distance(ref_tree, inferred_tree, taxon_namespace) -> int:
    """Robinson-Foulds (symmetric difference) between two trees.

    Both trees must share `taxon_namespace`. We update bipartitions in-place
    on copies that live inside the same namespace.
    """
    # `symmetric_difference` requires bipartitions to be encoded.
    for tr in (ref_tree, inferred_tree):
        tr.encode_bipartitions()
    return int(treecompare.symmetric_difference(ref_tree, inferred_tree))


def snj_sweep_for_params(cfg: StructuredConfig, n_taxa: int, seq_len: int,
                         run_dir: str, *, bifurcating: bool = False,
                         snj_alpha: float = 1.0) -> Dict:
    """Run the SNJ-subsampling sweep at fixed (n_taxa, seq_len).

    Writes ``snj_results.json`` into ``run_dir`` and returns the same dict.
    """
    os.makedirs(run_dir, exist_ok=True)
    log_info("snj", f"Building ground truth (n={n_taxa}, L={seq_len})...", force=True)
    tree, observations, R, taxa_metadata = _build_truth(cfg, n_taxa, seq_len)
    # Encode bipartitions once on the reference so RF can compare against it.
    tree.encode_bipartitions()

    # We don't ask the SimilarityMatrixBuilder to re-build R from observations
    # since we already have it; instead we sample R directly via its sampler.
    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="similarity")
    sampler = builder.sampler

    p_values: List[float] = list(cfg.experiment.p_values)
    bootstrap_reps: int = int(cfg.experiment.bootstrap_reps)
    seed_base: int = int(cfg.experiment.seed)

    per_p: List[Dict] = []
    t0 = time.time()
    for p_idx, p in enumerate(p_values):
        log_info("snj", f"  p = {p:.4g} ({bootstrap_reps} reps)", force=True)
        spec_norms: List[float] = []
        rf_values: List[int] = []
        adj_all: List[np.ndarray] = []
        nonadj_all: List[np.ndarray] = []

        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            R_hat = sampler.sample(R, p, seed=seed)

            # Metric 1: spectral-norm error
            spec_norm = float(estimate_operator_norm_diff(R_hat, R, n_iter=20))
            spec_norms.append(spec_norm)

            # Metric 2: σ₂ separation at step 0
            adj, nonadj = compute_sigma2_separation(R_hat, tree, taxa_metadata)
            adj_all.append(adj)
            nonadj_all.append(nonadj)

            # Metric 3: RF distance via SNJ
            try:
                inferred = _run_snj(R_hat, taxa_metadata, bifurcating=bifurcating,
                                    alpha=snj_alpha)
                # Clone of reference tree under the *same* namespace, so RF works.
                rf = _rf_distance(
                    ref_tree=dendropy.Tree(tree),
                    inferred_tree=inferred,
                    taxon_namespace=taxa_metadata.taxon_namespace,
                )
            except Exception as e:
                # SNJ can in principle fail at very low p (e.g., R̂ is all-zero
                # off-diagonal, ties produce degenerate merges). Record as
                # the worst-case RF.
                log_info("snj", f"    rep {rep}: SNJ failed ({e!r}); RF=2(m-3)", force=True)
                rf = 2 * (n_taxa - 3)
            rf_values.append(int(rf))

        # Aggregate σ₂ separation across reps: store concatenated samples.
        adj_concat = np.concatenate(adj_all) if adj_all else np.array([], dtype=np.float64)
        nonadj_concat = np.concatenate(nonadj_all) if nonadj_all else np.array([], dtype=np.float64)

        per_p.append({
            "p": float(p),
            "spec_norm_mean": float(np.mean(spec_norms)),
            "spec_norm_std": float(np.std(spec_norms)),
            "spec_norm_per_rep": [float(x) for x in spec_norms],
            "rf_mean": float(np.mean(rf_values)),
            "rf_std": float(np.std(rf_values)),
            "rf_per_rep": [int(x) for x in rf_values],
            "rf_success_rate": float(np.mean([1.0 if r == 0 else 0.0 for r in rf_values])),
            "_sigma2_adj": adj_concat,        # arrays kept out of meta JSON
            "_sigma2_nonadj": nonadj_concat,
        })

    elapsed = time.time() - t0
    # Split layout: small meta.json + compressed .npz for σ₂ arrays.
    sigma2_adj_stack = np.stack([entry.pop("_sigma2_adj") for entry in per_p])
    sigma2_nonadj_stack = np.stack([entry.pop("_sigma2_nonadj") for entry in per_p])

    meta = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": [float(p) for p in p_values],
        "bootstrap_reps": bootstrap_reps,
        "tree_model": cfg.tree.model,
        "seq_model": cfg.sequence.model,
        "mutation_rate": float(cfg.sequence.params.get("mutation_rate", 1.0)),
        "sampling_method": "uniform_no_completion",
        "snj_bifurcating": bifurcating,
        "snj_alpha": snj_alpha,
        "per_p": per_p,
        "elapsed_seconds": elapsed,
    }
    run_path = Path(run_dir)
    meta_path = run_path / "snj_meta.json"
    npz_path = run_path / "snj_arrays.npz"
    save_json(meta, str(meta_path))
    np.savez_compressed(npz_path,
                        sigma2_adj=sigma2_adj_stack,
                        sigma2_nonadj=sigma2_nonadj_stack)
    log_info("snj", f"Wrote {meta_path} + {npz_path.name} ({elapsed:.1f}s)", force=True)
    # Return same-shape dict for backward compatibility with notebook code.
    for idx, entry in enumerate(per_p):
        entry["sigma2_adj"] = sigma2_adj_stack[idx].tolist()
        entry["sigma2_nonadj"] = sigma2_nonadj_stack[idx].tolist()
    return meta
