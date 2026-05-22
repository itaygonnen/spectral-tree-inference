"""Classical-NJ-with-subsampling sweep over ``p`` for a fixed (n_taxa, seq_len).

Mirrors ``snj_sweep.snj_sweep_for_params`` but runs the *classical* Neighbor
Joining algorithm on a sub-sampled JC distance matrix:

  1. ``||D - D̂||_2``  (spec-norm error, NJ analog of SNJ Thm 4.2)
  2. Q-criterion separation at step 0 (Fig. 3 analog, adj vs. non-adj)
  3. Robinson-Foulds via NJ on the sub-sampled distance matrix

Sampling is plain uniform; ``imputation`` controls what is fed to NJ:

  * ``"zero"`` (IPW unbiased estimator): unsampled entries are 0,
    sampled entries are scaled by ``1/p``. Faithful to the SNJ
    no-completion baseline, but catastrophic for NJ — zeros look like
    "identical taxa" and dominate ``argmin Q``, producing spurious
    cherries everywhere even at high p. RF saturates to max for any
    p<1.

  * ``"mean"`` (default): same IPW sampling, then unsampled entries are
    overwritten with the mean of the sampled+scaled entries. Curves
    degrade gracefully; the spectral-norm metric is unchanged in spirit
    (sampling is identical).
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
from spectraltree import utils as st_utils

from ..config import StructuredConfig
from ..core.similarity_builder import SimilarityMatrixBuilder
from ..models import get_tree_factory, get_sequence_factory
from ..utils.logging import log_info
from ..utils.metrics import estimate_operator_norm_diff
from ..utils.nj_q_separation import compute_q_separation
from ..utils.summaries import save_json


_SIM_FLOOR = 1e-12  # protects -log(0) when computing D = -log(R)


def _impute_mean(D_hat: np.ndarray) -> np.ndarray:
    """Overwrite exactly-zero off-diagonal entries with the mean of the
    strictly-positive (i.e. sampled, IPW-scaled) entries. Diagonal stays 0.

    This is a min-change fix to the zero-fill catastrophe described in the
    module docstring: without it, NJ's first merge is always a phantom
    cherry on an unsampled pair (D̂ = 0 → very negative Q), producing
    maximum RF for any p < 1.
    """
    D = D_hat.copy()
    off_mask = ~np.eye(D.shape[0], dtype=bool)
    pos = (D > 0) & off_mask
    zero_off = (D == 0) & off_mask
    if zero_off.any() and pos.any():
        mean_val = float(D[pos].mean())
        D[zero_off] = mean_val
    return D


def _build_truth(cfg: StructuredConfig, n_taxa: int, seq_len: int
                 ) -> Tuple[object, np.ndarray, np.ndarray, object]:
    """Generate tree + observations + full JC distance matrix D + TaxaMetadata.

    D = -log(JC_similarity_matrix(obs)) — matches the LaTeX algorithm's
    Step 2 (JC distance formula). The (d-1) prefactor in the paper is a
    constant, so we drop it; it scales D uniformly and does not affect Q
    ordering or RF.
    """
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
    R = np.clip(R, _SIM_FLOOR, 1.0)
    np.fill_diagonal(R, 1.0)
    D = -np.log(R)
    np.fill_diagonal(D, 0.0)
    return tree, char_matrix, D, taxa_metadata


def _run_nj(D_hat: np.ndarray, taxa_metadata):
    """Run classical NJ on a precomputed distance matrix.

    Bypasses ``spectraltree.NeighborJoining.reconstruct_from_similarity``
    (which would apply an extra -log) and feeds D_hat straight into the
    dendropy NJ path used at ``spectraltree/reconstruct_tree.py:174``.
    """
    dm = st_utils.array2distance_matrix(D_hat, taxa_metadata)
    return dm.nj_tree()


def _rf_distance(ref_tree, inferred_tree, taxon_namespace) -> int:
    """Robinson-Foulds (symmetric difference). Both trees must share
    ``taxon_namespace``. Bipartitions are encoded in-place.
    """
    for tr in (ref_tree, inferred_tree):
        tr.encode_bipartitions()
    return int(treecompare.symmetric_difference(ref_tree, inferred_tree))


def nj_sweep_for_params(cfg: StructuredConfig, n_taxa: int, seq_len: int,
                        run_dir: str, *, imputation: str = "mean") -> Dict:
    """Run the NJ-subsampling sweep at fixed (n_taxa, seq_len).

    Writes ``nj_meta.json`` + ``nj_arrays.npz`` into ``run_dir`` and
    returns the in-memory dict (with Q arrays re-attached) for
    notebook-style use.

    Parameters
    ----------
    imputation : {"mean", "zero"}
        How to handle unsampled entries (see module docstring). Default
        ``"mean"`` produces informative curves; ``"zero"`` reproduces the
        IPW-unbiased estimator and shows the failure mode.
    """
    if imputation not in ("mean", "zero"):
        raise ValueError(f"imputation must be 'mean' or 'zero', got {imputation!r}")
    os.makedirs(run_dir, exist_ok=True)
    log_info("nj", f"Building ground truth (n={n_taxa}, L={seq_len})...", force=True)
    tree, observations, D, taxa_metadata = _build_truth(cfg, n_taxa, seq_len)
    tree.encode_bipartitions()

    # Use only the sampler from the builder so we bypass build_subsampled's
    # exp(-α·D̂) post-processing (we want raw D̂, not similarity-shaped).
    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="distance")
    sampler = builder.sampler  # UniformSampler with self_value=0.0

    p_values: List[float] = list(cfg.experiment.p_values)
    bootstrap_reps: int = int(cfg.experiment.bootstrap_reps)
    seed_base: int = int(cfg.experiment.seed)

    per_p: List[Dict] = []
    t0 = time.time()
    for p_idx, p in enumerate(p_values):
        log_info("nj", f"  p = {p:.4g} ({bootstrap_reps} reps)", force=True)
        spec_norms: List[float] = []
        rf_values: List[int] = []
        adj_all: List[np.ndarray] = []
        nonadj_all: List[np.ndarray] = []

        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            D_hat = sampler.sample(D, p, seed=seed)
            if imputation == "mean":
                D_hat = _impute_mean(D_hat)

            # Metric 1: spectral-norm error on D
            spec_norm = float(estimate_operator_norm_diff(D_hat, D, n_iter=20))
            spec_norms.append(spec_norm)

            # Metric 2: Q-criterion separation at step 0
            adj, nonadj = compute_q_separation(D_hat, tree, taxa_metadata)
            adj_all.append(adj)
            nonadj_all.append(nonadj)

            # Metric 3: RF distance via classical NJ
            try:
                inferred = _run_nj(D_hat, taxa_metadata)
                rf = _rf_distance(
                    ref_tree=dendropy.Tree(tree),
                    inferred_tree=inferred,
                    taxon_namespace=taxa_metadata.taxon_namespace,
                )
            except Exception as e:
                log_info("nj", f"    rep {rep}: NJ failed ({e!r}); RF=2(m-3)", force=True)
                rf = 2 * (n_taxa - 3)
            rf_values.append(int(rf))

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
            "_q_adj": adj_concat,
            "_q_nonadj": nonadj_concat,
        })

    elapsed = time.time() - t0
    q_adj_stack = np.stack([entry.pop("_q_adj") for entry in per_p])
    q_nonadj_stack = np.stack([entry.pop("_q_nonadj") for entry in per_p])

    meta = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": [float(p) for p in p_values],
        "bootstrap_reps": bootstrap_reps,
        "tree_model": cfg.tree.model,
        "seq_model": cfg.sequence.model,
        "mutation_rate": float(cfg.sequence.params.get("mutation_rate", 1.0)),
        "sampling_method": f"uniform_imputation={imputation}",
        "method": "classical_nj",
        "distance_kind": "jc_minus_log",
        "imputation": imputation,
        "per_p": per_p,
        "elapsed_seconds": elapsed,
    }
    run_path = Path(run_dir)
    meta_path = run_path / "nj_meta.json"
    npz_path = run_path / "nj_arrays.npz"
    save_json(meta, str(meta_path))
    np.savez_compressed(npz_path, q_adj=q_adj_stack, q_nonadj=q_nonadj_stack)
    log_info("nj", f"Wrote {meta_path} + {npz_path.name} ({elapsed:.1f}s)", force=True)
    # Re-attach arrays for notebook callers.
    for idx, entry in enumerate(per_p):
        entry["q_adj"] = q_adj_stack[idx].tolist()
        entry["q_nonadj"] = q_nonadj_stack[idx].tolist()
    return meta
