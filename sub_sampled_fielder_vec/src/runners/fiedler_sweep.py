"""First-layer Fiedler-recovery sweep over ``p`` for a fixed (n_taxa, seq_len).

Similarity-matrix analogue of ``griffing_sweep.griffing_sweep_for_params``:

  1. Build the JC similarity matrix ``S`` from sequences.
  2. Reference partition vector: Fiedler of ``L(S) = D_S - S`` (full data).
  3. Per (p, rep): sub-sample ``S → Ŝ`` (uniform IPW, diagonal preserved at 1.0),
     compute Fiedler of ``L(Ŝ)``, record sign-agnostic agreement.
  4. Aggregate.

This is the "first layer" of STDR's recursive partition phase — one binary
cut of the n taxa. Recursive (final) partition is intentionally out of
scope here; see ``stdr_partition_recovery.ipynb`` for the recursive
machinery.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import spectraltree

from ..config import StructuredConfig
from ..core.similarity_builder import SimilarityMatrixBuilder
from ..models import get_tree_factory, get_sequence_factory
from ..utils.logging import log_info
from ..utils.metrics import estimate_operator_norm_diff, compute_sign_agreement
from ..utils.summaries import save_json
from ..core.utils import compute_laplacian
from ..core.fiedler_computer import FiedlerVectorComputer


def compute_fiedler_of_S(S: np.ndarray, sampling_prob: float | None = None) -> np.ndarray:
    """Fiedler vector of ``L = D - S`` with the project's sign convention.

    Inlined from ``analysis.theoretical_interpretation.utils.spectral``, which this
    module used to reach via a ``sys.path`` hack -- a src -> analysis dependency
    inversion. Both pieces it needs already live in ``src/core``, so the import is
    unnecessary as well as backwards.
    """
    return FiedlerVectorComputer().compute(
        compute_laplacian(S), sampling_prob=sampling_prob
    )


def _build_truth_similarity(cfg: StructuredConfig, n_taxa: int, seq_len: int
                            ) -> Tuple[object, np.ndarray, np.ndarray, object]:
    """Generate tree + observations + full JC similarity matrix S + TaxaMetadata.

    Mirrors ``nj_sweep._build_truth`` but returns S instead of D (skips the
    -log step). Diagonal of S is enforced to 1.0.
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
    S = spectraltree.JC_similarity_matrix(char_matrix)
    S = np.clip(S, 0.0, 1.0)
    np.fill_diagonal(S, 1.0)
    return tree, char_matrix, S, taxa_metadata


def fiedler_first_layer_sweep_for_params(cfg: StructuredConfig, n_taxa: int,
                                         seq_len: int, run_dir: str) -> Dict:
    """Run the first-layer Fiedler-recovery sweep at fixed (n_taxa, seq_len).

    Writes ``fiedler_meta.json`` into ``run_dir`` and returns the in-memory
    dict. Uses uniform IPW sampling on the similarity matrix; the diagonal
    is preserved at 1.0 (self-similarity is always observed).
    """
    os.makedirs(run_dir, exist_ok=True)
    log_info("fiedler", f"Building ground truth (n={n_taxa}, L={seq_len})...", force=True)
    tree, observations, S, taxa_metadata = _build_truth_similarity(cfg, n_taxa, seq_len)

    log_info("fiedler", "Computing reference v_ref = Fiedler(L(S))...", force=True)
    v_ref = compute_fiedler_of_S(S)

    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="similarity")
    sampler = builder.sampler  # self_value=1.0 for similarity diagonals

    p_values: List[float] = list(cfg.experiment.p_values)
    bootstrap_reps: int = int(cfg.experiment.bootstrap_reps)
    seed_base: int = int(cfg.experiment.seed)

    per_p: List[Dict] = []
    t0 = time.time()
    for p_idx, p in enumerate(p_values):
        log_info("fiedler", f"  p = {p:.4g} ({bootstrap_reps} reps)", force=True)
        spec_norms: List[float] = []
        recoveries: List[float] = []
        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            S_hat = sampler.sample(S, p, seed=seed)
            spec_norm = float(estimate_operator_norm_diff(S_hat, S, n_iter=20))
            spec_norms.append(spec_norm)
            try:
                v_hat = compute_fiedler_of_S(S_hat, sampling_prob=float(p))
                s = compute_sign_agreement(v_ref, v_hat) / 100.0
                recovery = float(max(s, 1.0 - s))
            except Exception as e:
                log_info("fiedler", f"    rep {rep}: Fiedler failed ({e!r}); recovery=0.5",
                         force=True)
                recovery = 0.5
            recoveries.append(recovery)

        per_p.append({
            "p": float(p),
            "spec_norm_mean": float(np.mean(spec_norms)),
            "spec_norm_std": float(np.std(spec_norms)),
            "spec_norm_per_rep": [float(x) for x in spec_norms],
            "recovery_mean": float(np.mean(recoveries)),
            "recovery_std": float(np.std(recoveries)),
            "recovery_per_rep": [float(x) for x in recoveries],
        })

    elapsed = time.time() - t0
    meta = {
        "n_taxa": n_taxa,
        "seq_len": seq_len,
        "p_values": [float(p) for p in p_values],
        "bootstrap_reps": bootstrap_reps,
        "tree_model": cfg.tree.model,
        "seq_model": cfg.sequence.model,
        "mutation_rate": float(cfg.sequence.params.get("mutation_rate", 1.0)),
        "tree_params": dict(cfg.tree.params),
        "sampling_method": "uniform_similarity_ipw",
        "method": "fiedler_first_layer",
        "per_p": per_p,
        "elapsed_seconds": elapsed,
    }
    run_path = Path(run_dir)
    meta_path = run_path / "fiedler_meta.json"
    save_json(meta, str(meta_path))
    log_info("fiedler", f"Wrote {meta_path} ({elapsed:.1f}s)", force=True)
    return meta
