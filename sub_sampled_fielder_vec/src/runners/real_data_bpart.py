"""Real-dataset B-matrix clan-partition sweep (600-tree benchmark).

For each (newick, fasta) pair:
  1. Parse FASTA → character matrix → JC distance D.
  2. Parse Newick → ground-truth top bipartition (root's two leaf sets).
  3. Apply uniform sub-sampling of D at fractions p ∈ cfg.experiment.p_values.
  4. Predict clan by sign of leading eigvec of B̂ = (I-11ᵀ/m) D̂ (I-11ᵀ/m).
  5. Measure NMI against ground-truth (step 2) and against full-D griffing (step 1).

Output: per-tree JSON rows aggregated into a summary JSON with mean±std over all
trees, binned by imbalance η.
"""
from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import dendropy
import scipy.spatial.distance

from ..core.similarity_builder import SimilarityMatrixBuilder
from ..utils.griffing import griffing_leading_eigvec
from ..utils.logging import log_info
from ..utils.summaries import save_json
from .nj_sweep import _impute_mean

_SIM_FLOOR = 1e-12


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def fasta_to_distance(fasta_path: str) -> Tuple[np.ndarray, List[str]]:
    """Read a FASTA alignment and return JC distance matrix + taxon names.

    Returns
    -------
    D : (n, n) float64 distance matrix (JC corrected, diagonal = 0)
    names : list of n taxon label strings
    """
    char_matrix = dendropy.DnaCharacterMatrix.get(
        path=fasta_path, schema="fasta"
    )
    names = [str(t.label) for t in char_matrix.taxon_namespace]
    # Build (n, L) character array as integers (A=0,C=1,G=2,T=3, gaps→4)
    _base_map = {c: i for i, c in enumerate("ACGT")}
    seqs = []
    for t in char_matrix.taxon_namespace:
        row = char_matrix[t]
        seq = np.array(
            [_base_map.get(str(s).upper(), 4) for s in row],
            dtype=np.int8,
        )
        seqs.append(seq)
    obs = np.vstack(seqs)  # (n, L)
    k = 4  # DNA
    hamming = scipy.spatial.distance.squareform(
        scipy.spatial.distance.pdist(obs.astype(np.float64), metric="hamming")
    )
    inside_log = np.clip(1.0 - hamming * k / (k - 1), _SIM_FLOOR, 1.0)
    R = inside_log ** (k - 1)
    np.fill_diagonal(R, 1.0)
    D = -np.log(R)
    np.fill_diagonal(D, 0.0)
    return D.astype(np.float64), names


def newick_top_bipartition(newick_path: str, names: List[str]) -> np.ndarray:
    """Return a boolean array (length n) representing the top bipartition.

    For a rooted tree: the root's two child subtrees define the bipartition.
    For an unrooted tree: we reroot at the first internal node.

    Parameters
    ----------
    names : taxon labels in the same order as rows/cols of D.

    Returns
    -------
    partition : bool array, True = first child's leaf set.
    """
    tree = dendropy.Tree.get(path=newick_path, schema="newick",
                             preserve_underscores=True)
    # Collect seed node's two child subtrees
    root = tree.seed_node
    children = list(root.child_nodes())
    if len(children) < 2:
        # Degenerate; reroot at first non-leaf
        for node in tree.preorder_node_iter():
            if len(list(node.child_nodes())) >= 2:
                root = node
                children = list(root.child_nodes())
                break

    def leaf_labels(node) -> set:
        return {leaf.taxon.label for leaf in node.leaf_iter()}

    left_labels = leaf_labels(children[0])
    name_set = {n: i for i, n in enumerate(names)}
    partition = np.array([n in left_labels for n in names], dtype=bool)
    # Degenerate: one side empty → use all-False (will give NMI=0 gracefully)
    n_left = int(partition.sum())
    if n_left == 0 or n_left == len(names):
        partition[:len(names) // 2] = True
    return partition


# ---------------------------------------------------------------------------
# Per-tree metrics
# ---------------------------------------------------------------------------

def _nmi(labels_true: np.ndarray, labels_pred: np.ndarray) -> float:
    from sklearn.metrics import normalized_mutual_info_score
    return float(normalized_mutual_info_score(labels_true.astype(int),
                                              labels_pred.astype(int)))


def _ari(labels_true: np.ndarray, labels_pred: np.ndarray) -> float:
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(labels_true.astype(int),
                                     labels_pred.astype(int)))


def _agreement(a: np.ndarray, b: np.ndarray) -> float:
    """Orientation-invariant % matching clan labels."""
    matches = int(np.sum(a == b))
    matches = max(matches, len(a) - matches)
    return 100.0 * matches / len(a)


def _eta_from_partition(partition: np.ndarray) -> float:
    """Imbalance ratio max/min of the two clan sizes."""
    n1 = int(partition.sum())
    n2 = len(partition) - n1
    return float(max(n1, n2)) / float(max(min(n1, n2), 1))


def _sweep_one_tree(
    D: np.ndarray,
    gt_partition: np.ndarray,
    p_values: List[float],
    bootstrap_reps: int,
    seed_base: int,
    imputation: str = "mean",
) -> Dict:
    """Run the sub-sampling sweep for one tree, return per-p metrics."""
    builder = SimilarityMatrixBuilder(method="uniform", matrix_kind="distance")
    sampler = builder.sampler

    # Reference: full-D griffing prediction
    v_ref = griffing_leading_eigvec(D)
    ref_partition = v_ref >= 0

    per_p = []
    for p_idx, p in enumerate(p_values):
        nmis_gt, aris_gt, agr_gt = [], [], []
        nmis_ref, aris_ref, agr_ref = [], [], []

        for rep in range(bootstrap_reps):
            seed = int(seed_base + 10_000 * p_idx + rep)
            D_hat = sampler.sample(D, p, seed=seed)
            if imputation == "mean":
                D_hat = _impute_mean(D_hat)

            v_hat = griffing_leading_eigvec(D_hat)
            pred_partition = v_hat >= 0

            # vs. ground-truth Newick bipartition
            nmis_gt.append(_nmi(gt_partition, pred_partition))
            aris_gt.append(_ari(gt_partition, pred_partition))
            agr_gt.append(_agreement(gt_partition, pred_partition))

            # vs. full-D reference
            nmis_ref.append(_nmi(ref_partition, pred_partition))
            aris_ref.append(_ari(ref_partition, pred_partition))
            agr_ref.append(_agreement(ref_partition, pred_partition))

        per_p.append({
            "p": float(p),
            "nmi_gt_mean": float(np.mean(nmis_gt)),
            "nmi_gt_std": float(np.std(nmis_gt)),
            "ari_gt_mean": float(np.mean(aris_gt)),
            "ari_gt_std": float(np.std(aris_gt)),
            "agr_gt_mean": float(np.mean(agr_gt)),
            "agr_gt_std": float(np.std(agr_gt)),
            "nmi_ref_mean": float(np.mean(nmis_ref)),
            "nmi_ref_std": float(np.std(nmis_ref)),
            "ari_ref_mean": float(np.mean(aris_ref)),
            "ari_ref_std": float(np.std(aris_ref)),
            "agr_ref_mean": float(np.mean(agr_ref)),
            "agr_ref_std": float(np.std(agr_ref)),
        })
    return {"per_p": per_p}


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

def run_real_data_benchmark(
    newick_dir: str,
    fasta_dir: str,
    run_dir: str,
    p_values: Optional[List[float]] = None,
    bootstrap_reps: int = 5,
    seed: int = 42,
    imputation: str = "mean",
    max_trees: Optional[int] = None,
    fasta_ext: str = ".fasta",
) -> Dict:
    """Run the 600-tree real-data benchmark.

    Iterates over all ``random_tree_*.nwk`` files in ``newick_dir``,
    loads the matching FASTA from ``fasta_dir``, runs the sub-sampling
    sweep, and writes two output files to ``run_dir``:

    - ``per_tree_results.json``  — one record per tree
    - ``aggregate_results.json`` — mean ± std NMI/ARI curves, η-binned

    Parameters
    ----------
    fasta_ext : extension of FASTA files (default ``.fasta``; also try ``.fa``)
    """
    if p_values is None:
        p_values = [round(0.1 * i, 1) for i in range(1, 11)]

    os.makedirs(run_dir, exist_ok=True)
    newick_paths = sorted(Path(newick_dir).glob("random_tree_*.nwk"))
    if max_trees is not None:
        newick_paths = newick_paths[:max_trees]

    log_info("real_data", f"Found {len(newick_paths)} trees in {newick_dir}", force=True)
    per_tree = []
    t0 = time.time()

    for i, nwk_path in enumerate(newick_paths):
        tree_id = nwk_path.stem  # e.g. "random_tree_42"
        # Match FASTA: try same stem + ext, or numbered variants
        fasta_path = Path(fasta_dir) / (tree_id + fasta_ext)
        if not fasta_path.exists():
            # try .fa
            fasta_path = Path(fasta_dir) / (tree_id + ".fa")
        if not fasta_path.exists():
            log_info("real_data", f"[{i+1}] MISSING fasta for {tree_id}, skipping", force=True)
            continue

        try:
            D, names = fasta_to_distance(str(fasta_path))
            gt_partition = newick_top_bipartition(str(nwk_path), names)
            eta = _eta_from_partition(gt_partition)
            n_taxa = len(names)

            result = _sweep_one_tree(
                D, gt_partition, p_values,
                bootstrap_reps=bootstrap_reps,
                seed_base=seed + i * 100_000,
                imputation=imputation,
            )
            result.update({
                "tree_id": tree_id,
                "n_taxa": n_taxa,
                "eta": eta,
                "gt_sizes": [int(gt_partition.sum()),
                             int((~gt_partition).sum())],
            })
            per_tree.append(result)

        except Exception as exc:
            log_info("real_data", f"[{i+1}] ERROR on {tree_id}: {exc}", force=True)
            continue

        if (i + 1) % 50 == 0 or i == 0:
            elapsed = time.time() - t0
            log_info("real_data",
                     f"  {i+1}/{len(newick_paths)} trees done ({elapsed:.0f}s)",
                     force=True)

    save_json({"trees": per_tree}, str(Path(run_dir) / "per_tree_results.json"))
    agg = _aggregate(per_tree, p_values)
    save_json(agg, str(Path(run_dir) / "aggregate_results.json"))
    log_info("real_data",
             f"Done. {len(per_tree)} trees processed in {time.time()-t0:.1f}s",
             force=True)
    return agg


def _aggregate(per_tree: List[Dict], p_values: List[float]) -> Dict:
    """Compute mean±std NMI/ARI over all trees per p, plus η-binned curves."""
    n = len(per_tree)
    if n == 0:
        return {"n_trees": 0, "p_values": p_values}

    metrics = ["nmi_gt", "ari_gt", "agr_gt", "nmi_ref", "ari_ref", "agr_ref"]
    global_per_p = []
    for pi, p in enumerate(p_values):
        row: Dict = {"p": float(p)}
        for m in metrics:
            vals = [t["per_p"][pi][f"{m}_mean"] for t in per_tree]
            row[f"{m}_mean"] = float(np.mean(vals))
            row[f"{m}_std"] = float(np.std(vals))
            row[f"{m}_median"] = float(np.median(vals))
        global_per_p.append(row)

    # η-stratified: bin by log2(η) into ~3 bands
    eta_values = np.array([t["eta"] for t in per_tree])
    # Bins: balanced [1, 2), moderate [2, 5), imbalanced [5, ∞)
    eta_bins = [(1.0, 2.0, "balanced"), (2.0, 5.0, "moderate"),
                (5.0, np.inf, "imbalanced")]
    binned = []
    for lo, hi, label in eta_bins:
        mask = (eta_values >= lo) & (eta_values < hi)
        subset = [t for t, m in zip(per_tree, mask) if m]
        if not subset:
            continue
        bin_per_p = []
        for pi, p in enumerate(p_values):
            brow: Dict = {"p": float(p)}
            for met in metrics:
                vals = [t["per_p"][pi][f"{met}_mean"] for t in subset]
                brow[f"{met}_mean"] = float(np.mean(vals))
                brow[f"{met}_std"] = float(np.std(vals))
            bin_per_p.append(brow)
        binned.append({
            "label": label,
            "eta_range": [lo, hi if np.isfinite(hi) else None],
            "n_trees": len(subset),
            "per_p": bin_per_p,
        })

    return {
        "n_trees": n,
        "p_values": [float(p) for p in p_values],
        "global_per_p": global_per_p,
        "eta_binned": binned,
    }
