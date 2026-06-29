"""Shared screen / analysis / plot logic for the distance-vs-similarity notebooks.

Two notebooks consume this module — `simulation_distance_vs_similarity.ipynb`
(in-memory generator, cache files prefixed ``gen_``) and
`real_data_distance_vs_similarity.ipynb` (FASTA+Newick loader, un-prefixed
caches). They differ only in (1) the data ``loader`` and (2) the cache prefix,
so every function here takes the loader / cohort / cache_prefix as arguments and
is otherwise data-source agnostic.

Loader contract: ``loader(tid) -> (S, labels, tree, D)`` with ``D`` aligned to
``labels`` order, or ``None`` to skip a missing tree.

Cache filenames are preserved exactly. The base names below combine with
``cache_prefix`` ("gen_" for sim, "" for real) to reproduce each notebook's
current filenames:
  screen 1 (Fiedler-S):   gen_screen_cache.npz   / screen_cache_600.npz
  screen 2 (Griffing-D):  gen_screen_griffing.npz/ screen_griffing_600.npz
  screen 3 (L_sym):       gen_screen_lsym.npz    / screen_lsym_600.npz
The real-data screens use a "_600" suffix, so callers pass the full base name
explicitly via ``screen_basenames`` (see :func:`run_operator_screens`).

All plotting fns keep the exact colors / thresholds / labels of the originals;
pure-compute fns return arrays/dicts. ``cache_dir`` defaults to ``Path('.')`` so
the thin-wrapper notebooks keep reading the .npz caches sitting next to them.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

import scipy.linalg
from sklearn.metrics import normalized_mutual_info_score

from src.core.utils import (
    compute_laplacian,
    compute_fiedler_from_laplacian,
    compute_normalized_laplacian,
)
from src.utils.metrics import compute_reference_partition_and_quality
from src.utils.partition_validity import check_partition_valid_in_tree
from src.runners.p_sweep_inner import (
    bootstrap_p_sweep_simple,
    _kmeans_bipartition,
    _subsample_matrix_entries,
)
from src.utils.screening import run_screen
from src.utils.griffing import griffing_leading_eigvec, griffing_centered
from src.utils.bpart_sweep_cache import bpart_sweep_raw
from spectraltree.spectral_tree_reconstruction import partition_taxa

# --- shared semantics (sign-partition colors, method palette) ---
COLOR_POS, COLOR_NEG = "#d97706", "#1d4ed8"          # orange / blue (sign partition)
COL_S, COL_D, COL_LS = "#1d4ed8", "#065f46", "#b45309"  # method colors: S / D / L_sym

# default screen base names (un-prefixed); real notebook overrides with "_600" variants
SCREEN_BASENAMES = {
    "fiedler": "screen_cache.npz",
    "griffing": "screen_griffing.npz",
    "lsym": "screen_lsym.npz",
}


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def tree_leaf_index(tree, labels) -> np.ndarray:
    """Permutation mapping tree leaf order -> ``labels`` order."""
    idx = {l: i for i, l in enumerate(labels)}
    return np.array([idx[str(l.taxon.label)] for l in tree.leaf_node_iter()])


def _cache(cache_dir, cache_prefix: str, basename: str) -> Path:
    return Path(cache_dir) / f"{cache_prefix}{basename}"


def patristic_distance_matrix(tree) -> np.ndarray:
    leaf_taxa = [leaf.taxon for leaf in tree.leaf_node_iter()]
    n = len(leaf_taxa)
    pdm = tree.phylogenetic_distance_matrix()
    Dp = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = pdm.distance(leaf_taxa[i], leaf_taxa[j])
            Dp[i, j] = Dp[j, i] = d
    return Dp


def sign_orient(v) -> np.ndarray:
    """Sign partition, with the majority side forced to 1 (orange) for consistency."""
    s = (np.asarray(v) >= 0).astype(int)
    if s.sum() < len(s) - s.sum():
        s = 1 - s
    return s


def _fast_fiedler_Lsym(M: np.ndarray) -> np.ndarray:
    deg = np.asarray(M.sum(axis=1)).ravel()
    dis = 1.0 / np.sqrt(np.maximum(deg, 1e-12))
    Lsym = np.eye(len(deg)) - (M * dis[:, None]) * dis[None, :]
    Lsym = 0.5 * (Lsym + Lsym.T)
    _, vecs = scipy.linalg.eigh(Lsym, subset_by_index=[0, 1])
    return vecs[:, 1]


def p_star_per_tree(nmi_arr, p_values, thresh: float = 0.95) -> np.ndarray:
    """Smallest p with NMI >= thresh per tree (inf if never reached)."""
    out = []
    for row in np.asarray(nmi_arr, float):
        ok = np.where(row >= thresh)[0]
        out.append(float(p_values[ok[0]]) if ok.size else np.inf)
    return np.array(out)


# ----------------------------------------------------------------------------
# operator screens (cells 3/5/6)
# ----------------------------------------------------------------------------
def run_operator_screens(
    loader: Callable,
    tree_ids: Sequence[str],
    *,
    num_gaps: int,
    min_split: int,
    cache_dir=Path("."),
    cache_prefix: str = "",
    screen_basenames: Optional[Dict[str, str]] = None,
) -> dict:
    """Run the three operator screens (Fiedler-on-S, Griffing-on-D, L_sym).

    Returns a dict with keys ``fiedler`` / ``griffing`` / ``lsym``, each a
    ``(valid_ids, etas, rows)`` tuple from :func:`run_screen`, plus the rows
    convenience aliases ``rows`` / ``griff_rows`` / ``lsym_rows``.
    """
    bn = dict(SCREEN_BASENAMES)
    if screen_basenames:
        bn.update(screen_basenames)

    fiedler = run_screen(
        tree_ids, loader,
        lambda S, D: compute_reference_partition_and_quality(
            compute_fiedler_from_laplacian(compute_laplacian(S)), S,
            num_gaps=num_gaps, min_split=min_split)[0],
        _cache(cache_dir, cache_prefix, bn["fiedler"]), label="Fiedler-on-S")

    griffing = run_screen(
        tree_ids, loader,
        lambda S, D: (griffing_leading_eigvec(D) >= 0),
        _cache(cache_dir, cache_prefix, bn["griffing"]), label="Griffing-on-D")

    lsym = run_screen(
        tree_ids, loader,
        lambda S, D: compute_reference_partition_and_quality(
            compute_fiedler_from_laplacian(compute_normalized_laplacian(S)), S,
            num_gaps=num_gaps, min_split=min_split)[0],
        _cache(cache_dir, cache_prefix, bn["lsym"]), label="L_sym")

    return {
        "fiedler": fiedler, "griffing": griffing, "lsym": lsym,
        "rows": fiedler[2], "griff_rows": griffing[2], "lsym_rows": lsym[2],
    }


# ----------------------------------------------------------------------------
# eta histograms (cell 7)
# ----------------------------------------------------------------------------
def plot_eta_histograms(etas_valid, griff_etas, lsym_etas):
    """Side-by-side eta histograms for the three operators, each on its own x-range."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    panels = [(axes[0], etas_valid, r"Fiedler of $L=\mathrm{Deg}(S)-S$"),
              (axes[1], griff_etas, r"Griffing of $B=HDH$"),
              (axes[2], lsym_etas,  r"Fiedler of $L_{\rm sym}=I-D^{-1/2}SD^{-1/2}$")]
    for ax, data, name in panels:
        ax.hist(data, bins=30, range=(1, max(data.max() * 1.05, 2.0)),
                color="#9ca3af", edgecolor="white")
        ax.set_xlabel(r"$\eta$ = larger clan / smaller clan", fontsize=12)
        ax.set_title(f"{name}\n{len(data)} valid trees, median $\\eta$={np.median(data):.1f}",
                     fontsize=12)
    axes[0].set_ylabel("# valid trees", fontsize=12)
    fig.suptitle("Reference-partition imbalance decides recoverability", fontsize=13)
    fig.tight_layout()
    plt.show()


# ----------------------------------------------------------------------------
# operator x threshold sweep (cell 9) + grid (cell 10)
# ----------------------------------------------------------------------------
def _operators():
    return {
        "L":     lambda S, D: compute_fiedler_from_laplacian(compute_laplacian(S)),
        "L_sym": lambda S, D: compute_fiedler_from_laplacian(compute_normalized_laplacian(S)),
        "B":     lambda S, D: griffing_leading_eigvec(D),
    }


def _thresholds(num_gaps: int, min_split: int):
    return {
        "$\\sigma_2$": lambda v, S: partition_taxa(v, S, num_gaps, min_split),
        "k-means":     lambda v, S: _kmeans_bipartition(v, min_split),
        "sign":        lambda v, S: (v >= 0),
    }


def run_threshold_sweep(
    loader: Callable,
    tree_ids: Sequence[str],
    *,
    num_gaps: int,
    min_split: int,
    meta: dict,
    cache_basename: str,
    cache_dir=Path("."),
    cache_prefix: str = "",
    progress: int = 50,
) -> dict:
    """Valid-only operator x threshold partition-imbalance sweep.

    ``meta`` is the cache validity key (sim: ``dict(ids, m, L, min_split,
    num_gaps)``; real: ``dict(ids, min_split, num_gaps)``) — passed through so
    each notebook reproduces its own cache. Trees the loader returns ``None``
    for are skipped (real fasta gaps); ``n_screened`` counts those actually used.

    Returns ``{thr_etas, n_screened, op_names, thr_names}`` where ``thr_etas``
    maps ``(op, thr) -> eta array``.
    """
    operators = _operators()
    thresholds = _thresholds(num_gaps, min_split)
    op_names, thr_names = list(operators), list(thresholds)

    cache_path = _cache(cache_dir, cache_prefix, cache_basename)
    thr_etas = None
    n_screened = 0
    if cache_path.exists():
        z = np.load(cache_path, allow_pickle=True)
        if dict(z["meta"].item()) == meta:
            thr_etas = {tuple(k.split("|")): np.asarray(v, float)
                        for k, v in z["etas"].item().items()}
            n_screened = int(z["n_screened"]) if "n_screened" in z else len(tree_ids)
            print(f"loaded threshold-eta cache ({cache_path})")

    if thr_etas is None:
        acc = {(op, thr): [] for op in op_names for thr in thr_names}
        n_screened = 0
        for k, tid in enumerate(tree_ids, 1):
            loaded = loader(tid)
            if loaded is None:
                continue
            S, labels, tree, D = loaded
            n_screened += 1
            sidx = tree_leaf_index(tree, labels)
            for op in op_names:
                v = operators[op](S, D)
                for thr in thr_names:
                    try:
                        part = np.asarray(thresholds[thr](v, S)).astype(bool)
                    except Exception:
                        continue
                    if not bool(check_partition_valid_in_tree(tree, part[sidx])):
                        continue
                    n1 = int(part.sum()); n2 = len(part) - n1
                    acc[(op, thr)].append(max(n1, n2) / max(min(n1, n2), 1))
            if k % progress == 0:
                print(f"  screened {k}/{len(tree_ids)}")
        thr_etas = {key: np.asarray(val, float) for key, val in acc.items()}
        np.savez(cache_path, meta=meta, n_screened=n_screened,
                 etas={f"{op}|{thr}": thr_etas[(op, thr)]
                       for op in op_names for thr in thr_names})
        print(f"saved threshold-eta cache ({cache_path})")

    print(f"\nvalid partitions / {n_screened} screened:")
    print(f"  {'op':>6} | " + " | ".join(f"{t:>8}" for t in thr_names))
    for op in op_names:
        print(f"  {op:>6} | " +
              " | ".join(f"{len(thr_etas[(op, t)]):>8d}" for t in thr_names))

    return {"thr_etas": thr_etas, "n_screened": n_screened,
            "op_names": op_names, "thr_names": thr_names}


def plot_operator_threshold_grid(sweep: dict, *, dataset_label: str):
    """3x3 histogram grid: rows = operator, cols = threshold, own x-range per panel."""
    thr_etas = sweep["thr_etas"]; n = sweep["n_screened"]
    op_names, thr_names = sweep["op_names"], sweep["thr_names"]
    fig, axes = plt.subplots(len(op_names), len(thr_names), figsize=(15, 12))
    for r, op in enumerate(op_names):
        for c, thr in enumerate(thr_names):
            ax = axes[r, c]
            data = np.asarray(thr_etas[(op, thr)], float)
            if data.size:
                ax.hist(data, bins=30, range=(1, max(data.max() * 1.05, 2.0)),
                        color="#9ca3af", edgecolor="white")
                med = np.median(data)
            else:
                med = float("nan")
            ax.set_title(f"{op} × {thr}\n{data.size}/{n} valid, median $\\eta$={med:.1f}",
                         fontsize=11)
            if r == len(op_names) - 1:
                ax.set_xlabel(r"$\eta$ = larger clan / smaller clan", fontsize=11)
            if c == 0:
                ax.set_ylabel(f"{op}\n# valid trees", fontsize=12)
    fig.suptitle(f"Partition imbalance by operator × threshold rule — {dataset_label}",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    plt.show()


# ----------------------------------------------------------------------------
# ranked bars (cell 12)
# ----------------------------------------------------------------------------
def _median_ci(x, rng, n_boot=2000, alpha=0.05):
    x = np.asarray(x, float)
    if x.size == 0:
        return np.nan, np.nan, np.nan
    med = float(np.median(x))
    if x.size == 1:
        return med, med, med
    boots = np.median(rng.choice(x, size=(n_boot, x.size), replace=True), axis=1)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return med, float(lo), float(hi)


def plot_ranked_bars(sweep: dict, *, dataset_label: str):
    """Two barplots: (left) median eta + 95% bootstrap CI; (right) # valid partitions."""
    thr_etas = sweep["thr_etas"]; n_total = sweep["n_screened"]
    op_names, thr_names = sweep["op_names"], sweep["thr_names"]
    rng = np.random.default_rng(0)

    thr_color = {thr_names[0]: "#1d4ed8", thr_names[1]: "#065f46", thr_names[2]: "#b45309"}
    recs = []
    for op in op_names:
        for thr in thr_names:
            x = thr_etas[(op, thr)]
            med, lo, hi = _median_ci(x, rng)
            recs.append({"lab": f"{op} × {thr}", "thr": thr,
                         "med": med, "lo": lo, "hi": hi, "k": len(x)})

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(18, 6))

    sL = sorted([r for r in recs if np.isfinite(r["med"])],
                key=lambda r: r["med"], reverse=True)
    xL = np.arange(len(sL))
    medsL = np.array([r["med"] for r in sL])
    losL = np.array([r["lo"] for r in sL]); hisL = np.array([r["hi"] for r in sL])
    axL.bar(xL, medsL, yerr=np.vstack([medsL - losL, hisL - medsL]), capsize=4,
            color=[thr_color[r["thr"]] for r in sL], edgecolor="#374151")
    _padL = 0.05 * float(hisL.max())
    for xi, r in zip(xL, sL):
        axL.text(xi, r["hi"] + _padL, f"{r['med']:.1f}", ha="center", va="bottom",
                 fontsize=14, fontweight="bold")
    axL.set_ylim(0, float(hisL.max()) * 1.28)
    axL.set_xticks(xL); axL.set_xticklabels([r["lab"] for r in sL], rotation=30, ha="right")
    axL.set_ylabel(r"median partition imbalance $\eta$")
    axL.set_title("Median $\\eta$ (95% bootstrap CI)", pad=12)
    axL.grid(True, axis="y", alpha=0.3)

    sR = sorted(recs, key=lambda r: r["k"], reverse=True)
    xR = np.arange(len(sR))
    ksR = np.array([r["k"] for r in sR])
    axR.bar(xR, ksR, color=[thr_color[r["thr"]] for r in sR], edgecolor="#374151")
    _padR = 0.02 * n_total
    for xi, r in zip(xR, sR):
        axR.text(xi, r["k"] + _padR, f"{r['k']}\n({100 * r['k'] / n_total:.0f}%)",
                 ha="center", va="bottom", fontsize=14, fontweight="bold")
    axR.set_ylim(0, n_total * 1.25)
    axR.set_xticks(xR); axR.set_xticklabels([r["lab"] for r in sR], rotation=30, ha="right")
    axR.set_ylabel(f"# valid partitions (of {n_total})")
    axR.set_title("Valid partitions", pad=12)
    axR.grid(True, axis="y", alpha=0.3)

    fig.legend(handles=[Patch(color=thr_color[t], label=t) for t in thr_names],
               title="threshold", loc="upper center", bbox_to_anchor=(0.5, 0.92),
               ncol=3, frameon=False)
    fig.suptitle(f"Operator × threshold on {dataset_label}", fontsize=14, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    plt.show()


# ----------------------------------------------------------------------------
# validity agreement matrix (cell 14)
# ----------------------------------------------------------------------------
def plot_validity_agreement_matrix(screens: dict, tree_ids: Sequence[str]) -> dict:
    """3x3 validity-overlap heatmap across the three methods.

    Returns the downstream cohorts ``{both_valid_ids, only_S_ids, only_D_ids,
    neither_count}`` (S vs D) consumed by the recovery / demo panels.
    """
    fiedler_valid_ids = list(screens["fiedler"][0])
    griff_valid_ids = list(screens["griffing"][0])
    lsym_valid_ids = list(screens["lsym"][0])

    S_set = set(fiedler_valid_ids); D_set = set(griff_valid_ids); L_set = set(lsym_valid_ids)
    N_TOTAL = len(tree_ids)

    both_valid_ids = sorted(S_set & D_set)
    only_S_ids = sorted(S_set - D_set)
    only_D_ids = sorted(D_set - S_set)
    neither_count = N_TOTAL - len(S_set | D_set)

    methods = [("Fiedler-$S$", S_set), ("Griffing-$D$", D_set), ("$L_{\\rm sym}$", L_set)]
    names = [m for m, _ in methods]
    M = np.array([[len(si & sj) for _, sj in methods] for _, si in methods])
    print("pairwise both-valid overlap (diagonal = total valid per method):")
    print(pd.DataFrame(M, index=names, columns=names))
    print(f"\nvalid in all THREE methods: {len(S_set & D_set & L_set)}")
    print(f"valid in at least one:      {len(S_set | D_set | L_set)} / {N_TOTAL}")

    fig, ax = plt.subplots(figsize=(5.5, 5))
    im = ax.imshow(M, cmap="Blues")
    ax.set_xticks(range(3)); ax.set_xticklabels(names)
    ax.set_yticks(range(3)); ax.set_yticklabels(names)
    _thr = M.max() * 0.5 if M.max() else 1
    for i in range(3):
        for j in range(3):
            ax.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=14,
                    color="white" if M[i, j] > _thr else "#111827")
    ax.set_title("Validity agreement matrix\n(both-valid overlap; diagonal = total valid)",
                 fontsize=12)
    fig.colorbar(im, fraction=0.046, pad=0.04); fig.tight_layout(); plt.show()

    print(f"\ncohort used by every downstream P*/recovery panel: "
          f"{len(both_valid_ids)} trees valid in both Fiedler-$S$ and Griffing-$D$")
    print(f"  only-S: {len(only_S_ids)}   only-D: {len(only_D_ids)}   neither: {neither_count}")

    return {"both_valid_ids": both_valid_ids, "only_S_ids": only_S_ids,
            "only_D_ids": only_D_ids, "neither_count": neither_count}


# ----------------------------------------------------------------------------
# one tree, three partitions (cells 16/17)
# ----------------------------------------------------------------------------
def pick_demo_tree(screens: dict):
    """Tree valid in both Fiedler-S and Griffing-D, nearest the griffing-median eta."""
    fiedler_valid_ids = list(screens["fiedler"][0])
    griff_valid_ids = list(screens["griffing"][0])
    griff_etas = np.asarray(screens["griffing"][1], float)
    _fset = set(fiedler_valid_ids)
    both = [t for t in griff_valid_ids if t in _fset]
    _gmap = {t: e for t, e in zip(griff_valid_ids, griff_etas)}
    _med = np.median([_gmap[t] for t in both])
    tid = min(both, key=lambda t: abs(_gmap[t] - _med))
    return tid, _gmap[tid], len(both)


def plot_demo_tree_partitions(loader: Callable, tid: str, *, title_id: Optional[str] = None):
    """One tree, three partitions: patristic dendrogram + Fiedler-S / Griffing-D / Fiedler-L(D) strips."""
    S, labels, tree, D = loader(tid)
    vS = compute_fiedler_from_laplacian(compute_laplacian(S))     # labels order
    vGD = griffing_leading_eigvec(D)                              # D order (== labels via loader)
    vFD = compute_fiedler_from_laplacian(compute_laplacian(D))    # D order

    to_tree = tree_leaf_index(tree, labels)
    strips = [("Fiedler\n$S$",    sign_orient(vS[to_tree])),
              ("Griffing\n$D$",   sign_orient(vGD[to_tree])),
              ("Fiedler\n$L(D)$", sign_orient(vFD[to_tree]))]

    Dpat = patristic_distance_matrix(tree)
    Z = hierarchy.linkage(squareform(Dpat, checks=False), method="average")

    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(1, 4, width_ratios=[4.5, 1, 1, 1], wspace=0.30)
    ax_tree = fig.add_subplot(gs[0, 0])
    dendro = hierarchy.dendrogram(Z, ax=ax_tree, orientation="left", no_labels=True,
                                  color_threshold=0, above_threshold_color="#262626")
    leaf_order = dendro["leaves"]
    ax_tree.set_title("patristic tree", fontsize=12)
    ax_tree.set_xticks([]); ax_tree.set_yticks([])
    for sp in ax_tree.spines.values():
        sp.set_visible(False)

    cmap = ListedColormap([COLOR_NEG, COLOR_POS])
    for j, (name, sgn) in enumerate(strips, start=1):
        ax = fig.add_subplot(gs[0, j])
        ax.imshow(sgn[leaf_order][:, None], aspect="auto", cmap=cmap, vmin=0, vmax=1,
                  origin="lower")
        n1 = int(sgn.sum()); n2 = len(sgn) - n1
        eta = max(n1, n2) / max(min(n1, n2), 1)
        ax.set_title(name, fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel(f"{max(n1, n2)}/{min(n1, n2)}\n$\\eta$={eta:.1f}", fontsize=10)
    fig.suptitle(f"Tree {title_id if title_id is not None else tid}: partition by method",
                 fontsize=13)
    plt.show()


# ----------------------------------------------------------------------------
# what each operator sees (cells 18/19)
# ----------------------------------------------------------------------------
def plot_operator_matrices(loader: Callable, tid: str):
    """Heatmaps of raw inputs (S, D) and operators (L, B, L_sym), sorted along the Griffing axis."""
    S_m, labels_m, _t_m, D_m = loader(tid)
    L_m = compute_laplacian(S_m)
    B_m = griffing_centered(D_m)
    _deg = np.asarray(S_m.sum(axis=1)).ravel()
    _dis = 1.0 / np.sqrt(np.maximum(_deg, 1e-12))
    Lsym_m = np.eye(len(_deg)) - (S_m * _dis[:, None]) * _dis[None, :]
    Lsym_m = 0.5 * (Lsym_m + Lsym_m.T)

    order = np.argsort(griffing_leading_eigvec(D_m))
    _ro = lambda M: M[order][:, order]
    Sr, Dr, Lr, Br, LSr = _ro(S_m), _ro(D_m), _ro(L_m), _ro(B_m), _ro(Lsym_m)

    DIV = "RdYlBu_r"
    def _sym_range(M):
        off = M[~np.eye(M.shape[0], dtype=bool)]
        return float(np.percentile(np.abs(off), 95)) or 1.0

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for ax, Mx, ttl in [(axes[0, 0], Sr, r"$S$  (similarity)"),
                        (axes[0, 1], Dr, r"$D$  (distance)")]:
        im = ax.imshow(Mx, cmap="viridis")
        ax.set_title(ttl, fontsize=12); fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    axes[0, 2].axis("off")

    for ax, Mx, ttl in [(axes[1, 0], Lr,  r"$L=\mathrm{Deg}(S)-S$"),
                        (axes[1, 1], Br,  r"$B=HDH$"),
                        (axes[1, 2], LSr, r"$L_{\rm sym}=I-D^{-1/2}SD^{-1/2}$")]:
        v = _sym_range(Mx)
        im = ax.imshow(Mx, cmap=DIV, vmin=-v, vmax=v)
        ax.set_title(ttl, fontsize=12); fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    axes[0, 0].set_ylabel("raw inputs", fontsize=13)
    axes[1, 0].set_ylabel("operators", fontsize=13)
    fig.suptitle(f"Tree {tid}", fontsize=14)
    fig.tight_layout(); plt.show()


# ----------------------------------------------------------------------------
# recovery sweep (cells 20/21) + recovery curve
# ----------------------------------------------------------------------------
def _lsym_sweep_one_tree(S, fr_ref, pr_ref, p_values, reps, seed_base):
    pr_int = np.asarray(pr_ref).astype(int)
    nmi = []
    for p_idx, p in enumerate(p_values):
        if p >= 0.9999:
            nmi.append(1.0); continue
        aligned = []
        for r in range(reps):
            S_hat = _subsample_matrix_entries(
                S, p, seed=seed_base + 10_000 * p_idx + r, builder=None)
            try:
                f_est = _fast_fiedler_Lsym(S_hat)
            except Exception:
                continue
            if np.dot(f_est, fr_ref) < 0:
                f_est = -f_est
            aligned.append(f_est)
        if not aligned:
            nmi.append(float("nan")); continue
        v_avg = np.mean(aligned, axis=0)
        nmi.append(float(normalized_mutual_info_score(pr_int, (v_avg > 0).astype(int))))
    return np.asarray(nmi, float)


def run_recovery_sweep(
    loader: Callable,
    cohort: Sequence[str],
    p_values: Sequence[float],
    *,
    reps: int,
    num_gaps: int,
    min_split: int,
    cache_dir=Path("."),
    cache_prefix: str = "",
    cmp_basename: str = "compare_sweep.npz",
    cmp_lsym_basename: str = "compare_sweep_lsym.npz",
) -> dict:
    """Per-tree NMI-vs-p recovery for the 3 methods over a fixed cohort.

    Two caches (S/Griffing together, L_sym separate) keyed on (p_values, cohort).
    Returns ``{nmi_S, nmi_G, nmi_LS, p_values}``.
    """
    cmp = _cache(cache_dir, cache_prefix, cmp_basename)
    cmp_lsym = _cache(cache_dir, cache_prefix, cmp_lsym_basename)
    cohort = list(cohort)
    _pvc = np.round(p_values, 6)

    def _cache_ok(path, key):
        if not Path(path).exists():
            return None
        z = np.load(path, allow_pickle=True)
        if list(np.round(z["p_values"], 6)) == list(_pvc) and list(z["tree_ids"]) == cohort:
            return z[key]
        return None

    nmi_S = _cache_ok(cmp, "nmi_S")
    nmi_G = _cache_ok(cmp, "nmi_G")
    nmi_LS = _cache_ok(cmp_lsym, "nmi_LS")
    need_SG = nmi_S is None or nmi_G is None
    need_LS = nmi_LS is None
    if not (need_SG or need_LS):
        print("loaded recovery sweeps from cache")

    if need_SG or need_LS:
        print(f"computing recovery sweep: {len(cohort)} trees x {len(p_values)} p x {reps} reps "
              f"(S/G={need_SG}, L_sym={need_LS})...")
        _S, _G, _LS = [], [], []
        for ti, tid in enumerate(cohort):
            S, labels, tree, D = loader(tid)
            fr = compute_fiedler_from_laplacian(compute_laplacian(S))
            pr, _, _ = compute_reference_partition_and_quality(
                fr, S, num_gaps=num_gaps, min_split=min_split)
            if need_SG:
                out = bootstrap_p_sweep_simple(
                    S, fr, p_values, bootstrap_reps=reps, seed=1000 * ti,
                    num_gaps=num_gaps, min_split=min_split,
                    partition_ref=pr, partition_method="sigma2")
                _S.append(np.asarray(out["partition_nmi_M"], float))
                raw = bpart_sweep_raw(D, p_values, reps=reps, seed_base=1000 * ti)
                _G.append(np.array([float(np.mean(pp["nmi"])) for pp in raw["per_p"]]))
            if need_LS:
                fr_LS = _fast_fiedler_Lsym(S); pr_LS = (fr_LS > 0).astype(int)
                _LS.append(_lsym_sweep_one_tree(
                    S, fr_LS, pr_LS, p_values, reps, seed_base=1000 * ti + 7777))
            if (ti + 1) % 5 == 0:
                print(f"  swept {ti+1}/{len(cohort)}")
        if need_SG:
            nmi_S = np.vstack(_S); nmi_G = np.vstack(_G)
            np.savez(cmp, p_values=_pvc, tree_ids=np.array(cohort, dtype=object),
                     nmi_S=nmi_S, nmi_G=nmi_G)
        if need_LS:
            nmi_LS = np.vstack(_LS)
            np.savez(cmp_lsym, p_values=_pvc, tree_ids=np.array(cohort, dtype=object),
                     nmi_LS=nmi_LS)
        print("saved recovery caches")

    return {"nmi_S": nmi_S, "nmi_G": nmi_G, "nmi_LS": nmi_LS, "p_values": _pvc}


def plot_recovery_curve(recovery: dict, *, n_taxa: int):
    """Median NMI vs p for the 3 methods (band = +-1 std)."""
    P = np.asarray(recovery["p_values"], float)
    nmi_S, nmi_G, nmi_LS = recovery["nmi_S"], recovery["nmi_G"], recovery["nmi_LS"]
    nt = nmi_S.shape[0]
    fig, ax = plt.subplots(figsize=(7, 7))
    for data, name, col in [
        (nmi_S,  "Fiedler-on-$S$ (similarity)",          COL_S),
        (nmi_G,  "Griffing-on-$D$ (distance)",           COL_D),
        (nmi_LS, "Fiedler-on-$L_{\\rm sym}$ (normalised)", COL_LS),
    ]:
        med = np.nanmedian(data, 0); sd = np.nanstd(data, 0)
        ax.plot(P, med, color=col, lw=2, marker="o", ms=3, label=f"{name} — median")
        ax.fill_between(P, med - sd, med + sd, color=col, alpha=0.15)
    ax.set_xscale("log")
    ax.set_xlabel("sub-sampling fraction $p$ (log)", fontsize=12)
    ax.set_ylabel("NMI vs full-matrix reference", fontsize=12)
    ax.set_title(f"Sub-sampling recovery: similarity vs distance vs normalised\n"
                 f"median NMI over {nt} trees valid in both (band = $\\pm$1 std), n={n_taxa}",
                 fontsize=12)
    ax.set_xlim(P.min() * 0.8, 1.2); ax.set_ylim(-0.05, 1.05)
    ax.legend(fontsize=10, loc="center left")
    fig.tight_layout(); plt.show()


# ----------------------------------------------------------------------------
# p* scatters (cells 21/22 sim, 22/23 real)
# ----------------------------------------------------------------------------
# method palette + per-method markers (real-data style)
_METHOD_MARKERS = [("Fiedler-on-$S$", COL_S, "o"),
                   ("Griffing-on-$D$", COL_D, "s"),
                   ("Fiedler-on-$L_{\\rm sym}$", COL_LS, "^")]


def _pstars(recovery: dict):
    P = np.asarray(recovery["p_values"], float)
    return (P,
            p_star_per_tree(recovery["nmi_S"], P),
            p_star_per_tree(recovery["nmi_G"], P),
            p_star_per_tree(recovery["nmi_LS"], P))


def _panel_pstar(ax, x, methods, *, gen_per_cohort=None, gen_markers=None, legend=False):
    """Scatter + tertile-median lines for p* vs x. If gen_markers given, marker = generator
    and color = method (sim style); else marker per method (real style)."""
    use_gen = gen_per_cohort is not None and gen_markers is not None

    def _scatter(ax, x, y, base_mask, col, mk):
        if use_gen:
            for g, gm in gen_markers.items():
                m = base_mask & (gen_per_cohort == g)
                ax.scatter(x[m], y[m], color=col, marker=gm, alpha=0.55, s=35)
        else:
            ax.scatter(x[base_mask], y[base_mask], color=col, marker=mk, alpha=0.55, s=35,
                       label=None if legend is False else None)

    for entry in methods:
        ps, col, mk = entry[0], entry[1], entry[2]
        msk = np.isfinite(ps) & (ps > 0)
        _scatter(ax, x, ps, msk, col, mk)
    edges = np.quantile(x, [0, 1 / 3, 2 / 3, 1.0]); centers = 0.5 * (edges[:-1] + edges[1:])
    for entry in methods:
        ps, col, mk = entry[0], entry[1], entry[2]
        m = np.isfinite(ps) & (ps > 0)
        meds = [np.median(ps[m & (x >= lo) & (x <= hi)])
                if (m & (x >= lo) & (x <= hi)).any() else np.nan
                for lo, hi in zip(edges[:-1], edges[1:])]
        line_mk = "o" if use_gen else mk
        ax.plot(centers, meds, color=col, lw=2.5, marker=line_mk, ms=10, mec="white", mew=1.5)
    ax.set_yscale("log"); ax.grid(True, alpha=0.3, which="both")


def plot_pstar_vs_eta(
    recovery: dict,
    eta_S_per_cohort,
    eta_D_per_cohort,
    cohort,
    *,
    gen_per_cohort=None,
    gen_markers=None,
):
    """Two-panel p* vs eta_D / eta_S scatter for all three methods.

    Pass ``gen_per_cohort`` + ``gen_markers`` (sim) to mark points by generator;
    omit (real) for one marker per method.
    """
    _, pstar_S, pstar_G, pstar_LS = _pstars(recovery)
    methods = [(pstar_S, COL_S, "o", "Fiedler-on-$S$"),
               (pstar_G, COL_D, "s", "Griffing-on-$D$"),
               (pstar_LS, COL_LS, "^", "Fiedler-on-$L_{\\rm sym}$")]
    eta_S_per_cohort = np.asarray(eta_S_per_cohort, float)
    eta_D_per_cohort = np.asarray(eta_D_per_cohort, float)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    _panel_pstar(a1, eta_D_per_cohort, methods,
                 gen_per_cohort=gen_per_cohort, gen_markers=gen_markers)
    a1.set_xlabel(r"$\eta_D$  (Griffing partition imbalance)", fontsize=12)
    a1.set_ylabel(r"$p^\star$ per tree  (NMI$\geq 0.95$)", fontsize=12)
    a1.set_title(r"$p^\star$ vs $\eta_D$")

    if gen_per_cohort is not None and gen_markers is not None:
        _mh = [Line2D([0], [0], marker="s", ls="", color=col, label=lbl)
               for _, col, _, lbl in methods]
        _gh = [Line2D([0], [0], marker=m, ls="", color="#374151", label=g)
               for g, m in gen_markers.items()]
        _l1 = a1.legend(handles=_mh, fontsize=9, loc="upper left", title="method")
        a1.add_artist(_l1)
        a1.legend(handles=_gh, fontsize=9, loc="lower right", title="generator")
    else:
        _mh = [Line2D([0], [0], marker=mk, ls="", color=col, label=lbl)
               for _, col, mk, lbl in methods]
        a1.legend(handles=_mh, fontsize=9, loc="upper left")

    _panel_pstar(a2, eta_S_per_cohort, methods,
                 gen_per_cohort=gen_per_cohort, gen_markers=gen_markers)
    a2.set_xscale("log")
    a2.set_xlabel(r"$\eta_S$  (Fiedler partition imbalance)", fontsize=12)
    a2.set_title(r"Same data, x-axis = $\eta_S$")
    fig.suptitle(f"$p^\\star$ per tree vs $\\eta$, both-valid cohort (n={len(cohort)})",
                 fontsize=13)
    fig.tight_layout(); plt.show()


def plot_pstar_vs_diameter(
    loader: Callable,
    recovery: dict,
    cohort,
    *,
    cache_dir=Path("."),
    cache_prefix: str = "",
    rt_basename: str = "rT_cohort.npz",
    gen_per_cohort=None,
    gen_markers=None,
):
    """p* vs r(T) = tree diameter under JC distance, for all three methods (cached r(T))."""
    _, pstar_S, pstar_G, pstar_LS = _pstars(recovery)
    methods = [(pstar_S, COL_S, "o", "Fiedler-on-$S$"),
               (pstar_G, COL_D, "s", "Griffing-on-$D$"),
               (pstar_LS, COL_LS, "^", "Fiedler-on-$L_{\\rm sym}$")]
    cohort = list(cohort)

    rt_cache = _cache(cache_dir, cache_prefix, rt_basename)
    rT_per_tree = None
    if rt_cache.exists():
        zr = np.load(rt_cache, allow_pickle=True)
        if list(zr["cohort"]) == cohort:
            rT_per_tree = np.asarray(zr["rT"], float)
            print(f"loaded {rt_cache}")
    if rT_per_tree is None:
        print(f"computing r(T) on {len(cohort)} cohort trees...")
        rT_list = []
        for tid in cohort:
            _, _, _, Dt = loader(tid)
            rT_list.append(float(Dt.max()))
        rT_per_tree = np.asarray(rT_list, float)
        np.savez(rt_cache, cohort=np.array(cohort, dtype=object), rT=rT_per_tree)
        print(f"saved {rt_cache}")

    use_gen = gen_per_cohort is not None and gen_markers is not None
    fig, ax = plt.subplots(figsize=(8, 6))
    for ps, col, mk, lbl in methods:
        m = np.isfinite(ps) & (ps > 0)
        if use_gen:
            for g, gm in gen_markers.items():
                mm = m & (gen_per_cohort == g)
                ax.scatter(rT_per_tree[mm], ps[mm], color=col, marker=gm, alpha=0.55, s=35)
        else:
            ax.scatter(rT_per_tree[m], ps[m], color=col, marker=mk, alpha=0.55, s=35, label=lbl)

    if any((np.isfinite(ps) & (ps > 0)).sum() >= 6 for ps, _, _, _ in methods):
        edges = np.quantile(rT_per_tree, [0, 1 / 3, 2 / 3, 1.0])
        centers = 0.5 * (edges[:-1] + edges[1:])
        for ps, col, mk, _ in methods:
            m = np.isfinite(ps) & (ps > 0)
            meds = [np.median(ps[m & (rT_per_tree >= lo) & (rT_per_tree <= hi)])
                    if (m & (rT_per_tree >= lo) & (rT_per_tree <= hi)).any() else np.nan
                    for lo, hi in zip(edges[:-1], edges[1:])]
            ax.plot(centers, meds, color=col, lw=2.5, marker="o" if use_gen else mk,
                    ms=10, mec="white", mew=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$r(T)$  (tree diameter under JC distance)", fontsize=12)
    ax.set_ylabel(r"$p^\star$  per tree (smallest $p$ with NMI$\geq 0.95$)", fontsize=12)
    ax.set_title(f"$p^\\star$ vs $r(T)$ on the both-valid cohort (n={len(cohort)})\n"
                 f"(thick markers: tertile medians; thin: per-tree scatter)", fontsize=12)
    if use_gen:
        _mh = [Line2D([0], [0], marker="s", ls="", color=col, label=lbl)
               for _, col, _, lbl in methods]
        _gh = [Line2D([0], [0], marker=m, ls="", color="#374151", label=g)
               for g, m in gen_markers.items()]
        _l1 = ax.legend(handles=_mh, fontsize=10, loc="upper left", title="method")
        ax.add_artist(_l1)
        ax.legend(handles=_gh, fontsize=10, loc="lower right", title="generator")
    else:
        ax.legend(fontsize=11, loc="best")
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout(); plt.show()

    print(f"r(T) range: {rT_per_tree.min():.3f} – {rT_per_tree.max():.3f}    "
          f"median = {np.median(rT_per_tree):.3f}")
    for ps, _, _, name in methods:
        m = np.isfinite(ps) & (ps > 0)
        print(f"{name}: p* finite on {int(m.sum())}/{len(cohort)} trees, "
              f"median={np.nanmedian(ps[m]):.3f}" if m.any()
              else f"{name}: p* finite on 0/{len(cohort)} trees")

    return rT_per_tree
