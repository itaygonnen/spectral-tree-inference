"""One generated tree, drawn: its two operators, their spectra, their vectors.

The screens (``eta_screen``, ``real_eta_screen``) keep scalars only. This module
rebuilds a single tree from its id and exposes the objects behind those scalars --
``S``, ``L(S)``, ``D``, ``B = HDH``, the Fiedler vector and the leading-|lambda|
vector -- plus the three drawing routines the notebooks use to look at them.

``build_exhibit`` mirrors ``eta_screen.screen_one`` step for step; the assertions in
the calling notebook check that it reproduces the screened row, which is what makes
the picture evidence about the screen rather than a second, unrelated computation.

Not re-exported from ``analysis.utils.__init__``: notebooks import it by module path,
the same policy as ``sweep_plots``.

``eta_by_operator.ipynb`` still carries its own copy of this code inline. It is a
paper-figure producer and editing a cell clears the stored outputs that hold the
published numbers, so it was left alone; fold it into this module the next time that
notebook has to be re-run anyway.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
import scipy.linalg
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from analysis.utils.eta_screen import MIN_SPLIT, _eta
from src.core.utils import compute_fiedler_from_laplacian, compute_laplacian
from src.models.generated_trees import make_generated
from src.runners.p_sweep_inner import _kmeans_bipartition
from src.utils.griffing import griffing_centered, griffing_leading_eigvec
from src.utils.operator_comparison import patristic_distance_matrix
from src.utils.partition_validity import check_partition_valid_in_tree
from src.utils.screening import _tree_leaf_index

COL_A, COL_B, COL_MIX = "#1d4ed8", "#d97706", "#b8bcc4"   # clan A / clan B / mixed clade
C_MOD = {"kingman": "#0f766e", "bd": "#7c3aed"}           # tree model colours


# --- orientation conventions -------------------------------------------------
def orient(v: np.ndarray, part: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Display convention: clan B (``part == True``) is the minority side.

    Flipping ``v`` and ``part`` together is the same bipartition, so eta and validity
    are unchanged. They MUST flip together -- for arm B the mask *is* sign(v), so for
    that arm this also puts the minority clan on the positive side.
    """
    return (-v, ~part) if part.sum() > part.size / 2 else (v, part)


def orient_minority(v: np.ndarray, part: np.ndarray) -> np.ndarray:
    """Unit norm, scaled by sqrt(m), minority clan on the positive side."""
    v, part = orient(v, part)                     # part == True is now the minority
    if v[part].mean() < v[~part].mean():          # k-means labels carry no sign
        v = -v
    return v / np.linalg.norm(v) * np.sqrt(v.size)


def link_colors(Z: np.ndarray, mask: np.ndarray):
    """Colour every merge in ``Z`` by the clan of the leaves below it.

    A clade is COL_A or COL_B when all of its leaves are on one side of the
    bipartition and COL_MIX when they are not, so a partition that IS a tree edge
    shows as two monochrome subtrees joined by a single grey path to the root.
    ``mask`` must be in the linkage's own leaf order, i.e. tree-leaf order.
    """
    n_b, total, cols = list(mask.astype(int)), [1] * mask.size, {}
    for i, (a, b) in enumerate(Z[:, :2].astype(int)):
        n_b.append(n_b[a] + n_b[b])
        total.append(total[a] + total[b])
        j = i + mask.size
        cols[j] = COL_B if n_b[j] == total[j] else (COL_A if n_b[j] == 0 else COL_MIX)
    return lambda k: cols.get(k, COL_MIX)


# --- the exhibit -------------------------------------------------------------
@dataclass
class Exhibit:
    """Everything one tree needs to be drawn, in two leaf orders.

    ``sidx`` maps labels order -> tree-leaf order (that is the order the validity
    check and the dendrogram want); ``perm`` maps labels order -> display order, so
    ``A[np.ix_(perm, perm)]`` puts a clan of the true tree in a contiguous block.
    """

    tree_id: str
    S: np.ndarray
    D: np.ndarray
    labels: List[str]
    tree: object
    sidx: np.ndarray
    L_S: np.ndarray
    B: np.ndarray
    v_S: np.ndarray
    part_S: np.ndarray
    v_B: np.ndarray
    part_B: np.ndarray
    Z: np.ndarray
    leaf_order: np.ndarray
    perm: np.ndarray

    @property
    def m(self) -> int:
        return int(self.S.shape[0])

    def eta(self, tag: str) -> float:
        return _eta(self.part_S if tag == "S" else self.part_B)

    def valid(self, tag: str) -> bool:
        part = self.part_S if tag == "S" else self.part_B
        return bool(check_partition_valid_in_tree(self.tree, part[self.sidx]))


def build_exhibit(tree_id: str, seq_len: int) -> Exhibit:
    """Rebuild one tree and both operators' vectors; ``make_generated`` is
    deterministic in the id, so this is the same tree the screen saw."""
    S, labels, tree, D = make_generated(tree_id, seq_len=seq_len)
    sidx = _tree_leaf_index(tree, labels)          # labels order -> tree-leaf order

    L_S = compute_laplacian(S)
    v_S = compute_fiedler_from_laplacian(L_S)
    v_S, part_S = orient(v_S, np.asarray(_kmeans_bipartition(v_S, MIN_SPLIT)).astype(bool))
    v_B = griffing_leading_eigvec(D, solver="lm_k1")
    v_B, part_B = orient(v_B, v_B >= 0)

    # Display order: the dendrogram's, from the tree's own patristic distances.
    Z = hierarchy.linkage(squareform(patristic_distance_matrix(tree), checks=False),
                          method="average")
    leaf_order = hierarchy.leaves_list(Z)          # tree-leaf order -> bottom-to-top
    return Exhibit(tree_id, S, D, labels, tree, sidx, L_S, griffing_centered(D),
                   v_S, part_S, v_B, part_B, Z, leaf_order, sidx[leaf_order])


# --- drawing -----------------------------------------------------------------
def draw_tree_pair(ex: Exhibit, figsize=(16, 8)):
    """The tree drawn twice -- once per operator -- with its vector alongside.

    Every clade is coloured by the split that operator read off, so a partition is a
    real tree edge exactly when the grey reduces to a single path to the root.
    """
    n_leaf, y = ex.m, np.arange(ex.m) * 10 + 5     # scipy puts leaf i at 10i+5
    fig = plt.figure(figsize=figsize)
    # explicit margins rather than tight_layout: the dendrogram axes have every spine
    # and tick removed, which is exactly the case tight_layout warns it cannot lay out
    gs = fig.add_gridspec(1, 5, width_ratios=[1, 1, 0.15, 1, 1], wspace=0.07,
                          top=0.80, bottom=0.12, left=0.04, right=0.99)
    handles = ([], [])
    for pair, (v, part, opname, rule, tag) in zip((0, 3), [
            (ex.v_S, ex.part_S, r"Fiedler of $L(S)$", "k-means, k=2", "S"),
            (ex.v_B, ex.part_B, r"leading-$|\lambda|$ of $B=H\mathcal{D}H$",
             r"sign$(v)\geq 0$", "B")]):
        pt = part[ex.sidx]                         # labels order -> tree order

        ax_t = fig.add_subplot(gs[pair])
        with plt.rc_context({"lines.linewidth": 0.4}):
            dendro = hierarchy.dendrogram(ex.Z, ax=ax_t, orientation="left",
                                          no_labels=True,
                                          link_color_func=link_colors(ex.Z, pt))
        lo = np.asarray(dendro["leaves"])          # tree-leaf order -> bottom-to-top
        ax_t.scatter(np.zeros(n_leaf), y, s=4, marker="s", clip_on=False,
                     c=np.where(pt[lo], COL_B, COL_A))      # tip strip
        ax_t.set_ylim(0, 10 * n_leaf)
        ax_t.set_xticks([]); ax_t.set_yticks([])
        for s in ax_t.spines.values():
            s.set_visible(False)
        ax_t.set_title(f"{opname}\n{rule}", fontsize=12)
        ax_t.set_xlabel("the tree, tips and clades by this split", fontsize=10)

        ax_v = fig.add_subplot(gs[pair + 1], sharey=ax_t)
        vv, pp = v[ex.sidx][lo], pt[lo]
        for mask, col, name in ((~pp, COL_A, "clan A"), (pp, COL_B, "clan B")):
            ax_v.scatter(vv[mask], y[mask], s=7, c=col, linewidths=0, alpha=0.8,
                         label=name)
        ax_v.axvline(0, color="#9ca3af", lw=0.6)
        ax_v.set_title(f"$\\eta$={ex.eta(tag):.2f}  "
                       f"({int((~pp).sum())} / {int(pp.sum())})\n"
                       f"{'valid edge' if ex.valid(tag) else 'NOT an edge'}", fontsize=12)
        ax_v.set_xlabel("eigenvector entry", fontsize=10)
        ax_v.tick_params(labelleft=False, labelsize=9)
        for s in ("top", "right", "left"):
            ax_v.spines[s].set_visible(False)
        # a few extreme entries would flatten the bulk onto one vertical line; symlog
        # keeps the bulk linear and compresses those into a tail
        q_lo, q_hi = np.percentile(vv, [1, 99])
        if np.ptp(vv) > 3 * (q_hi - q_lo):
            ax_v.set_xscale("symlog", linthresh=max(abs(q_lo), abs(q_hi)), linscale=2.5)
        handles = ax_v.get_legend_handles_labels()
    fig.legend(*handles, loc="lower center", ncol=2, fontsize=10, markerscale=2.5,
               frameon=False)
    return fig


def matrix_panels(ex: Exhibit):
    """The four operators with the norm each one needs, as draw specs.

    Each panel gets its OWN norm: S lives in (0, 1], L's diagonal is O(m) while its
    off-diagonal is in [-1, 0], D is a distance and B is signed and doubly centred.
    On one shared scale three of the four would be a flat square.
    """
    off = ~np.eye(ex.m, dtype=bool)
    deg = ex.L_S.diagonal()
    # L's diagonal is the degree, orders of magnitude above every other entry; drawn
    # on one scale it turns the panel into a red line on white. Mask it instead.
    L_off = np.where(off, ex.L_S, np.nan)
    lo_S, hi_S = np.percentile(ex.S[off], [1, 99])
    lo_D, hi_D = np.percentile(ex.D[off], [1, 99])
    lo_L, hi_L = np.percentile(ex.L_S[off], [1, 99])
    vmax_B = float(np.percentile(np.abs(ex.B), 99))
    return off, [
        ("S", ex.S, r"$S$   similarity",
         f"linear {lo_S:.3f}-{hi_S:.3f} (1-99 pct off-diag); diag $=1$, clipped",
         "viridis", Normalize(lo_S, hi_S)),
        ("L", L_off, r"$L(S) = \mathrm{Deg}(S) - S$",
         f"linear [{lo_L:.3f}, {hi_L:.3f}]; diag (degree, "
         f"{deg.min():.0f}-{deg.max():.0f}) masked",
         "viridis_r", Normalize(lo_L, hi_L)),
        ("D", ex.D, r"$\mathcal{D} = -\log S$",
         f"linear {lo_D:.2f}-{hi_D:.2f} (1-99 pct off-diag); diag $=0$, clipped",
         "magma_r", Normalize(lo_D, hi_D)),
        ("B", ex.B, r"$B = H\mathcal{D}H$",
         f"diverging $\\pm${vmax_B:.2f} (99 pct of $|B|$); row/col means removed",
         "RdBu_r", Normalize(-vmax_B, vmax_B)),
    ]


def draw_matrices(ex: Exhibit, axes: Sequence, fig, title_size=10):
    """Draw the four operators into four axes, leaves in display order."""
    off, panels = matrix_panels(ex)
    stats = []
    for ax, (key, A, name, how, cmap, norm) in zip(axes, panels):
        cm = plt.get_cmap(cmap).copy()
        cm.set_bad("#e5e7eb")                      # the masked diagonal of L
        im = ax.imshow(A[np.ix_(ex.perm, ex.perm)], cmap=cm, norm=norm,
                       interpolation="nearest", origin="upper")
        ax.set_title(f"{name}\n{how}", fontsize=title_size)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02).ax.tick_params(labelsize=7)
        d = A[off]
        stats.append({"matrix": key, "min": float(np.nanmin(A)), "max": float(np.nanmax(A)),
                      "off-diag mean": float(np.nanmean(d)),
                      "off-diag sd": float(np.nanstd(d))})
    return stats


def leading_spectra(ex: Exhibit, k: int = 50):
    """``(ev_L, ev_B)``: L(S) ascending, B ordered by decreasing |lambda|.

    Each arm is given its eigenvalues in the order it reads them -- L(S) wants the
    second smallest (the Fiedler value), B the largest |lambda|, which for a raw
    distance matrix is the most negative.
    """
    ev_L = scipy.linalg.eigvalsh(ex.L_S, subset_by_index=[0, k - 1])
    ev_B = scipy.linalg.eigvalsh(ex.B)
    return ev_L, ev_B[np.argsort(-np.abs(ev_B))][:k]


def draw_spectra(ex: Exhibit, axes: Sequence, k: int = 50, k_plot: int = 30,
                 zoom: Tuple[int, int] = (1, 10)):
    """The two leading spectra, each on its own scale, with a head zoom.

    ``lambda_1 = 0`` and the bulk sit orders of magnitude away from the gaps that
    matter, so the head of each spectrum needs its own inset to be readable.
    """
    ev_L, ev_B = leading_spectra(ex, k)
    spectra = [
        (ev_L, 1, rf"$L(S)$:  $\lambda_1 \leq \dots \leq \lambda_{{{k_plot}}}$",
         r"$\lambda_2$ (Fiedler)"),
        (ev_B, 0, rf"$B = H\mathcal{{D}}H$:  {k_plot} largest $|\lambda|$",
         r"largest $|\lambda|$"),
    ]
    for ax, (ev, used, name, note) in zip(axes, spectra):
        evp, idx = ev[:k_plot], np.arange(1, k_plot + 1)
        ax.plot(idx, evp, "o-", ms=4, lw=0.8, color="#6b7280", zorder=1)
        ax.scatter([idx[used]], [evp[used]], s=90, facecolor="none",
                   edgecolor="#b91c1c", lw=2, zorder=3)
        ax.axhline(0, color="#9ca3af", lw=0.6)
        ax.set_title(f"{name}\nred circle: {note}", fontsize=11)
        ax.set_xlabel("index", fontsize=10)
        ax.set_ylabel(r"$\lambda$", fontsize=10)
        ax.grid(alpha=0.25)

        a, b = zoom
        axi = ax.inset_axes([0.58, 0.10, 0.40, 0.34])   # flush bottom-right
        axi.plot(idx[a:b], evp[a:b], "o-", ms=4, lw=0.8, color="#6b7280")
        if a <= used < b:                          # B's leading value is outside the zoom
            axi.scatter([idx[used]], [evp[used]], s=70, facecolor="none",
                        edgecolor="#b91c1c", lw=2, zorder=3)
        axi.set_title(rf"zoom: $\lambda_{{{a+1}}}\dots\lambda_{{{b}}}$", fontsize=8)
        axi.set_xticks(np.arange(a + 1, b + 1, 2))
        axi.tick_params(labelsize=7)
        axi.grid(alpha=0.25)
    return ev_L, ev_B


def relative_gaps(ev_L: np.ndarray, ev_B: np.ndarray) -> Tuple[float, float]:
    """``((l3-l2)/l3`` for L(S), ``(|l1|-|l2|)/|l1|`` for B``) -- how isolated each
    eigenvector is, which is the whole story of how stable that split is."""
    gap_L = (ev_L[2] - ev_L[1]) / ev_L[2] if ev_L[2] > 0 else float("nan")
    return float(gap_L), float((abs(ev_B[0]) - abs(ev_B[1])) / abs(ev_B[0]))


__all__ = ["COL_A", "COL_B", "COL_MIX", "C_MOD", "Exhibit", "build_exhibit",
           "draw_matrices", "draw_spectra", "draw_tree_pair", "leading_spectra",
           "link_colors", "matrix_panels", "orient", "orient_minority",
           "relative_gaps"]
