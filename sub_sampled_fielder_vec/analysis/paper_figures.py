"""Single source of truth for the v10 paper's figure provenance.

Every paper figure maps to exactly one producing notebook, the claim it supports,
the cache scope it consumes, and the command that rebuilds it. ``PAPER_MAP.md`` is
*generated* from this file and ``sync_paper_figures.py`` validates against it, so
the docs, the checker and the manifest cannot drift apart.

The previous hand-maintained mapping (inline in ``scripts/sync_paper_figures.py``)
had drifted badly and its own author flagged it as "a guess"
(``open-items/16-C3.md``): it listed the retired eta set {1,2,4,8} instead of the
paper's {1,5,10,15}, omitted the four ``pstar_gen_kmeans_eta*`` panels the paper
actually prints, and attributed ``identity_scatter_kingman.png`` to
``simulation_distance_vs_similarity.ipynb`` when the real producer is
``kingman_threshold_vs_theory.ipynb``. Every entry below was re-derived from the
notebooks' own ``savefig`` calls and cross-checked against ``\\includegraphics`` in
``v10/sections/*.tex``.

When a figure changes owner, edit ONLY this file, then run:
    python scripts/sync_paper_figures.py --check      # validate
    python scripts/sync_paper_figures.py --write-map  # regenerate PAPER_MAP.md
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Figure:
    """One paper figure (or one panel group sharing a float)."""

    number: str          # as printed in the PDF, e.g. "2 (a-d)"
    files: tuple[str, ...]   # paths relative to v10/figures/
    section: str         # the .tex that \includegraphics it
    label: str           # the float's \label
    notebook: str        # producer, relative to analysis/
    claims: tuple[str, ...]  # theorem/prop/cor labels the figure speaks to
    cache: str           # cache scope / results dir it consumes
    rebuild: str         # command(s) to regenerate the inputs
    notes: str = ""


ETAS = (1, 5, 10, 15)

PAPER_FIGURES: tuple[Figure, ...] = (
    Figure(
        number="3 (a-d)",
        files=tuple(f"pstar_synth_nmi_eta{e}.png" for e in ETAS),
        section="sections/8-empirical-study.tex",
        label="fig:pstar_synth",
        notebook="paper/fig02_pstar_synth_cbm.ipynb",
        claims=("thm:main-sim",),
        cache="cache/full_matrix + cache/sweep_trial (via utils/cache.py)",
        rebuild="run the notebook; closed-form CBM, no sweep prerequisite",
        notes="p-hat-star read off at NMI >= 0.90. m must be divisible by (1+eta), so "
              "each eta needs its own m grid.",
    ),
    Figure(
        number="4 (a-d)",
        files=tuple(f"pstar_gen_kmeans_eta{e}.png" for e in ETAS),
        section="sections/8-empirical-study.tex",
        label="fig:pstar_gen",
        notebook="paper/fig03_pstar_gen_kingman.ipynb",
        claims=("thm:main-sim", "cor:nmi"),
        cache="cache/pool_sample + cache/bootstrap_sweep",
        rebuild="python scripts/build_eta_pool.py --n <n> ; "
                "python scripts/build_sweeps.py --ns <n> --methods kmeans",
        notes="k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.",
    ),
    Figure(
        number="2 (a-c)",
        files=tuple(f"Fiedler_Bipartitions/fiedler_tree_partition_eta{e:02d}.png"
                    for e in (1, 5, 15)),
        section="sections/4-assumptions.tex",
        label="fig:fiedler_partitions",
        notebook="paper/fig05_fiedler_partitions.ipynb",
        claims=("thm:stdr-split",),
        cache="cache/pool_sample (n=500 Kingman)",
        rebuild="python scripts/build_eta_pool.py --n 500 ; then run the notebook",
        notes="Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These "
              "used to live outside the manuscript tree and resolved through "
              "\\graphicspath{{../}}, which made the paper non-self-contained. The "
              "eta=10 panel is RETIRED: v10 prints three panels, not four. This is "
              "also where the two plateaus of eq:two_plateau are visible in data.",
    ),
    # Figures 8-11 lived in sections/appendix-emp.tex, deleted 2026-08-04 when the
    # supplementary-figures appendix was removed. Their PNGs are in RETIRED below;
    # the notebooks that build them are untouched and still listed in Supporting.
    Figure(
        number="5 (a-b)",
        files=("eta_hist_kingman_S.png", "eta_hist_bd_S.png"),
        section="sections/E-split-choice.tex",
        label="fig:eta-hist",
        notebook="supporting/eta_by_operator.ipynb",
        claims=("thm:stdr-split",),
        cache="analysis/notebooks_cache/eta_by_operator/{kingman,bd}_n1000_L10000.npz",
        rebuild="run the notebook (it fills the cache on a miss via "
                "analysis/utils/eta_screen.py, ~15 min on 6 workers, resumable); "
                "set FIG_DIR=<scratch> to preview without overwriting v10/figures/",
        notes="No sub-sampling -- this is what the Fiedler rule returns on the FULL "
              "matrix. 1000 Kingman + 1000 birth-death trees, m=1000, L=1e4, k-means on "
              "L_sym as thm:main-sim prescribes. Both panels share bins and axes. The "
              "notebook also emits the B-operator panels; v10 drops the distance route, "
              "so those two PNGs are RETIRED below.",
    ),
)

# Figures that are TikZ/tabular in the .tex -- no image file, so the checker must
# not expect one.
NON_IMAGE_FLOATS = {
    "fig:tree_topology": "Fig 1, TikZ in sections/2-background.tex",
}

# Present in v10/figures/ but deliberately unused. Listed so the checker can tell
# "retired on purpose" from "someone forgot to reference this".
RETIRED: dict[str, str] = {
    # The distance route (B = H*D*H) was cut from the manuscript at v10; its appendix
    # had already been removed, leaving 6 dangling labels in v9.
    **{f"pstar_synth_dist_eta{e}.png": "distance route cut from the manuscript at v10"
       for e in (1, 5, 10, 15)},
    **{f"pstar_gen_dist_eta{e}.png": "distance route cut from the manuscript at v10"
       for e in (1, 5, 10, 15)},
    "eta_hist_kingman_B.png": "distance route cut from the manuscript at v10",
    "eta_hist_bd_B.png": "distance route cut from the manuscript at v10",
    "pstar_balanced_nmi.png": "was Fig 8; sections/appendix-emp.tex deleted 2026-08-04",
    "identity_scatter_kingman.png": "was Fig 9; sections/appendix-emp.tex deleted "
                                    "2026-08-04 -- the generated-regime identity claim "
                                    "at empirical.tex now has no figure behind it "
                                    "(open-items/23-appG-removal.md)",
    "recovery_grid_kmeans.png": "was Fig 10; sections/appendix-emp.tex deleted 2026-08-04",
    "pstar_gen_3operators.png": "was Fig 11; sections/appendix-emp.tex deleted "
                                "2026-08-04 -- the only evidence for choosing k-means "
                                "over the sign and sigma2 rules",
    "pstar_synth_agr_eta1.png": "scored by sign-agreement; paper switched to NMI",
    "pstar_synth_agr_eta2.png": "scored by sign-agreement; paper switched to NMI",
    "pstar_synth_agr_eta4.png": "scored by sign-agreement; paper switched to NMI",
    "pstar_synth_agr_eta8.png": "scored by sign-agreement; paper switched to NMI",
    "pstar_synth_nmi_eta2.png": "retired eta set {1,2,4,8}; paper uses {1,5,10,15}",
    "pstar_synth_nmi_eta4.png": "retired eta set {1,2,4,8}; paper uses {1,5,10,15}",
    "pstar_synth_nmi_eta8.png": "retired eta set {1,2,4,8}; paper uses {1,5,10,15}",
    "pstar_gen_kmeans_2panel.png": "superseded by the four per-eta panels of Fig 3",
    "S_by_alpha.png": "was Fig 7; the HBM appendix was deleted 2026-08-26",
    "spectral_gap_bound.png": "was Fig 7; the HBM appendix was deleted 2026-08-26",
    "Fiedler_Bipartitions/fiedler_tree_partition_eta10.png":
        "Fig 2 went to one row of three etas (1, 5, 15) when it moved into "
        "sections/2-generative-model.tex, 2026-08-25",
    "coherence_vs_eta_by_n.png": "produced by fiedler_tree_partition_by_eta but never "
                                 "placed in the paper (open-items/16-C3.md)",
}


@dataclass(frozen=True)
class Supporting:
    """A notebook that produces no paper figure."""

    notebook: str
    backs: str
    status: str


SUPPORTING: tuple[Supporting, ...] = (
    # The distance route (B = H*D*H) was cut from the manuscript at v10.  These two
    # produced Figs 5-6 of v9; they still run, and their PNGs are in RETIRED.
    Supporting("supporting/pstar_flat_cbm_distance.ipynb",
               "distance-route p* on the closed-form CBM; no v10 figure",
               "frozen; produced v9 Fig 5"),
    Supporting("supporting/pstar_eta_pool_distance.ipynb",
               "distance-route p* on the Kingman eta pool; no v10 figure",
               "frozen; produced v9 Fig 6"),
    Supporting("supporting/distance_vs_similarity_generated.ipynb",
               "App F / thm:main-dist beyond the sweeps of Figs 4-5",
               "live evidence; generated twin of the real-data notebook"),
    Supporting("supporting/distance_vs_similarity_real.ipynb",
               "App F / thm:main-dist on the 600-tree real benchmark",
               "live evidence; needs data/real_datasets (gitignored)"),
    # Deleted 2026-08-08, all three uncited and none reproducible in a fresh clone:
    # nj_subsampling_balanced (snj_subsampling with s/snj/nj/), snj_subsampling
    # (both read gitignored results/runs/ sweeps), and nj_subsampling_nonbalanced --
    # which was not NJ at all but the B=HDH clan-partition prototype, superseded by
    # pstar_{eta_pool,flat_cbm}_distance. subsample_S_vs_subsample_D went too.
    # src/runners/nj_sweep.py and src/plots/plot_{s,}nj_* stay: scripts/ imports them.
    # Was paper/fig04 (Fig 7, S_by_alpha + spectral_gap_bound) until the HBM
    # appendix was deleted on 2026-08-26; nothing in the paper cites the HBM now.
    Supporting("supporting/hbm_spectral_gap.ipynb",
               "the geometric floor S_in^min >= S_in*alpha^Dmax and the gap it induces",
               "frozen; produces no paper figure since 2026-08-26"),
    Supporting("supporting/flat_cbm_variance_recovery.ipynb",
               "asm:margin / cross-clan variance intuition",
               "supporting; frozen"),
    Supporting("supporting/decay_cbm_features.ipynb",
               "HBM identities: mu(U), lambda2, geometric floor (no longer in the paper)",
               "supporting; frozen"),
    Supporting("supporting/nnm_vs_ipw.ipynb",
               "asm:sampling -- why IPW rather than nuclear-norm completion",
               "supporting; outputs cleared, writes to results/notebooks/"),
    Supporting("supporting/stdr_partition_recovery.ipynb",
               "context: one split is not the full recursion (open-items/12-A3.md)",
               "exploratory - not cited; outputs cleared"),
    # Both were paper notebooks (Figs 8 and 9) until sections/appendix-emp.tex was
    # deleted on 2026-08-04. Kept, and moved out of paper/, because the claims they
    # back are still made in the text -- see open-items/23-appG-removal.md.
    Supporting("supporting/pstar_balanced_nmi.ipynb",
               "thm:main-sim at eta=1; was Fig 8 of the removed App G",
               "frozen; its refitted C=32.69 was never carried into any caption"),
    # The birth-death half of eta_by_operator, looked at rather than summarised: three
    # trees across the imbalance range with their matrices, spectra and vectors, then
    # the pooled 200-tree distributions. Reuses eta_by_operator's caches; writes no PNG.
    Supporting("supporting/bd_model_exploration.ipynb",
               "why B returns a valid edge on only 55% of birth-death trees: the "
               "vectors, matrices and spectra behind that number",
               "exploratory - not cited; reads notebooks_cache/eta_by_operator/"),
    Supporting("supporting/identity_checks_kingman.ipynb",
               "prop:coherence / lem:gap in the generated regime; was Fig 9 of the "
               "removed App G, and empirical.tex still asserts those identities",
               "frozen; asserts a gitignored results/ run in its setup cell"),
)


def all_files() -> list[str]:
    """Every figure file the paper prints, relative to v10/figures/."""
    return [f for fig in PAPER_FIGURES for f in fig.files]


def by_file() -> dict[str, Figure]:
    return {f: fig for fig in PAPER_FIGURES for f in fig.files}


def by_notebook() -> dict[str, list[Figure]]:
    out: dict[str, list[Figure]] = {}
    for fig in PAPER_FIGURES:
        out.setdefault(fig.notebook, []).append(fig)
    return out
