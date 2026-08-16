"""Single source of truth for the v9 paper's figure provenance.

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
``v9/sections/*.tex``.

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
    files: tuple[str, ...]   # paths relative to v9/figures/
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
        number="2 (a-d)",
        files=tuple(f"pstar_synth_nmi_eta{e}.png" for e in ETAS),
        section="sections/empirical.tex",
        label="fig:pstar_synth",
        notebook="paper/fig02_pstar_synth_cbm.ipynb",
        claims=("thm:main-sim",),
        cache="cache/full_matrix + cache/sweep_trial (via utils/cache.py)",
        rebuild="run the notebook; closed-form CBM, no sweep prerequisite",
        notes="p-hat-star read off at NMI >= 0.90. m must be divisible by (1+eta), so "
              "each eta needs its own m grid.",
    ),
    Figure(
        number="3 (a-d)",
        files=tuple(f"pstar_gen_kmeans_eta{e}.png" for e in ETAS),
        section="sections/empirical.tex",
        label="fig:pstar_gen",
        notebook="paper/fig03_pstar_gen_kingman.ipynb",
        claims=("thm:main-sim", "cor:nmi"),
        cache="cache/pool_sample + cache/bootstrap_sweep",
        rebuild="python scripts/build_eta_pool.py --n <n> ; "
                "python scripts/build_sweeps.py --ns <n> --methods kmeans",
        notes="k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.",
    ),
    Figure(
        number="4 (a-d)",
        files=tuple(f"pstar_synth_dist_eta{e}.png" for e in ETAS),
        section="sections/empirical.tex",
        label="fig:pstar_synth_dist",
        notebook="supporting/pstar_flat_cbm_distance.ipynb",
        claims=("thm:main-dist", "cor:dist-balanced"),
        cache="cache/bpart_sweep (eigsolver=lm_k1 keyspace)",
        rebuild="run the notebook; closed-form CBM in distance space, no sweep "
                "prerequisite",
        notes="Distance analogue of Figure 2 and LIKE-FOR-LIKE with it (same sign rule, "
              "same per-trial accounting). Only eta=1 is inside thm:main-dist's verified "
              "configuration (cor:dist-balanced). Open markers are p*=1 read-offs, i.e. "
              "never recovered below full data -- censored, excluded from the fitted C. "
              "At eta=15 only 2 of 7 sizes are uncensored, so no C is quoted.",
    ),
    Figure(
        number="5 (a-d)",
        files=tuple(f"pstar_gen_dist_eta{e}.png" for e in ETAS),
        section="sections/empirical.tex",
        label="fig:pstar_gen_dist",
        notebook="supporting/pstar_eta_pool_distance.ipynb",
        claims=("thm:main-dist",),
        cache="cache/pool_sample + cache/bpart_sweep "
              "(eigsolver=lm_k1 + aggregation=avg_vector, 33-point 1e-4 keyspace)",
        rebuild="python scripts/build_eta_pool.py --n <n> ; then run the notebook, which "
                "writes these four PNGs here directly (D = -log S is recovered from the "
                "cached S -- no re-simulation)",
        notes="Replicate accounting is ALIGNED with Figure 3: both sign-align the 10 "
              "sub-sampled eigenvectors, average them, and partition the average once "
              "(aggregation='avg_vector'). Sweeps the WIDE grid geomspace(1e-4, 1, 33), an "
              "exact superset of Figure 3's geomspace(1e-3, 1, 25), so the two stay "
              "readable at matched p; the extra decade was needed because aligning the "
              "accounting pushed eta~1's large-n thresholds below the old 1e-3 floor. "
              "EVERY threshold is now resolved strictly inside the grid -- no floor "
              "censoring anywhere, one ceiling bound (n=500 at eta=10). n=3000 is excluded "
              "(NS_INCLUDE) and is the only size never swept on this grid. C = 0.475 "
              "(spread 2.5x) at eta=1 on 6 of 6 sizes; 2.11 / 19.7 / 24.2 with spreads "
              "166x / 30.5x / 184x at eta=5/10/15, where eta=10,15 rest on 3 pool samples "
              "so those spreads are sample noise. Rounding rule differs from Fig 3 "
              "(k-means there, sign here -- each as its own theorem states it). Kingman "
              "trees are outside thm:main-dist's verified configuration at every eta.",
    ),
    Figure(
        number="6",
        files=("S_by_alpha.png", "spectral_gap_bound.png"),
        section="sections/appendix-D.tex",
        label="fig:hbm_spectral_verification",
        notebook="paper/fig04_hbm_spectral_gap.ipynb",
        claims=("prop:D-gap", "prop:D-floor"),
        cache="none (synthesized HBM matrices, built in-memory)",
        rebuild="run the notebook",
        notes="Both panels MUST come from one run -- a v8 version of this figure "
              "shipped with two panels from two different runs (open-items/00-R1.md), "
              "which is why sync_paper_figures checks for mixed runs.",
    ),
    Figure(
        number="7 (a-d)",
        files=tuple(f"Fiedler_Bipartitions/fiedler_tree_partition_eta{e:02d}.png"
                    for e in ETAS),
        section="sections/appendix-D.tex",
        label="fig:fiedler_partitions",
        notebook="paper/fig05_fiedler_partitions.ipynb",
        claims=("prop:coherence", "lem:gap"),
        cache="cache/pool_sample (n=500 Kingman)",
        rebuild="python scripts/build_eta_pool.py --n 500 ; then run the notebook",
        notes="Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These "
              "used to live outside v9/ and resolved through \\graphicspath{{../}}, "
              "which made the paper non-self-contained.",
    ),
    # Figures 8-11 lived in sections/appendix-emp.tex, deleted 2026-08-04 when the
    # supplementary-figures appendix was removed. Their PNGs are in RETIRED below;
    # the notebooks that build them are untouched and still listed in Supporting.
    Figure(
        number="8 (a-d)",
        files=("eta_hist_kingman_S.png", "eta_hist_kingman_B.png",
               "eta_hist_bd_S.png", "eta_hist_bd_B.png"),
        section="sections/appendix-eta.tex",
        label="fig:eta-hist",
        notebook="supporting/eta_by_operator.ipynb",
        claims=("lem:dist-edge",),
        cache="analysis/notebooks_cache/eta_by_operator/{kingman,bd}_n1000_L10000.npz",
        rebuild="run the notebook (it fills the cache on a miss via "
                "analysis/utils/eta_screen.py, ~15 min on 6 workers, resumable); "
                "set FIG_DIR=<scratch> to preview without overwriting v9/figures/",
        notes="No sub-sampling -- this is what each operator returns on the FULL matrix. "
              "1000 Kingman + 1000 birth-death trees, m=1000, L=1e4. Rules are each "
              "theorem's own: k-means on L(S), sign on B. All four panels share bins and "
              "axes. The birth-death row is the load-bearing one: B's validity collapses "
              "there because it cannot return an imbalanced split -- see app:eta item 3.",
    ),
)

# Figures that are TikZ/tabular in the .tex -- no image file, so the checker must
# not expect one.
NON_IMAGE_FLOATS = {
    "fig:structural_split_balanced": "Fig 1, TikZ in sections/problem.tex",
    "tab:assumption_chain": "Table 1, tabular in sections/problem.tex",
}

# Present in v9/figures/ but deliberately unused. Listed so the checker can tell
# "retired on purpose" from "someone forgot to reference this".
RETIRED: dict[str, str] = {
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
    Supporting("supporting/flat_cbm_variance_recovery.ipynb",
               "asm:margin / cross-clan variance intuition",
               "supporting; frozen"),
    Supporting("supporting/decay_cbm_features.ipynb",
               "App D HBM identities: mu(U), lambda2, geometric floor",
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
    Supporting("supporting/identity_checks_kingman.ipynb",
               "prop:coherence / lem:gap in the generated regime; was Fig 9 of the "
               "removed App G, and empirical.tex still asserts those identities",
               "frozen; asserts a gitignored results/ run in its setup cell"),
)


def all_files() -> list[str]:
    """Every figure file the paper prints, relative to v9/figures/."""
    return [f for fig in PAPER_FIGURES for f in fig.files]


def by_file() -> dict[str, Figure]:
    return {f: fig for fig in PAPER_FIGURES for f in fig.files}


def by_notebook() -> dict[str, list[Figure]]:
    out: dict[str, list[Figure]] = {}
    for fig in PAPER_FIGURES:
        out.setdefault(fig.notebook, []).append(fig)
    return out
