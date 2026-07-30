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
        notebook="sec5_empirical/synthesized/nonbalanced_flat_cbm.ipynb",
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
        notebook="sec5_empirical/generated/eta_pool_sweep.ipynb",
        claims=("thm:main-sim", "cor:nmi"),
        cache="cache/pool_sample + cache/bootstrap_sweep",
        rebuild="python scripts/build_eta_pool.py --n <n> ; "
                "python scripts/build_sweeps.py --ns <n> --methods kmeans",
        notes="k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.",
    ),
    Figure(
        number="4",
        files=("S_by_alpha.png", "spectral_gap_bound.png"),
        section="sections/appendix-D.tex",
        label="fig:hbm_spectral_verification",
        notebook="appD_hbm/synthesized/hbm_spectral_gap_verification.ipynb",
        claims=("prop:D-gap", "prop:D-floor"),
        cache="none (synthesized HBM matrices, built in-memory)",
        rebuild="run the notebook",
        notes="Both panels MUST come from one run -- a v8 version of this figure "
              "shipped with two panels from two different runs (open-items/00-R1.md), "
              "which is why sync_paper_figures checks for mixed runs.",
    ),
    Figure(
        number="5 (a-d)",
        files=tuple(f"Fiedler_Bipartitions/fiedler_tree_partition_eta{e:02d}.png"
                    for e in ETAS),
        section="sections/appendix-D.tex",
        label="fig:fiedler_partitions",
        notebook="appD_hbm/generated/fiedler_tree_partition_by_eta.ipynb",
        claims=("prop:coherence", "lem:gap"),
        cache="cache/pool_sample (n=500 Kingman)",
        rebuild="python scripts/build_eta_pool.py --n 500 ; then run the notebook",
        notes="Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These "
              "used to live outside v9/ and resolved through \\graphicspath{{../}}, "
              "which made the paper non-self-contained.",
    ),
    Figure(
        number="6",
        files=("pstar_balanced_nmi.png",),
        section="sections/appendix-emp.tex",
        label="fig:pstar_balanced",
        notebook="appG_supplementary/synthesized/balanced_binary_threshold.ipynb",
        claims=("thm:main-sim",),
        cache="results/notebooks/01_cbm_theory/balanced_binary_threshold/ "
              "(historical cache key -- do not rename)",
        rebuild="run the notebook; reads trials.csv/agg.csv from the scope above",
        notes="Caption still prints C=[TBD, C3]; the refitted value lives only in "
              "open-items/16-C3.md. See open-items/19-cleanup.md.",
    ),
    Figure(
        number="7",
        files=("identity_scatter_kingman.png",),
        section="sections/appendix-emp.tex",
        label="fig:identity_gen",
        notebook="appG_supplementary/generated/kingman_threshold_vs_theory.ipynb",
        claims=("prop:coherence", "lem:gap"),
        cache="results/runs/kingman_mean/uniform/"
              "20260501-194101-kingman_mean_n500-8000_mu_0p1_uniform/",
        rebuild="pipeline A (scripts/run_experiment.py, kingman_mean + uniform) -- "
                "results/ is gitignored, so a fresh clone MUST regenerate this first",
        notes="The notebook asserts that run exists in its setup cell, so it cannot "
              "even start in a fresh clone. Uses kingman_mu0.1 of 3 datasets.",
    ),
    Figure(
        number="8",
        files=("recovery_grid_kmeans.png",),
        section="sections/appendix-emp.tex",
        label="fig:recovery_grid",
        notebook="sec5_empirical/generated/eta_pool_sweep.ipynb",
        claims=("thm:main-sim",),
        cache="cache/pool_sample + cache/bootstrap_sweep",
        rebuild="same as Figure 3 -- one notebook produces Figs 3, 8 and 9",
    ),
    Figure(
        number="9",
        files=("pstar_gen_3operators.png",),
        section="sections/appendix-emp.tex",
        label="fig:operator_sensitivity",
        notebook="sec5_empirical/generated/eta_pool_sweep.ipynb",
        claims=(),
        cache="cache/pool_sample + cache/bootstrap_sweep (all three operators)",
        rebuild="python scripts/build_sweeps.py --ns <n> --methods sign sigma2 kmeans",
        notes="Rounding-rule sensitivity: needs sigma2, the expensive operator. Pass it "
              "as cached-only so a new n cannot trigger a fresh expensive run.",
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
    Supporting("supporting/generated/distance_vs_similarity/simulation_distance_vs_similarity.ipynb",
               "App F / thm:main-dist (the distance route, which has NO figure)",
               "live evidence; generated twin of the real-data notebook"),
    Supporting("supporting/real/distance_vs_similarity/real_data_distance_vs_similarity.ipynb",
               "App F / thm:main-dist on the 600-tree real benchmark",
               "live evidence; needs data/real_datasets (gitignored)"),
    Supporting("supporting/generated/nj_distance/nj_subsampling_balanced.ipynb",
               "distance-route context: NJ under sub-sampling",
               "exploratory - not cited; outputs cleared"),
    Supporting("supporting/generated/nj_distance/nj_subsampling_nonbalanced.ipynb",
               "distance-route context: NJ on unbalanced birth-death trees",
               "exploratory - not cited; cell 15 LAUNCHES an n=4000 sweep"),
    Supporting("supporting/generated/tree_reconstruction/snj_subsampling.ipynb",
               "SNJ sigma2 separation / RF distance vs p",
               "exploratory - not cited; outputs cleared"),
    Supporting("supporting/generated/sampling_methods/distance_vs_similarity.ipynb",
               "sub-sampling S directly vs sub-sampling D then S=exp(-alpha*D)",
               "exploratory - not cited; outputs cleared"),
    Supporting("supporting/synthesized/cbm_theory/flat_cbm_variance_recovery.ipynb",
               "asm:margin / cross-clan variance intuition",
               "supporting; frozen"),
    Supporting("supporting/synthesized/cbm_theory/decay_cbm_features.ipynb",
               "App D HBM identities: mu(U), lambda2, geometric floor",
               "supporting; frozen"),
    Supporting("supporting/synthesized/sampling_methods/nnm_vs_ipw.ipynb",
               "asm:sampling -- why IPW rather than nuclear-norm completion",
               "supporting; outputs cleared, writes to results/notebooks/"),
    Supporting("supporting/synthesized/tree_reconstruction/stdr_partition_recovery.ipynb",
               "context: one split is not the full recursion (open-items/12-A3.md)",
               "exploratory - not cited; outputs cleared"),
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
