# PAPER_MAP — figure provenance for the v10 manuscript

<!-- GENERATED FILE. Edit paper_figures.py, then run:
     python scripts/sync_paper_figures.py --write-map -->

Where every figure in `docs/overleafs/v10/` comes from: which notebook produces it,
which claim it speaks to, which cache or results dir it consumes, and how to rebuild
it. Validated by `python scripts/sync_paper_figures.py --check`.

## Paper figures

| Fig | File(s) | Section | Label | Notebook | Claim(s) |
|---|---|---|---|---|---|
| **3 (a-d)** | `pstar_synth_nmi_eta1.png`<br>`pstar_synth_nmi_eta5.png`<br>`pstar_synth_nmi_eta10.png`<br>`pstar_synth_nmi_eta15.png` | `sections/8-empirical-study.tex` | `fig:pstar_synth` | `paper/fig02_pstar_synth_cbm.ipynb` | `thm:main-sim` |
| **4 (a-d)** | `pstar_gen_kmeans_eta1.png`<br>`pstar_gen_kmeans_eta5.png`<br>`pstar_gen_kmeans_eta10.png`<br>`pstar_gen_kmeans_eta15.png` | `sections/8-empirical-study.tex` | `fig:pstar_gen` | `paper/fig03_pstar_gen_kingman.ipynb` | `thm:main-sim`, `cor:nmi` |
| **2 (a-c)** | `Fiedler_Bipartitions/fiedler_tree_partition_eta01.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta05.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta15.png` | `sections/4-assumptions.tex` | `fig:fiedler_partitions` | `paper/fig05_fiedler_partitions.ipynb` | `thm:stdr-split` |
| **5 (a-b)** | `eta_hist_kingman_S.png`<br>`eta_hist_bd_S.png` | `sections/E-split-choice.tex` | `fig:eta-hist` | `supporting/eta_by_operator.ipynb` | `thm:stdr-split` |

## Inputs and rebuild

### Figure 3 (a-d) — `fig:pstar_synth`

- **Notebook**: `paper/fig02_pstar_synth_cbm.ipynb`
- **Consumes**: cache/full_matrix + cache/sweep_trial (via utils/cache.py)
- **Rebuild**: `run the notebook; closed-form CBM, no sweep prerequisite`
- **Note**: p-hat-star read off at NMI >= 0.90. m must be divisible by (1+eta), so each eta needs its own m grid.

### Figure 4 (a-d) — `fig:pstar_gen`

- **Notebook**: `paper/fig03_pstar_gen_kingman.ipynb`
- **Consumes**: cache/pool_sample + cache/bootstrap_sweep
- **Rebuild**: `python scripts/build_eta_pool.py --n <n> ; python scripts/build_sweeps.py --ns <n> --methods kmeans`
- **Note**: k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.

### Figure 2 (a-c) — `fig:fiedler_partitions`

- **Notebook**: `paper/fig05_fiedler_partitions.ipynb`
- **Consumes**: cache/pool_sample (n=500 Kingman)
- **Rebuild**: `python scripts/build_eta_pool.py --n 500 ; then run the notebook`
- **Note**: Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These used to live outside the manuscript tree and resolved through \graphicspath{{../}}, which made the paper non-self-contained. The eta=10 panel is RETIRED: v10 prints three panels, not four. This is also where the two plateaus of eq:two_plateau are visible in data.

### Figure 5 (a-b) — `fig:eta-hist`

- **Notebook**: `supporting/eta_by_operator.ipynb`
- **Consumes**: analysis/notebooks_cache/eta_by_operator/{kingman,bd}_n1000_L10000.npz
- **Rebuild**: `run the notebook (it fills the cache on a miss via analysis/utils/eta_screen.py, ~15 min on 6 workers, resumable); set FIG_DIR=<scratch> to preview without overwriting v10/figures/`
- **Note**: No sub-sampling -- this is what the Fiedler rule returns on the FULL matrix. 1000 Kingman + 1000 birth-death trees, m=1000, L=1e4, k-means on L_sym as thm:main-sim prescribes. Both panels share bins and axes. The notebook also emits the B-operator panels; v10 drops the distance route, so those two PNGs are RETIRED below.

## Floats with no image file

- `fig:tree_topology` — Fig 1, TikZ in sections/2-background.tex

## Retired figures

Present in `v10/figures/` (or `figures/retired/`) but deliberately unreferenced.
Listed so a stray file reads as "retired on purpose", not "someone forgot".

| File | Why retired |
|---|---|
| `pstar_synth_dist_eta1.png` | distance route cut from the manuscript at v10 |
| `pstar_synth_dist_eta5.png` | distance route cut from the manuscript at v10 |
| `pstar_synth_dist_eta10.png` | distance route cut from the manuscript at v10 |
| `pstar_synth_dist_eta15.png` | distance route cut from the manuscript at v10 |
| `pstar_gen_dist_eta1.png` | distance route cut from the manuscript at v10 |
| `pstar_gen_dist_eta5.png` | distance route cut from the manuscript at v10 |
| `pstar_gen_dist_eta10.png` | distance route cut from the manuscript at v10 |
| `pstar_gen_dist_eta15.png` | distance route cut from the manuscript at v10 |
| `eta_hist_kingman_B.png` | distance route cut from the manuscript at v10 |
| `eta_hist_bd_B.png` | distance route cut from the manuscript at v10 |
| `pstar_balanced_nmi.png` | was Fig 8; sections/appendix-emp.tex deleted 2026-08-04 |
| `identity_scatter_kingman.png` | was Fig 9; sections/appendix-emp.tex deleted 2026-08-04 -- the generated-regime identity claim at empirical.tex now has no figure behind it (open-items/23-appG-removal.md) |
| `recovery_grid_kmeans.png` | was Fig 10; sections/appendix-emp.tex deleted 2026-08-04 |
| `pstar_gen_3operators.png` | was Fig 11; sections/appendix-emp.tex deleted 2026-08-04 -- the only evidence for choosing k-means over the sign and sigma2 rules |
| `pstar_synth_agr_eta1.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta2.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta4.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta8.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_nmi_eta2.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_synth_nmi_eta4.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_synth_nmi_eta8.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_gen_kmeans_2panel.png` | superseded by the four per-eta panels of Fig 3 |
| `S_by_alpha.png` | was Fig 7; the HBM appendix was deleted 2026-08-26 |
| `spectral_gap_bound.png` | was Fig 7; the HBM appendix was deleted 2026-08-26 |
| `Fiedler_Bipartitions/fiedler_tree_partition_eta10.png` | Fig 2 went to one row of three etas (1, 5, 15) when it moved into sections/2-generative-model.tex, 2026-08-25 |
| `coherence_vs_eta_by_n.png` | produced by fiedler_tree_partition_by_eta but never placed in the paper (open-items/16-C3.md) |

## Supporting notebooks (no paper figure)

All 12 are kept. This table exists so nobody mistakes them for figure sources,
and so the evidence behind the distance route stays findable. `thm:main-dist`
does carry figures — Figs 4 and 5 of §6 — but they are produced from
`supporting/`, not from `paper/`, so this table is where their siblings live.

| Notebook | Backs | Status |
|---|---|---|
| `supporting/pstar_flat_cbm_distance.ipynb` | distance-route p* on the closed-form CBM; no v10 figure | frozen; produced v9 Fig 5 |
| `supporting/pstar_eta_pool_distance.ipynb` | distance-route p* on the Kingman eta pool; no v10 figure | frozen; produced v9 Fig 6 |
| `supporting/distance_vs_similarity_generated.ipynb` | App F / thm:main-dist beyond the sweeps of Figs 4-5 | live evidence; generated twin of the real-data notebook |
| `supporting/distance_vs_similarity_real.ipynb` | App F / thm:main-dist on the 600-tree real benchmark | live evidence; needs data/real_datasets (gitignored) |
| `supporting/hbm_spectral_gap.ipynb` | the geometric floor S_in^min >= S_in*alpha^Dmax and the gap it induces | frozen; produces no paper figure since 2026-08-26 |
| `supporting/flat_cbm_variance_recovery.ipynb` | asm:margin / cross-clan variance intuition | supporting; frozen |
| `supporting/decay_cbm_features.ipynb` | HBM identities: mu(U), lambda2, geometric floor (no longer in the paper) | supporting; frozen |
| `supporting/nnm_vs_ipw.ipynb` | asm:sampling -- why IPW rather than nuclear-norm completion | supporting; outputs cleared, writes to results/notebooks/ |
| `supporting/stdr_partition_recovery.ipynb` | context: one split is not the full recursion (open-items/12-A3.md) | exploratory - not cited; outputs cleared |
| `supporting/pstar_balanced_nmi.ipynb` | thm:main-sim at eta=1; was Fig 8 of the removed App G | frozen; its refitted C=32.69 was never carried into any caption |
| `supporting/bd_model_exploration.ipynb` | why B returns a valid edge on only 55% of birth-death trees: the vectors, matrices and spectra behind that number | exploratory - not cited; reads notebooks_cache/eta_by_operator/ |
| `supporting/identity_checks_kingman.ipynb` | prop:coherence / lem:gap in the generated regime; was Fig 9 of the removed App G, and empirical.tex still asserts those identities | frozen; asserts a gitignored results/ run in its setup cell |

