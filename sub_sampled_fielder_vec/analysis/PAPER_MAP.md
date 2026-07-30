# PAPER_MAP — figure provenance for the v9 manuscript

<!-- GENERATED FILE. Edit paper_figures.py, then run:
     python scripts/sync_paper_figures.py --write-map -->

Where every figure in `docs/overleafs/v9/` comes from: which notebook produces it,
which claim it speaks to, which cache or results dir it consumes, and how to rebuild
it. Validated by `python scripts/sync_paper_figures.py --check`.

## Paper figures

| Fig | File(s) | Section | Label | Notebook | Claim(s) |
|---|---|---|---|---|---|
| **2 (a-d)** | `pstar_synth_nmi_eta1.png`<br>`pstar_synth_nmi_eta5.png`<br>`pstar_synth_nmi_eta10.png`<br>`pstar_synth_nmi_eta15.png` | `sections/empirical.tex` | `fig:pstar_synth` | `sec5_empirical/synthesized/nonbalanced_flat_cbm.ipynb` | `thm:main-sim` |
| **3 (a-d)** | `pstar_gen_kmeans_eta1.png`<br>`pstar_gen_kmeans_eta5.png`<br>`pstar_gen_kmeans_eta10.png`<br>`pstar_gen_kmeans_eta15.png` | `sections/empirical.tex` | `fig:pstar_gen` | `sec5_empirical/generated/eta_pool_sweep.ipynb` | `thm:main-sim`, `cor:nmi` |
| **4** | `S_by_alpha.png`<br>`spectral_gap_bound.png` | `sections/appendix-D.tex` | `fig:hbm_spectral_verification` | `appD_hbm/synthesized/hbm_spectral_gap_verification.ipynb` | `prop:D-gap`, `prop:D-floor` |
| **5 (a-d)** | `Fiedler_Bipartitions/fiedler_tree_partition_eta01.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta05.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta10.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta15.png` | `sections/appendix-D.tex` | `fig:fiedler_partitions` | `appD_hbm/generated/fiedler_tree_partition_by_eta.ipynb` | `prop:coherence`, `lem:gap` |
| **6** | `pstar_balanced_nmi.png` | `sections/appendix-emp.tex` | `fig:pstar_balanced` | `appG_supplementary/synthesized/balanced_binary_threshold.ipynb` | `thm:main-sim` |
| **7** | `identity_scatter_kingman.png` | `sections/appendix-emp.tex` | `fig:identity_gen` | `appG_supplementary/generated/kingman_threshold_vs_theory.ipynb` | `prop:coherence`, `lem:gap` |
| **8** | `recovery_grid_kmeans.png` | `sections/appendix-emp.tex` | `fig:recovery_grid` | `sec5_empirical/generated/eta_pool_sweep.ipynb` | `thm:main-sim` |
| **9** | `pstar_gen_3operators.png` | `sections/appendix-emp.tex` | `fig:operator_sensitivity` | `sec5_empirical/generated/eta_pool_sweep.ipynb` | — |

## Inputs and rebuild

### Figure 2 (a-d) — `fig:pstar_synth`

- **Notebook**: `sec5_empirical/synthesized/nonbalanced_flat_cbm.ipynb`
- **Consumes**: cache/full_matrix + cache/sweep_trial (via utils/cache.py)
- **Rebuild**: `run the notebook; closed-form CBM, no sweep prerequisite`
- **Note**: p-hat-star read off at NMI >= 0.90. m must be divisible by (1+eta), so each eta needs its own m grid.

### Figure 3 (a-d) — `fig:pstar_gen`

- **Notebook**: `sec5_empirical/generated/eta_pool_sweep.ipynb`
- **Consumes**: cache/pool_sample + cache/bootstrap_sweep
- **Rebuild**: `python scripts/build_eta_pool.py --n <n> ; python scripts/build_sweeps.py --ns <n> --methods kmeans`
- **Note**: k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.

### Figure 4 — `fig:hbm_spectral_verification`

- **Notebook**: `appD_hbm/synthesized/hbm_spectral_gap_verification.ipynb`
- **Consumes**: none (synthesized HBM matrices, built in-memory)
- **Rebuild**: `run the notebook`
- **Note**: Both panels MUST come from one run -- a v8 version of this figure shipped with two panels from two different runs (open-items/00-R1.md), which is why sync_paper_figures checks for mixed runs.

### Figure 5 (a-d) — `fig:fiedler_partitions`

- **Notebook**: `appD_hbm/generated/fiedler_tree_partition_by_eta.ipynb`
- **Consumes**: cache/pool_sample (n=500 Kingman)
- **Rebuild**: `python scripts/build_eta_pool.py --n 500 ; then run the notebook`
- **Note**: Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These used to live outside v9/ and resolved through \graphicspath{{../}}, which made the paper non-self-contained.

### Figure 6 — `fig:pstar_balanced`

- **Notebook**: `appG_supplementary/synthesized/balanced_binary_threshold.ipynb`
- **Consumes**: results/notebooks/01_cbm_theory/balanced_binary_threshold/ (historical cache key -- do not rename)
- **Rebuild**: `run the notebook; reads trials.csv/agg.csv from the scope above`
- **Note**: Caption still prints C=[TBD, C3]; the refitted value lives only in open-items/16-C3.md. See open-items/19-cleanup.md.

### Figure 7 — `fig:identity_gen`

- **Notebook**: `appG_supplementary/generated/kingman_threshold_vs_theory.ipynb`
- **Consumes**: results/runs/kingman_mean/uniform/20260501-194101-kingman_mean_n500-8000_mu_0p1_uniform/
- **Rebuild**: `pipeline A (scripts/run_experiment.py, kingman_mean + uniform) -- results/ is gitignored, so a fresh clone MUST regenerate this first`
- **Note**: The notebook asserts that run exists in its setup cell, so it cannot even start in a fresh clone. Uses kingman_mu0.1 of 3 datasets.

### Figure 8 — `fig:recovery_grid`

- **Notebook**: `sec5_empirical/generated/eta_pool_sweep.ipynb`
- **Consumes**: cache/pool_sample + cache/bootstrap_sweep
- **Rebuild**: `same as Figure 3 -- one notebook produces Figs 3, 8 and 9`

### Figure 9 — `fig:operator_sensitivity`

- **Notebook**: `sec5_empirical/generated/eta_pool_sweep.ipynb`
- **Consumes**: cache/pool_sample + cache/bootstrap_sweep (all three operators)
- **Rebuild**: `python scripts/build_sweeps.py --ns <n> --methods sign sigma2 kmeans`
- **Note**: Rounding-rule sensitivity: needs sigma2, the expensive operator. Pass it as cached-only so a new n cannot trigger a fresh expensive run.

## Floats with no image file

- `fig:structural_split_balanced` — Fig 1, TikZ in sections/problem.tex
- `tab:assumption_chain` — Table 1, tabular in sections/problem.tex

## Retired figures

Present in `v9/figures/` (or `figures/retired/`) but deliberately unreferenced.
Listed so a stray file reads as "retired on purpose", not "someone forgot".

| File | Why retired |
|---|---|
| `pstar_synth_agr_eta1.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta2.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta4.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_agr_eta8.png` | scored by sign-agreement; paper switched to NMI |
| `pstar_synth_nmi_eta2.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_synth_nmi_eta4.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_synth_nmi_eta8.png` | retired eta set {1,2,4,8}; paper uses {1,5,10,15} |
| `pstar_gen_kmeans_2panel.png` | superseded by the four per-eta panels of Fig 3 |
| `coherence_vs_eta_by_n.png` | produced by fiedler_tree_partition_by_eta but never placed in the paper (open-items/16-C3.md) |

## Supporting notebooks (no paper figure)

All 10 are kept. This table exists so nobody mistakes them for figure sources,
and so the evidence behind the distance route stays findable — App F /
`thm:main-dist` carries no figure of its own.

| Notebook | Backs | Status |
|---|---|---|
| `supporting/generated/distance_vs_similarity/simulation_distance_vs_similarity.ipynb` | App F / thm:main-dist (the distance route, which has NO figure) | live evidence; generated twin of the real-data notebook |
| `supporting/real/distance_vs_similarity/real_data_distance_vs_similarity.ipynb` | App F / thm:main-dist on the 600-tree real benchmark | live evidence; needs data/real_datasets (gitignored) |
| `supporting/generated/nj_distance/nj_subsampling_balanced.ipynb` | distance-route context: NJ under sub-sampling | exploratory - not cited; outputs cleared |
| `supporting/generated/nj_distance/nj_subsampling_nonbalanced.ipynb` | distance-route context: NJ on unbalanced birth-death trees | exploratory - not cited; cell 15 LAUNCHES an n=4000 sweep |
| `supporting/generated/tree_reconstruction/snj_subsampling.ipynb` | SNJ sigma2 separation / RF distance vs p | exploratory - not cited; outputs cleared |
| `supporting/generated/sampling_methods/distance_vs_similarity.ipynb` | sub-sampling S directly vs sub-sampling D then S=exp(-alpha*D) | exploratory - not cited; outputs cleared |
| `supporting/synthesized/cbm_theory/flat_cbm_variance_recovery.ipynb` | asm:margin / cross-clan variance intuition | supporting; frozen |
| `supporting/synthesized/cbm_theory/decay_cbm_features.ipynb` | App D HBM identities: mu(U), lambda2, geometric floor | supporting; frozen |
| `supporting/synthesized/sampling_methods/nnm_vs_ipw.ipynb` | asm:sampling -- why IPW rather than nuclear-norm completion | supporting; outputs cleared, writes to results/notebooks/ |
| `supporting/synthesized/tree_reconstruction/stdr_partition_recovery.ipynb` | context: one split is not the full recursion (open-items/12-A3.md) | exploratory - not cited; outputs cleared |

