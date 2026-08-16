# PAPER_MAP — figure provenance for the v9 manuscript

<!-- GENERATED FILE. Edit paper_figures.py, then run:
     python scripts/sync_paper_figures.py --write-map -->

Where every figure in `docs/overleafs/v9/` comes from: which notebook produces it,
which claim it speaks to, which cache or results dir it consumes, and how to rebuild
it. Validated by `python scripts/sync_paper_figures.py --check`.

## Paper figures

| Fig | File(s) | Section | Label | Notebook | Claim(s) |
|---|---|---|---|---|---|
| **2 (a-d)** | `pstar_synth_nmi_eta1.png`<br>`pstar_synth_nmi_eta5.png`<br>`pstar_synth_nmi_eta10.png`<br>`pstar_synth_nmi_eta15.png` | `sections/empirical.tex` | `fig:pstar_synth` | `paper/fig02_pstar_synth_cbm.ipynb` | `thm:main-sim` |
| **3 (a-d)** | `pstar_gen_kmeans_eta1.png`<br>`pstar_gen_kmeans_eta5.png`<br>`pstar_gen_kmeans_eta10.png`<br>`pstar_gen_kmeans_eta15.png` | `sections/empirical.tex` | `fig:pstar_gen` | `paper/fig03_pstar_gen_kingman.ipynb` | `thm:main-sim`, `cor:nmi` |
| **4 (a-d)** | `pstar_synth_dist_eta1.png`<br>`pstar_synth_dist_eta5.png`<br>`pstar_synth_dist_eta10.png`<br>`pstar_synth_dist_eta15.png` | `sections/empirical.tex` | `fig:pstar_synth_dist` | `supporting/pstar_flat_cbm_distance.ipynb` | `thm:main-dist`, `cor:dist-balanced` |
| **5 (a-d)** | `pstar_gen_dist_eta1.png`<br>`pstar_gen_dist_eta5.png`<br>`pstar_gen_dist_eta10.png`<br>`pstar_gen_dist_eta15.png` | `sections/empirical.tex` | `fig:pstar_gen_dist` | `supporting/pstar_eta_pool_distance.ipynb` | `thm:main-dist` |
| **6** | `S_by_alpha.png`<br>`spectral_gap_bound.png` | `sections/appendix-D.tex` | `fig:hbm_spectral_verification` | `paper/fig04_hbm_spectral_gap.ipynb` | `prop:D-gap`, `prop:D-floor` |
| **7 (a-d)** | `Fiedler_Bipartitions/fiedler_tree_partition_eta01.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta05.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta10.png`<br>`Fiedler_Bipartitions/fiedler_tree_partition_eta15.png` | `sections/appendix-D.tex` | `fig:fiedler_partitions` | `paper/fig05_fiedler_partitions.ipynb` | `prop:coherence`, `lem:gap` |
| **8 (a-d)** | `eta_hist_kingman_S.png`<br>`eta_hist_kingman_B.png`<br>`eta_hist_bd_S.png`<br>`eta_hist_bd_B.png` | `sections/appendix-eta.tex` | `fig:eta-hist` | `supporting/eta_by_operator.ipynb` | `lem:dist-edge` |

## Inputs and rebuild

### Figure 2 (a-d) — `fig:pstar_synth`

- **Notebook**: `paper/fig02_pstar_synth_cbm.ipynb`
- **Consumes**: cache/full_matrix + cache/sweep_trial (via utils/cache.py)
- **Rebuild**: `run the notebook; closed-form CBM, no sweep prerequisite`
- **Note**: p-hat-star read off at NMI >= 0.90. m must be divisible by (1+eta), so each eta needs its own m grid.

### Figure 3 (a-d) — `fig:pstar_gen`

- **Notebook**: `paper/fig03_pstar_gen_kingman.ipynb`
- **Consumes**: cache/pool_sample + cache/bootstrap_sweep
- **Rebuild**: `python scripts/build_eta_pool.py --n <n> ; python scripts/build_sweeps.py --ns <n> --methods kmeans`
- **Note**: k-means on L_sym (PAPER_METHOD). Kingman + JC69, l=1e4.

### Figure 4 (a-d) — `fig:pstar_synth_dist`

- **Notebook**: `supporting/pstar_flat_cbm_distance.ipynb`
- **Consumes**: cache/bpart_sweep (eigsolver=lm_k1 keyspace)
- **Rebuild**: `run the notebook; closed-form CBM in distance space, no sweep prerequisite`
- **Note**: Distance analogue of Figure 2 and LIKE-FOR-LIKE with it (same sign rule, same per-trial accounting). Only eta=1 is inside thm:main-dist's verified configuration (cor:dist-balanced). Open markers are p*=1 read-offs, i.e. never recovered below full data -- censored, excluded from the fitted C. At eta=15 only 2 of 7 sizes are uncensored, so no C is quoted.

### Figure 5 (a-d) — `fig:pstar_gen_dist`

- **Notebook**: `supporting/pstar_eta_pool_distance.ipynb`
- **Consumes**: cache/pool_sample + cache/bpart_sweep (eigsolver=lm_k1 + aggregation=avg_vector, 33-point 1e-4 keyspace)
- **Rebuild**: `python scripts/build_eta_pool.py --n <n> ; then run the notebook, which writes these four PNGs here directly (D = -log S is recovered from the cached S -- no re-simulation)`
- **Note**: Replicate accounting is ALIGNED with Figure 3: both sign-align the 10 sub-sampled eigenvectors, average them, and partition the average once (aggregation='avg_vector'). Sweeps the WIDE grid geomspace(1e-4, 1, 33), an exact superset of Figure 3's geomspace(1e-3, 1, 25), so the two stay readable at matched p; the extra decade was needed because aligning the accounting pushed eta~1's large-n thresholds below the old 1e-3 floor. EVERY threshold is now resolved strictly inside the grid -- no floor censoring anywhere, one ceiling bound (n=500 at eta=10). n=3000 is excluded (NS_INCLUDE) and is the only size never swept on this grid. C = 0.475 (spread 2.5x) at eta=1 on 6 of 6 sizes; 2.11 / 19.7 / 24.2 with spreads 166x / 30.5x / 184x at eta=5/10/15, where eta=10,15 rest on 3 pool samples so those spreads are sample noise. Rounding rule differs from Fig 3 (k-means there, sign here -- each as its own theorem states it). Kingman trees are outside thm:main-dist's verified configuration at every eta.

### Figure 6 — `fig:hbm_spectral_verification`

- **Notebook**: `paper/fig04_hbm_spectral_gap.ipynb`
- **Consumes**: none (synthesized HBM matrices, built in-memory)
- **Rebuild**: `run the notebook`
- **Note**: Both panels MUST come from one run -- a v8 version of this figure shipped with two panels from two different runs (open-items/00-R1.md), which is why sync_paper_figures checks for mixed runs.

### Figure 7 (a-d) — `fig:fiedler_partitions`

- **Notebook**: `paper/fig05_fiedler_partitions.ipynb`
- **Consumes**: cache/pool_sample (n=500 Kingman)
- **Rebuild**: `python scripts/build_eta_pool.py --n 500 ; then run the notebook`
- **Note**: Writes via PAPER_FIG_DIR into figures/Fiedler_Bipartitions/. These used to live outside v9/ and resolved through \graphicspath{{../}}, which made the paper non-self-contained.

### Figure 8 (a-d) — `fig:eta-hist`

- **Notebook**: `supporting/eta_by_operator.ipynb`
- **Consumes**: analysis/notebooks_cache/eta_by_operator/{kingman,bd}_n1000_L10000.npz
- **Rebuild**: `run the notebook (it fills the cache on a miss via analysis/utils/eta_screen.py, ~15 min on 6 workers, resumable); set FIG_DIR=<scratch> to preview without overwriting v9/figures/`
- **Note**: No sub-sampling -- this is what each operator returns on the FULL matrix. 1000 Kingman + 1000 birth-death trees, m=1000, L=1e4. Rules are each theorem's own: k-means on L(S), sign on B. All four panels share bins and axes. The birth-death row is the load-bearing one: B's validity collapses there because it cannot return an imbalanced split -- see app:eta item 3.

## Floats with no image file

- `fig:structural_split_balanced` — Fig 1, TikZ in sections/problem.tex
- `tab:assumption_chain` — Table 1, tabular in sections/problem.tex

## Retired figures

Present in `v9/figures/` (or `figures/retired/`) but deliberately unreferenced.
Listed so a stray file reads as "retired on purpose", not "someone forgot".

| File | Why retired |
|---|---|
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
| `coherence_vs_eta_by_n.png` | produced by fiedler_tree_partition_by_eta but never placed in the paper (open-items/16-C3.md) |

## Supporting notebooks (no paper figure)

All 8 are kept. This table exists so nobody mistakes them for figure sources,
and so the evidence behind the distance route stays findable. `thm:main-dist`
does carry figures — Figs 4 and 5 of §6 — but they are produced from
`supporting/`, not from `paper/`, so this table is where their siblings live.

| Notebook | Backs | Status |
|---|---|---|
| `supporting/distance_vs_similarity_generated.ipynb` | App F / thm:main-dist beyond the sweeps of Figs 4-5 | live evidence; generated twin of the real-data notebook |
| `supporting/distance_vs_similarity_real.ipynb` | App F / thm:main-dist on the 600-tree real benchmark | live evidence; needs data/real_datasets (gitignored) |
| `supporting/flat_cbm_variance_recovery.ipynb` | asm:margin / cross-clan variance intuition | supporting; frozen |
| `supporting/decay_cbm_features.ipynb` | App D HBM identities: mu(U), lambda2, geometric floor | supporting; frozen |
| `supporting/nnm_vs_ipw.ipynb` | asm:sampling -- why IPW rather than nuclear-norm completion | supporting; outputs cleared, writes to results/notebooks/ |
| `supporting/stdr_partition_recovery.ipynb` | context: one split is not the full recursion (open-items/12-A3.md) | exploratory - not cited; outputs cleared |
| `supporting/pstar_balanced_nmi.ipynb` | thm:main-sim at eta=1; was Fig 8 of the removed App G | frozen; its refitted C=32.69 was never carried into any caption |
| `supporting/identity_checks_kingman.ipynb` | prop:coherence / lem:gap in the generated regime; was Fig 9 of the removed App G, and empirical.tex still asserts those identities | frozen; asserts a gitignored results/ run in its setup cell |

