## CLEANUP — repository ordering pass (2026-07-30)

Fragment owner: repo-hygiene pass over `sub_sampled_fielder_vec/`, ahead of merge.

Scope: reorganized code, restored the paper to version control, built a machine-checked
figure→claim map. **No mathematics, figure, or plotted curve was changed** — all 23
pre-existing figure md5s are unchanged and the PDF still builds to 29 pages /
2,586,870 bytes. Every entry below is referred to the author, not decided by the pass.

### [CLEANUP/01] Two figure constants still print as `[TBD]` in the compiled PDF
- **Anchor:** `empirical.tex:64` (`$C=[\text{TBD}]$`, Fig 2), `appendix-emp.tex:45` (`$C=[\text{TBD, C3}]$`, Fig 6), `open-items/16-C3.md` `[C3/01]`
- **Type:** flaw
- **Item:** The refitted constant `C = 32.69` exists only as prose in the register and was never carried into the `.tex`. A reader of the current PDF sees the literal placeholder in two captions.
- **Blocking:** yes — a shipped PDF cannot contain `[TBD]`.

### [CLEANUP/02] Two non-comparable estimators of `C` are both live
- **Anchor:** `analysis/theoretical_interpretation/utils/sweep_plots.py` (`_theory_ref`), `utils/sweep_plots_two_panel.py` (`median_ratio_C`); register `[C2/xx]`
- **Type:** doubt
- **Item:** `_theory_ref` fits `C` in `p* = C log n / n` by least squares and feeds the appendix figures (Fig 8, Fig 9), while `median_ratio_C` uses the median of per-point ratios and feeds the main figures (Fig 2, Fig 3). A linear-scale LS fit is dominated by the largest `p*`, i.e. by the smallest trees, so constants from the two modules must not be quoted against each other. Relatedly, the free exponent `b` in Fig 9's `n^{b}` labels is scaffolding for the operator comparison and should not be presented as one of the paper's metrics.
- **Blocking:** no — default taken: both module headers now state which figures they serve and that the constants differ; neither curve was altered, since a plotted curve is a paper claim.

### [CLEANUP/03] `deferred/` is empty but the master says Related Work lives there
- **Anchor:** `thesis_v9.tex:27`, `sections/related.tex`, `v9/deferred/`
- **Type:** flaw
- **Item:** The directory contains no files; the outline is at `sections/related.tex` and is self-labelled "DELIVERABLE STATUS: POINTS OUTLINE, not finished prose".
- **Blocking:** no.

### [CLEANUP/04] The commented intro `\input` names a nonexistent file
- **Anchor:** `thesis_v9.tex:113`, `sections/intro_prev.tex`
- **Type:** flaw
- **Item:** `\input{sections/intro}` does not match the file on disk, `intro_prev.tex`, so uncommenting as written breaks the build. The paper currently ships with no Introduction and no Conclusion, the latter already on record in `open-items/13-B1.md`.
- **Blocking:** yes if an Introduction is intended for the submitted version.

### [CLEANUP/05] Appendix filenames do not match their compiled letters
- **Anchor:** `sections/appendix-G.tex` (compiles as App F), `sections/appendix-emp.tex` (compiles as App G)
- **Type:** dilemma
- **Item:** Anyone asked to edit "Appendix F" will open `appendix-G.tex` and edit the wrong appendix. The notebook directories deliberately follow the compiled letters, so `appG_supplementary/` currently corresponds to `appendix-emp.tex`.
- **Options:** (a) rename the `.tex` files to match compiled letters, then rename `appG_supplementary/` accordingly; (b) leave both and add a mapping comment in `thesis_v9.tex`; (c) drop letter-based names entirely in favour of topic names (`appendix-dist.tex`, `appendix-emp.tex`).
- **Blocking:** no.

### [CLEANUP/06] Master header comment describes the pre-deferral section numbering
- **Anchor:** `thesis_v9.tex:7-24`
- **Type:** flaw
- **Item:** The comment block still lists §1 Introduction through §6 Empirical, but the compiled numbering is §1 Problem through §5 Empirical.
- **Blocking:** no.

### [CLEANUP/07] The distance route has no figure
- **Anchor:** `sections/appendix-G.tex` (App F), `thm:main-dist`, `analysis/theoretical_interpretation/PAPER_MAP.md`
- **Type:** direction
- **Item:** Theorem 2 and its entire appendix carry no figure, while four supporting notebooks produce directly relevant evidence. Those notebooks are now labelled and indexed in `PAPER_MAP.md`, so the evidence is findable if a reviewer asks why the second main theorem is unillustrated.
- **Blocking:** no.

### [CLEANUP/08] Figure 7 cannot be regenerated from a clean checkout
- **Anchor:** `appG_supplementary/generated/kingman_threshold_vs_theory.ipynb` setup cell
- **Type:** flaw
- **Item:** The notebook hard-asserts `results/runs/kingman_mean/uniform/20260501-194101-…/n500_L10000/results.json`, and `results/` is gitignored, so in a fresh clone it fails before any compute. The 4.4 GB of results was deliberately left uncommitted, so the prerequisite sweep must be re-run first; this is recorded in `PAPER_MAP.md` and `docs/RUNBOOK.md`.
- **Blocking:** no — but it is the one paper figure not reproducible from a clean checkout.

### [CLEANUP/09] Register fragments reference pre-rename notebook paths
- **Anchor:** `open-items/15-C2.md`, `open-items/16-C3.md`
- **Type:** question
- **Item:** Notebooks moved from `theoretical_interpretation/{synthesized,generated,real}/` into `sec5_empirical/`, `appD_hbm/`, `appG_supplementary/` and `supporting/`, and the fragments were not rewritten because the register is a historical record of what each agent did. If it is meant to stay navigable, it needs a path-translation note or a sweep.
- **Blocking:** no — default taken: paths left as written; `PAPER_MAP.md` is the current authority for figure provenance.

### [CLEANUP/10] Restating the discarded-Neumann-remainder defect so the merge cannot lose it
- **Anchor:** `DISTANCE_REVIEW.md` §7.1, `open-items/13-B1.md`
- **Type:** flaw
- **Item:** The discarded Neumann remainder is reported to break the stated rate on both the similarity and the distance route. Nothing in this pass touched the mathematics.
- **Blocking:** yes, per the existing review.
