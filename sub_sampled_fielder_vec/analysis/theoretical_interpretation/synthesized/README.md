# synthesized

**Data source:** the similarity `S` (or distance `D`) is built **directly** from a closed-form
block model — flat CBM, hierarchical/decay CBM (HBM), or a balanced binary tree. No tree is
simulated and no sequences are evolved; the matrix is the analytic object under study.

| Topic | Notebook | What it produces |
|---|---|---|
| `cbm_theory/` | `balanced_binary_threshold` | p\* = C·log(n)/n on a symmetric binary tree (α=0.9) |
| `cbm_theory/` | `decay_cbm_features` | HBM identities: μ(U), λ₂, Lemma 0.4 slack across η, m |
| `cbm_theory/` | `flat_cbm_variance_recovery` | Flat-CBM cross-clan variance and sign recovery |
| `cbm_theory/` | `nonbalanced_flat_cbm` | p\* ∝ η(1+η)³ log(m) / [m(ρ − η·S_out)²] |
| `sampling_methods/` | `nnm_vs_ipw` | Soft-Impute (NNM) vs IPW: p\* scaling and runtime |
| `tree_reconstruction/` | `stdr_partition_recovery` | STDR recursive partition: first-layer Fiedler, bipartition Jaccard, ARI |
