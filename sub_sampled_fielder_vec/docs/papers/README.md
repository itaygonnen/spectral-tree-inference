# Reference papers

| File | Topic | Where it's used |
|---|---|---|
| `Spectral_top-down_recovery_of_latent_tree_models.pdf` | STDR / top-down spectral recovery | Backs `spectraltree/spectral_tree_reconstruction.py:STDR`. **The only PDF tracked in git.** |
| `Matrix Completion from a Few Entries copy` | leveraged matrix completion | Backs `src/core/sampling/leveraged/`. |
| `Universal Matrix Completion copy.pdf` | matrix completion bounds | ditto |
| `Noise Thresholds for Spectral Clustering copy.pdf.pdf` | spectral-clustering noise thresholds | context for the recovery threshold `p*` |
| `Unperturbed- spectral analysis beyond Davis-Kahan copy.pdf` | entry-wise eigenvector analysis | context for the Neumann route (App C) |
| `The Neighbor-joining Method copy.pdf` | NJ | context for the `nj_distance` supporting notebooks |

**Absent but previously listed here:** `SNJ_Jaffe_Kluger.pdf` and
`Leveraged_Matrix_Completion_With_Noise.pdf` are not in this directory and never were under
those names. The SNJ material the section below describes is real and load-bearing in the
code -- only the PDF is missing.

All PDFs except the STDR one are **gitignored** (third-party reference material; 10.9 MB of
them were tracked by mistake in commit `b53e2f5` and untracked again). Several still carry
Finder's `" copy"` suffix and one has a doubled `.pdf.pdf` extension.

## SNJ paper — what's load-bearing in this repo

The paper formulates partition recovery on a **similarity** matrix `M_{ij} = exp(-d_{ij})` and uses **σ₂(M[A, ¬A])** — the second singular value of the cross-block — as the rank-1 deviation score (lower = cleaner clade). Both are already implemented:

- σ₂ scorer: `spectraltree/snj.py:9` (`sv2`) → consumed via `svd2` in `src/utils/metrics.py:326, 430`.
- α power: `spectraltree/snj.py:28-32` (`SpectralNeighborJoining.alpha`). Exposed as `SamplingConfig.distance_alpha` for the distance sub-sampling path, since `M^α = exp(-α D)`.

The similarity↔distance bridge `D = -log(M)` (paper Eq. defining paralinear distance) is `spectraltree.similarity2distance` (`spectraltree/similarities.py:7`); `paralinear_distance(obs)` at `spectraltree/similarities.py:34` is what the distance-mode builder uses.
