# Reference papers

| File | Topic | Where it's used |
|---|---|---|
| `SNJ_Jaffe_Kluger.pdf` | Spectral Neighbor Joining (SNJ) for latent tree models | σ₂ rank-1 partition score (`spectraltree/snj.py:9-15`, used in `src/utils/metrics.py:326,430`); SNJ α parameter `M^α = exp(-α D)` exposed as `SamplingConfig.distance_alpha` in the distance-matrix sub-sampling path (`src/core/similarity_builder.py`). |
| `Leveraged_Matrix_Completion_With_Noise.pdf` | Leveraged matrix completion + IALM recovery | Backs the leveraged sampling method (`src/core/sampling/leveraged/`). |
| `Spectral_top-down_recovery_of_latent_tree_models.pdf` | STDR / top-down spectral recovery | Backs `spectraltree/spectral_tree_reconstruction.py:STDR`. |

## SNJ paper — what's load-bearing in this repo

The paper formulates partition recovery on a **similarity** matrix `M_{ij} = exp(-d_{ij})` and uses **σ₂(M[A, ¬A])** — the second singular value of the cross-block — as the rank-1 deviation score (lower = cleaner clade). Both are already implemented:

- σ₂ scorer: `spectraltree/snj.py:9` (`sv2`) → consumed via `svd2` in `src/utils/metrics.py:326, 430`.
- α power: `spectraltree/snj.py:28-32` (`SpectralNeighborJoining.alpha`). Exposed as `SamplingConfig.distance_alpha` for the distance sub-sampling path, since `M^α = exp(-α D)`.

The similarity↔distance bridge `D = -log(M)` (paper Eq. defining paralinear distance) is `spectraltree.similarity2distance` (`spectraltree/similarities.py:7`); `paralinear_distance(obs)` at `spectraltree/similarities.py:34` is what the distance-mode builder uses.
