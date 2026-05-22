"""Similarity / distance matrix construction and subsampling.

Supports two matrix kinds (selected via ``matrix_kind``):

* ``"similarity"`` — JC similarity matrix from spectraltree (default). Sub-sample
  similarity entries directly.
* ``"distance"``  — paralinear distance D. Sub-sample distance entries, then
  return ``S = exp(-α · D̂)`` so the downstream Fiedler / σ₂ pipeline keeps
  consuming a similarity-shaped matrix. The α parameter is from the SNJ paper
  (Jaffe & Kluger, ``docs/papers/SNJ_Jaffe_Kluger.pdf``; equivalent to
  ``M^α = exp(-α D)``). See also ``spectraltree/snj.py:28``.
"""
import hashlib
import numpy as np
import spectraltree

from ..utils.logging import log_info, suppress_warnings
from .sampling import get_sampler


class SimilarityMatrixBuilder:
    """
    Encapsulates matrix construction and subsampling for the sub-sampled Fiedler
    pipeline. Despite the historical name, supports both similarity and distance
    matrix kinds via the ``matrix_kind`` argument.
    """

    def __init__(self, method: str = "uniform", matrix_kind: str = "similarity",
                 distance_alpha: float = 1.0, **method_kwargs):
        """
        Initialize the matrix builder.

        Args:
            method: Sampling method ("uniform", "leveraged", "lds").
            matrix_kind: "similarity" (JC) or "distance" (paralinear). Controls
                which matrix is built and sub-sampled, and which self-value the
                sampler enforces on the diagonal (1.0 vs 0.0).
            distance_alpha: SNJ α used only when matrix_kind == "distance".
                Sub-sampled D̂ is post-transformed to ``exp(-α · D̂)`` so the
                downstream Fiedler / σ₂ machinery is unchanged.
            **method_kwargs: Method-specific parameters passed to sampler.
        """
        if matrix_kind not in ("similarity", "distance"):
            raise ValueError(f"matrix_kind must be 'similarity' or 'distance', got {matrix_kind!r}")
        self.matrix_kind = matrix_kind
        self.distance_alpha = float(distance_alpha)
        self._self_value = 1.0 if matrix_kind == "similarity" else 0.0
        self._cache = {}
        self.sampler = get_sampler(method, self_value=self._self_value, **method_kwargs)

    def build_full(self, observations: np.ndarray) -> np.ndarray:
        """
        Build full matrix (with caching). Returns similarity or distance
        depending on ``self.matrix_kind``.
        """
        obs_hash = self._get_observations_hash(observations)
        if obs_hash in self._cache:
            return self._cache[obs_hash]

        if self.matrix_kind == "similarity":
            log_info('cache', "Computing and caching full similarity matrix (JC)...")
            with suppress_warnings('similarity'):
                full = spectraltree.JC_similarity_matrix(observations)
        else:
            log_info('cache', "Computing and caching full distance matrix (paralinear)...")
            with suppress_warnings('similarity'):
                full = spectraltree.paralinear_distance(observations)
        self._cache[obs_hash] = full
        return full

    def build_subsampled(self, observations: np.ndarray, p: float, seed: int = None,
                        min_similarity: float = 0.0) -> np.ndarray:
        """
        Build a sub-sampled matrix in similarity-shape for downstream use.

        For ``matrix_kind == "similarity"`` this returns the sub-sampled
        similarity matrix directly (diag=1). For ``matrix_kind == "distance"``
        the distance matrix is sub-sampled (diag=0) and the result is then
        transformed via ``S = exp(-α · D̂)`` so callers always receive a
        similarity-shaped matrix with diag=1.
        """
        full = self.build_full(observations)
        subsampled = self.sampler.sample(full, p, seed)

        if self.matrix_kind == "distance":
            # Densify sparse outputs (LDS) so exp() doesn't fill them dense
            # with exp(0)=1 implicitly.
            if not isinstance(subsampled, np.ndarray):
                subsampled = subsampled.toarray()
            # CRITICAL: the sampler returns D̂ with sampled entries = D_ij/p
            # (IPW-debiased) and unsampled entries = 0 (placeholder for missing).
            # exp(-α·0) = 1 would invert the semantics, treating "missing"
            # as "maximally similar". So we must mask unsampled entries back
            # to 0 *after* the transform. The sampling mask is inferred from
            # D̂ ≠ 0 (true paralinear distance = 0 only for identical
            # sequences, which is vanishingly rare for non-trivial L; if it
            # does happen the entry gets zeroed, which is the conservative
            # choice).
            sampled_mask = subsampled != 0.0  # off-diagonal sampled entries
            similarity = np.exp(-self.distance_alpha * subsampled)
            similarity = np.where(sampled_mask, similarity, 0.0)
            np.fill_diagonal(similarity, 1.0)
            return similarity

        return subsampled
    

    def _get_observations_hash(self, observations: np.ndarray) -> str:
        """
        Generate a fast hash key for the observations to use as cache key.
        
        Uses array metadata + corner values + checksum for speed.
        Much faster than MD5 for large arrays (no full copy needed).
        
        Args:
            observations: Observation matrix
            
        Returns:
            Hash key string
        """
        # Fast hash using shape + dtype + corner values + checksum
        # This avoids creating a full copy with tobytes() for large arrays
        shape_tuple = observations.shape
        dtype_str = observations.dtype.str
        
        # Get corner values (if array is non-empty)
        if observations.size > 0:
            corner_0_0 = int(observations.flat[0]) if observations.size > 0 else 0
            corner_n_n = int(observations.flat[-1]) if observations.size > 0 else 0
            # Compute checksum (mod to prevent overflow)
            checksum = int(observations.sum()) % (2**31)
        else:
            corner_0_0 = 0
            corner_n_n = 0
            checksum = 0
        
        # Create hash tuple
        hash_tuple = (shape_tuple, dtype_str, corner_0_0, corner_n_n, checksum)
        
        # Use Python's built-in hash (fast and sufficient for cache keys)
        return str(hash(hash_tuple))
    
    def clear_cache(self):
        """Clear the similarity matrix cache."""
        self._cache.clear()
        log_info('cache', "Similarity matrix cache cleared")
    
    @property
    def cache_size(self) -> int:
        """Return the number of cached matrices."""
        return len(self._cache)

