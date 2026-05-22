"""Base abstract class for matrix samplers."""
from abc import ABC, abstractmethod
import numpy as np


class BaseSampler(ABC):
    """
    Abstract base class for matrix subsampling methods.

    All samplers must implement the sample() method which takes a full matrix
    and returns a subsampled or recovered version of it. Subclasses store the
    expected diagonal value (1.0 for similarity, 0.0 for distance) on
    ``self.self_value`` and use it when enforcing the diagonal after sampling.
    """

    def __init__(self, self_value: float = 1.0, **kwargs):
        self.self_value = self_value

    @abstractmethod
    def sample(self, matrix: np.ndarray, p: float, seed: int = None, **kwargs) -> np.ndarray:
        """
        Sample or recover a matrix from the full matrix.

        Args:
            matrix: Full symmetric matrix (n x n) — similarity or distance.
            p: Sampling probability or budget parameter (0 < p <= 1)
            seed: Random seed for reproducibility
            **kwargs: Method-specific additional parameters

        Returns:
            Subsampled or recovered matrix (n x n) with diagonal set to self.self_value.
        """
        pass

