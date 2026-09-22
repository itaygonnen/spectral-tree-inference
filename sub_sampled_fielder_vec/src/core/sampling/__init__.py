"""Sampling package registry.

Uniform sampling is the only method. The leveraged/IALM and LDS samplers were
removed on 2026-09-17: nothing in the pipeline had selected them since the
screen/gate/p-sweep flow landed, and every result the manuscript rests on is
uniform. Their history is in git (``src/core/sampling/leveraged/``).
"""
from .base import BaseSampler


def get_sampler(method: str = "uniform", **kwargs) -> BaseSampler:
    """
    Get a sampler instance based on method name.

    Args:
        method: Sampling method name. Only ``"uniform"`` is supported.
        **kwargs: Parameters passed to the sampler constructor.

    Returns:
        BaseSampler instance

    Raises:
        ValueError: If method is not "uniform".
    """
    if method != "uniform":
        raise ValueError(
            f"Unknown sampling method: {method!r}. Only 'uniform' is supported."
        )
    from .uniform.sampler import UniformSampler
    return UniformSampler(**kwargs)


__all__ = ["BaseSampler", "get_sampler"]
