"""Topology-agnostic spectral operations on a similarity matrix.

`compute_fiedler_of_S` builds the unnormalised Laplacian L = D - S and returns
the eigenvector of its second-smallest eigenvalue (the Fiedler vector). For an
S with constant row sums (e.g. balanced binary), the eigenvectors of S and L
coincide with reversed ordering, so this is equivalent to the second-largest
eigenvector of S used in the thesis.

`subsample_S` applies a symmetric Bernoulli(p) mask and rescales kept entries
by 1/p, leaving the diagonal untouched (matrix-completion convention).

`uniform_mask` / `ipw_from_mask` / `nnm_from_mask` are the figure-11 trio:
they let two reconstruction strategies (IPW, NNM via IALM) consume the *same*
uniform sample Ω, so the head-to-head is paired on sampling noise.
"""
from __future__ import annotations

import numpy as np
import scipy.linalg

from src.core.fiedler_computer import FiedlerVectorComputer
from src.core.utils import compute_laplacian


def compute_fiedler_of_S(S: np.ndarray, sampling_prob: float | None = None) -> np.ndarray:
    """Return the Fiedler vector of L = D - S with the project's sign convention."""
    laplacian = compute_laplacian(S)
    return FiedlerVectorComputer().compute(laplacian, sampling_prob=sampling_prob)


def subsample_S(S: np.ndarray, p: float, seed: int) -> np.ndarray:
    """Symmetric Bernoulli(p) sub-sample of S with 1/p rescaling on kept entries.

    The diagonal of S is preserved as-is (Sii = 1 in the population model and is
    always observed). Off-diagonal entries are kept independently with
    probability p across the upper triangle and mirrored to keep S_hat
    symmetric, then rescaled by 1/p.
    """
    if not 0.0 < p <= 1.0:
        raise ValueError(f"p must be in (0, 1], got {p}")
    n = S.shape[0]
    rng = np.random.default_rng(seed)

    upper_mask = rng.random((n, n), dtype=np.float32) < p
    upper_mask = np.triu(upper_mask, k=1)
    mask = upper_mask | upper_mask.T

    S_hat = np.zeros_like(S)
    S_hat[mask] = S[mask] / np.asarray(p, dtype=S.dtype)
    np.fill_diagonal(S_hat, np.diag(S))
    return S_hat


def uniform_mask(n: int, p: float, rng: np.random.Generator) -> np.ndarray:
    """Symmetric Bernoulli(p) mask on the strict upper triangle.

    Returns a boolean (n, n) array Ω with Ω = Ωᵀ, diagonal False. Sharing this
    mask across reconstruction strategies pairs the comparison on sampling
    noise (figure 11).
    """
    if not 0.0 < p <= 1.0:
        raise ValueError(f"p must be in (0, 1], got {p}")
    upper = rng.random((n, n)) < p
    upper = np.triu(upper, k=1)
    Omega = upper | upper.T
    np.fill_diagonal(Omega, False)
    return Omega


def ipw_from_mask(S: np.ndarray, Omega: np.ndarray, p: float) -> np.ndarray:
    """IPW (debiased) estimator: S[Ω]/p on Ω, 0 elsewhere; diagonal preserved.

    This is `subsample_S`'s estimator but driven by a pre-built Ω so the same
    sample can be fed to NNM as well.
    """
    S_hat = np.zeros_like(S)
    S_hat[Omega] = S[Omega] / np.asarray(p, dtype=S.dtype)
    np.fill_diagonal(S_hat, np.diag(S))
    return S_hat


def nnm_from_mask(
    S: np.ndarray,
    Omega: np.ndarray,
    lambda_param: float | None = None,
    max_iter: int = 200,
    tol: float = 1e-5,
    include_diagonal: bool = True,
) -> np.ndarray:
    """Nuclear-norm-regularised completion of S on Ω, via Soft-Impute.

    Solves   min_L  ½‖P_Ω(M − L)‖_F²  +  λ‖L‖_*
    by alternating: fill unobserved entries with the current L, then
    soft-threshold the singular values by λ. This is the Mazumder–Hastie–
    Tibshirani (2010) algorithm — the standard NNM solver for matrix
    completion, and the right tool for the Candès & Recht / Candès & Plan
    comparison under **uniform** sampling. (We deliberately do not call the
    RPCA-style `ialm_solve` here: that solver targets a different objective
    and is numerically fragile on this data.)

    For uniform sampling we default to λ = ‖P_Ω(M)‖_F / √(p·n²) — a scale-aware
    setting close to the "spectral noise" floor at sampling rate p. Override
    explicitly if you want a sensitivity sweep.

    The diagonal of M is always known a priori (= self-similarity = 1.0) for
    our similarity matrices, so `include_diagonal=True` augments Ω with the
    diagonal before running. This matches what `ipw_from_mask` does (it
    preserves `np.diag(S)` directly).
    """
    n = S.shape[0]
    Omega_eff = Omega.copy()
    if include_diagonal:
        np.fill_diagonal(Omega_eff, True)

    p_eff = Omega_eff.sum() / (n * n)
    if lambda_param is None:
        lambda_param = np.linalg.norm(S * Omega_eff, ord='fro') / np.sqrt(max(p_eff, 1e-10) * n * n)

    M_obs = S * Omega_eff  # observed entries (zero on Ω^c)
    L = np.zeros_like(S)
    prev_norm = 0.0
    for _ in range(max_iter):
        Y = M_obs + (~Omega_eff) * L  # fill unobserved with current L
        try:
            U, sv, Vt = scipy.linalg.svd(Y, full_matrices=False, lapack_driver='gesvd')
        except (np.linalg.LinAlgError, scipy.linalg.LinAlgError):
            # gesvd failed too — extremely ill-conditioned, e.g. essentially the
            # zero matrix at p·n < 1. Caller is expected to flag this case via a
            # degenerate-sample guard; we return the current L as a best effort.
            break
        sv_thr = np.maximum(sv - lambda_param, 0.0)
        L = (U * sv_thr) @ Vt
        norm = np.linalg.norm(L, ord='fro')
        if abs(norm - prev_norm) <= tol * (prev_norm + 1e-12):
            break
        prev_norm = norm
    L = (L + L.T) / 2
    if include_diagonal:
        np.fill_diagonal(L, np.diag(S))
    return L
