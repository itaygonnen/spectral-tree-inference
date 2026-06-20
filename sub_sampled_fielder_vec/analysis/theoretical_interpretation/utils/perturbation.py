"""Subspace-perturbation features for the Davis–Kahan verification of Lemma 3.4.

The "two failure modes" claim is a statement about a *ratio*: the sampling
perturbation ``E = L(Ŝ) − L(S)`` (the **statistical** threat, numerator) divided
by the population spectral gap ``Δλ = λ₃ − λ₂`` (the **algebraic** threat,
denominator). The Davis–Kahan sin-Θ theorem (Davis & Kahan 1970; the
statistician-friendly variant of Yu, Wang & Samworth 2015) bounds the rotation
of the rank-2 invariant subspace by exactly that ratio:

    ‖sin Θ(Û, U)‖ ≤ ‖E‖ / Δλ.

These helpers expose every piece of that inequality so a notebook can plot the
observed left-hand side against the predicted right-hand side and *attribute* the
growth to each mode. Topology-agnostic: plain dense ``np.ndarray`` in, scalars out.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
import scipy.linalg

from .linalg_features import compute_top_eigenpairs, cross_clan_variance


def _laplacian(S: np.ndarray) -> np.ndarray:
    """Unnormalised Laplacian L = D − S (same convention as linalg_features)."""
    return np.diag(S.sum(axis=1)) - S


def rank2_subspace(S: np.ndarray) -> Tuple[np.ndarray, float]:
    """Return the rank-2 bottom subspace U = [v₀, v₁] of L(S) and the gap Δλ.

    ``U`` is the orthonormal basis spanned by the two smallest-eigenvalue
    eigenvectors of L = D − S (the constant vector plus the Fiedler vector); the
    gap ``Δλ = λ₃ − λ₂`` separates it from the rest of the spectrum and is the
    Davis–Kahan denominator.
    """
    eigvals, eigvecs = compute_top_eigenpairs(_laplacian(S), k=3)
    U = eigvecs[:, :2]
    gap = float(eigvals[2] - eigvals[1])
    return U, gap


def laplacian_perturbation_norm(S: np.ndarray, S_hat: np.ndarray) -> float:
    """‖E‖₂ = ‖L(Ŝ) − L(S)‖_op — the Davis–Kahan numerator (sampling noise).

    This is the operator (spectral) norm of the Laplacian perturbation induced by
    sub-sampling. It inflates with the imbalance η through the cross-clan variance
    load ρ̃(η) = (1+η)/2·S_out², so it carries the *statistical* failure mode.
    """
    E = _laplacian(S_hat) - _laplacian(S)
    # E is symmetric: the operator norm is the largest-magnitude eigenvalue.
    ev = scipy.linalg.eigvalsh(E)
    return float(max(abs(ev[0]), abs(ev[-1])))


def subspace_sin_theta(U_ref: np.ndarray, U_hat: np.ndarray) -> float:
    """sin of the largest principal angle between span(U_ref) and span(U_hat).

    This is the operator-norm ‖sin Θ(Û, U)‖ — the observed left-hand side of the
    Davis–Kahan inequality, in [0, 1]: 0 = subspaces coincide, 1 = a direction of
    one is orthogonal to the other (the split is lost). Inputs need not be
    orthonormal; they are re-orthonormalised via QR first.
    """
    Qr, _ = np.linalg.qr(np.asarray(U_ref, dtype=float))
    Qh, _ = np.linalg.qr(np.asarray(U_hat, dtype=float))
    angles = scipy.linalg.subspace_angles(Qr, Qh)  # descending order
    return float(np.sin(angles).max())


def davis_kahan_bound(E_norm: float, gap: float) -> float:
    """Davis–Kahan / Yu–Wang–Samworth upper bound ‖E‖ / Δλ on sin Θ.

    Returns ``inf`` once the gap is non-positive: there the rank-2 subspace is no
    longer isolated and the bound is vacuous — precisely the algebraic-death
    regime η ≥ η* = ρ/S_out.
    """
    if gap <= 0:
        return float("inf")
    return float(E_norm) / float(gap)


def two_mode_decomposition(eta: float, S_out: float, p: float, gap: float) -> dict:
    """Closed-form split of the sin-Θ bound into its two compounding modes.

    Returns the *statistical* load ``ρ̃(η)/p = (1+η)/2·S_out²/p`` (numerator-side,
    grows with η and shrinks with p), the *algebraic* term ``1/Δλ``
    (denominator-side, grows as the gap collapses), and their product — the
    predictor whose super-linear blow-up in η is the quantitative signature of
    *compounding* failure. These are proportionality indicators, not the exact
    ‖E‖; pair them with the empirical :func:`laplacian_perturbation_norm`.
    """
    statistical = cross_clan_variance(eta, S_out) / float(p)
    algebraic = float("inf") if gap <= 0 else 1.0 / float(gap)
    return {
        "statistical": statistical,
        "algebraic": algebraic,
        "product": statistical * algebraic,
    }
