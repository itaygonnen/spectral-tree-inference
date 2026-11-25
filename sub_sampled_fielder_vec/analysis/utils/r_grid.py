"""The iso-``r`` grid: which ``(eta, m)`` cells are expected to recover alike.

``eq:main_rate`` makes ``p`` sufficient once

.. math:: p > 8(1+\\eta)^3\\beta_0^2 \\log m / [m(\\rho-\\eta S_{out}^{max})^2],

so the dimensionless part that carries both knobs is

.. math:: r = (1+\\eta)^a \\log m / m,  a = 3.

Cells sharing an ``r`` should share a recovery threshold, which is what makes a
larger ``m`` able to pay for a larger ``eta``. Reading that claim off measured
``p*`` needs cells matched in ``r`` and spread in ``m``, and that is what this
module lays out -- for the synthetic flat CBM, where ``eta`` is ours to choose,
and for simulated trees, where it is not.

The exponent is a parameter, not a constant, because it is the thing under test:
the grid under ``a=3`` and the grid under ``a=1`` put a cell's partner at
different sizes, so a band that is flat under one and sloped under the other
measures ``a``.
"""

from __future__ import annotations

from typing import Dict, List, NamedTuple, Sequence

import numpy as np
from scipy.optimize import brentq

#: Default exponent -- the one eq:main_rate states.
EXPONENT: float = 3.0

#: Default targets, spaced by roughly sqrt(2) and capped at R_MAX so that r stays
#: readable as a rate. The ladder is deliberately uneven at the ends: a target
#: below ~0.1 reaches only one or two cells before m leaves the affordable range,
#: while 1.0 sits exactly on the boundary where the bound stops being a
#: probability and carries the widest eta span of any target (2..8).
R_TARGETS: tuple = (0.10, 0.20, 0.283, 0.40, 0.80, 1.00)

#: Largest r worth a cell: above 1 the bound exceeds a probability and says
#: nothing (cor:infeasible), so a cell there tests nothing.
R_MAX: float = 1.0

#: Bounds on m for a synthetic cell. The floor is where a Fiedler vector of a
#: two-block matrix stops being a meaningful object; the ceiling is the measured
#: ~3 h/cell at m~11k (docs/RUNBOOK.md).
M_MIN: int = 90
M_MAX: int = 12288


def r_of(eta: np.ndarray | float, m: np.ndarray | float,
         exponent: float = EXPONENT) -> np.ndarray | float:
    """``r = (1+eta)^exponent * log(m) / m``."""
    return ((1.0 + np.asarray(eta, float)) ** exponent) * np.log(m) / np.asarray(m, float)


def eta_for(r: float, m: float, exponent: float = EXPONENT) -> float:
    """The ``eta`` a tree of ``m`` leaves needs to sit at ``r``.

    The inverse of :func:`r_of` in ``eta``, closed form. Returns a value below 1
    when ``r`` is too small for that ``m`` -- not an error, but a cell that no
    tree can fill, since ``eta >= 1`` by definition.
    """
    return float((r * m / np.log(m)) ** (1.0 / exponent) - 1.0)


def m_for(r: float, eta: float, exponent: float = EXPONENT,
          m_hi: float = 1e9) -> float:
    """The ``m`` at which a given ``eta`` sits at ``r``.

    Solved numerically: ``log(m)/m`` is not invertible in closed form. Raises
    ``ValueError`` when no ``m`` works, which happens whenever ``r`` exceeds
    ``(1+eta)^exponent / e`` -- ``log(m)/m`` peaks at ``1/e``, so a large enough
    ``r`` is out of reach at *every* size, not merely an expensive one.
    """
    f = lambda mm: r_of(eta, mm, exponent) - r
    lo = 3.0  # log(m)/m is maximal near here, so the peak is inside the bracket
    if f(lo) * f(m_hi) > 0:
        raise ValueError(
            f"no m in [{lo:g}, {m_hi:g}] puts eta={eta:g} at r={r:g} "
            f"(ceiling is r={(1.0 + eta) ** exponent / np.e:.4g})"
        )
    return float(brentq(f, lo, m_hi))


def eta_window(eta: float, r_tol: float, exponent: float = EXPONENT) -> float:
    """Half-width in ``eta`` of a ``+/- r_tol`` relative window in ``r``.

    ``r`` is a power of ``(1+eta)``, so a relative tolerance in ``r`` is the
    exponent times the relative tolerance in ``(1+eta)``: ``dr/r = a*d(1+eta)/(1+eta)``.
    The window therefore *narrows* by a factor of ``a`` relative to the naive
    reading, which is why an r-target is harder to hit than an eta-target of the
    same nominal tolerance.
    """
    return (1.0 + eta) * r_tol / exponent


class SyntheticCell(NamedTuple):
    """One flat-CBM cell: ``eta`` exact, ``m`` snapped to a multiple of ``1+eta``."""

    r_target: float
    eta: int
    m: int
    r_actual: float

    @property
    def n1(self) -> int:
        """Smaller block size. Exact, because ``m`` is a multiple of ``1+eta``."""
        return self.m // (1 + self.eta)


def synthetic_grid(r_targets: Sequence[float] = R_TARGETS,
                   etas: Sequence[int] = tuple(range(1, 9)),
                   exponent: float = EXPONENT,
                   m_min: int = M_MIN,
                   m_max: int = M_MAX) -> List[SyntheticCell]:
    """Cells for the synthetic matrix, one per ``(r_target, eta)`` that fits.

    ``m`` is snapped to a multiple of ``(1+eta)`` so that ``n1 = m/(1+eta)`` is
    exact -- the flat CBM builds the split by integer division, so an ``m`` that
    does not divide gives a realised ``eta`` other than the nominal one, and the
    cell would then sit at an ``r`` other than its target.

    ``eta`` is kept an ``int`` on purpose: the matrix builder floors
    ``int(eta)``, and the cache key formats ``int`` and ``float`` differently, so
    an integer ``eta`` is the one value that means the same thing to both.
    """
    cells: List[SyntheticCell] = []
    for r in r_targets:
        for eta in etas:
            try:
                m_exact = m_for(r, eta, exponent)
            except ValueError:
                continue
            step = 1 + eta
            m = int(round(m_exact / step)) * step
            if not (m_min <= m <= m_max):
                continue
            cells.append(SyntheticCell(r, int(eta), m, float(r_of(eta, m, exponent))))
    return cells


class TreeCell(NamedTuple):
    """One simulated-tree cell: ``m`` chosen, ``eta`` demanded of the sampler."""

    r_target: float
    m: int
    eta: float
    eta_lo: float
    eta_hi: float

    @property
    def feasible(self) -> bool:
        """``eta >= 1``, the floor a bipartition cannot go below."""
        return self.eta >= 1.0


def tree_grid(r_targets: Sequence[float] = R_TARGETS,
              ms: Sequence[int] = (1000, 2000, 4000, 6000, 8000),
              r_tol: float = 0.30,
              exponent: float = EXPONENT) -> List[TreeCell]:
    """Cells for simulated trees: the ``eta`` window each ``(r, m)`` asks for.

    Unlike :func:`synthetic_grid` there is nothing to snap -- ``eta`` is measured
    off whatever tree the simulator produces, so the grid states a target and a
    window and the rejection sampler does the rest.
    """
    cells: List[TreeCell] = []
    for r in r_targets:
        for m in ms:
            eta = eta_for(r, m, exponent)
            w = eta_window(eta, r_tol, exponent)
            cells.append(TreeCell(r, int(m), eta, max(1.0, eta - w), eta + w))
    return cells


def acceptance(cells: Sequence[TreeCell], eta_sample: np.ndarray) -> Dict[tuple, float]:
    """Fraction of ``eta_sample`` inside each cell's window.

    ``eta_sample`` is an empirical draw of measured ``eta`` -- the screened
    population of a tree model. The result is the rejection sampler's expected
    acceptance rate per cell, so a fill can be priced before it is started
    rather than after.
    """
    e = np.asarray(eta_sample, float)
    e = e[np.isfinite(e)]
    return {
        (c.r_target, c.m): float(((e >= c.eta_lo) & (e <= c.eta_hi)).mean())
        for c in cells
    }


#: Targets for the L_sym grid: a geometric ladder at ratio 1.43, which is twice
#: the +-18% cell width, so the windows tile without overlapping.
#:
#: The ladder deliberately runs past the band every size can reach (0.030-0.126).
#: Restricted to that band, a tree was discarded for being unreachable at some
#: OTHER size -- 71 of 91 unbinned trees sat above it, nearly all at m=500 and
#: m=1000, where L_sym's eta already exceeds the top target. Extending the ladder
#: keeps them, at the price of cells that span three sizes instead of five. The
#: ends stop where a target reaches fewer than two sizes.
LSYM_TARGETS: tuple = (0.021, 0.030, 0.043, 0.062, 0.088, 0.126, 0.180, 0.258)


def bin_by_r(eta, m, targets: Sequence[float] = LSYM_TARGETS,
             tol: float = 0.18, exponent: float = EXPONENT) -> np.ndarray:
    """Nearest target in log-r, or NaN for a tree too far from every one.

    The window is relative and applied in log space, because ``r`` spans decades
    and a fixed width would be generous at the top of the range and impossible at
    the bottom. ``tol`` is that relative half-width: 0.18 accepts a tree whose r is
    within 18% of its target, which is about half the gap between neighbouring
    targets here, so the windows do not overlap.

    Returns the target value rather than an index, so a binned table carries the
    number it was matched to and not a position in a tuple that may change.
    """
    r = np.asarray(r_of(eta, m, exponent), float)
    t = np.asarray(sorted(targets), float)
    lr, lt = np.log(r)[:, None], np.log(t)[None, :]
    k = np.abs(lr - lt).argmin(axis=1)
    out = t[k]
    return np.where(np.abs(np.log(r / out)) <= np.log1p(tol), out, np.nan)


def occupancy(eta, m, sizes: Sequence[int],
              targets: Sequence[float] = LSYM_TARGETS,
              tol: float = 0.18) -> "dict":
    """``{(target, m): count}`` plus ``("unbinned", m)`` for what matched nothing.

    The unbinned row is the point of reporting this at all: it is what the cell
    cut costs, and a table that hides it reads as though the population were the
    grid.
    """
    eta = np.asarray(eta, float)
    m = np.asarray(m, float)
    lab = bin_by_r(eta, m, targets, tol)
    out = {}
    for size in sizes:
        at = m == size
        for t in sorted(targets):
            out[(t, int(size))] = int((at & (lab == t)).sum())
        out[("unbinned", int(size))] = int((at & np.isnan(lab)).sum())
    return out
