"""The operators a screen or a sweep can read a bipartition off, in one table.

Three ways to turn one tree's matrices into a two-clan split:

    S      Fiedler vector of  L(S) = Deg(S) - S           (the similarity route)
    Lsym   Fiedler vector of  L_sym = I - Dg^-1/2 S Dg^-1/2
    B      leading-|lambda| eigenvector of  B = H D H     (the distance route)

and two ways to threshold a vector into a split: k-means (k=2) on its entries, or
its sign. The choice is made per tree, not per run: k-means alone routinely isolates
a single taxon on real data (1/999, eta=999) -- a real pendant edge, so a validity
gate passes it, but a split no sub-sample can recover. Both candidates are scored
and the more balanced one wins, with both etas recorded so the choice is auditable.

``B`` is the exception: there the sign pattern *is* the partition (that is what the
distance route claims), so it carries one rule and its ``eta_B`` stays comparable
with every screen already on disk.

This table is the single definition shared by the screen (``operator_screen``) and
the sweep (``operator_sweep``); adding a fourth operator means adding a row here.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Sequence, Tuple

import numpy as np

MIN_SPLIT = 5          # matches the STDR/test convention used throughout the repo


@dataclass(frozen=True)
class Operator:
    """One operator: how to build its reference vector, and how it may be cut."""

    key: str            # screen-row key and cache column: "S", "Lsym", "B"
    arm: str            # sweep-column suffix: "L", "Lsym", "B"
    label: str          # what a figure legend says
    kind: str           # "fiedler" | "griffing"
    laplacian: str = "unnormalized"      # fiedler only
    rules: Tuple[str, ...] = ("kmeans", "sign")

    def vector(self, S: np.ndarray, D: np.ndarray) -> np.ndarray:
        """The full-matrix reference vector this operator reads its split off."""
        from src.core.utils import (compute_fiedler_from_laplacian,
                                    compute_laplacian, compute_normalized_laplacian)
        from src.utils.griffing import griffing_leading_eigvec
        if self.kind == "griffing":
            # lm_k1 is ARPACK for the single eigenpair B needs: same vector as the
            # dense default at 0.39 s instead of 17.2 s per solve at m=6000.
            return griffing_leading_eigvec(D, solver="lm_k1")
        lap = (compute_normalized_laplacian(S) if self.laplacian == "normalized"
               else compute_laplacian(S))
        return compute_fiedler_from_laplacian(lap)


# Order matters: it is the column order in screening.csv and the legend order.
OPERATORS: Dict[str, Operator] = {
    "S":    Operator("S", "L", r"$L(S)$ Fiedler", "fiedler", "unnormalized"),
    "Lsym": Operator("Lsym", "Lsym", r"$L_{sym}$ Fiedler", "fiedler", "normalized"),
    "B":    Operator("B", "B", r"$B = H\mathcal{D}H$", "griffing", rules=("sign",)),
}
ALL_OPERATORS: Tuple[str, ...] = tuple(OPERATORS)
ARM_OF = {k: op.arm for k, op in OPERATORS.items()}
KEY_OF_ARM = {op.arm: k for k, op in OPERATORS.items()}


def resolve(operators: Sequence[str] | None) -> Tuple[str, ...]:
    """Validate an operator selection, defaulting to all of them."""
    if not operators:
        return ALL_OPERATORS
    out = tuple(operators)
    unknown = [o for o in out if o not in OPERATORS]
    if unknown:
        raise ValueError(f"unknown operator(s) {unknown}; expected {ALL_OPERATORS}")
    return out


def cut(vector: np.ndarray, rule: str, min_split: int = MIN_SPLIT) -> np.ndarray:
    """Threshold a vector into a boolean bipartition under one rule."""
    if rule == "sign":
        return np.asarray(vector >= 0).astype(bool)
    if rule == "kmeans":
        from src.runners.p_sweep_inner import _kmeans_bipartition
        return np.asarray(_kmeans_bipartition(vector, min_split)).astype(bool)
    raise ValueError(f"unknown cut rule {rule!r}")


def best_cut(op: Operator, vector: np.ndarray, min_split: int = MIN_SPLIT):
    """``(rule, partition, {rule: eta})`` -- the most balanced of the operator's rules.

    Ties go to k-means, which was the historical default, so a tree whose two rules
    agree keeps the verdict it already has on disk.
    """
    from src.utils.partition_metrics import eta as _eta
    parts = {r: cut(vector, r, min_split) for r in op.rules}
    etas = {r: _eta(p) for r, p in parts.items()}
    rule = min(parts, key=lambda r: (etas[r], r != "kmeans"))
    return rule, parts[rule], etas
