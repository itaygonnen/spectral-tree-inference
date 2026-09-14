"""The operators a screen or a sweep can read a bipartition off, in one table.

Three ways to turn one tree's matrices into a two-clan split:

    S      Fiedler vector of  L(S) = Deg(S) - S           (the similarity route)
    Lsym   Fiedler vector of  L_sym = I - Dg^-1/2 S Dg^-1/2
    B      leading-|lambda| eigenvector of  B = H D H     (the distance route)

and three ways to threshold a vector into a split: k-means (k=2) on its entries, its
sign, or the sigma2 gap search (``partition_taxa``). Which one is used is a *policy*,
not a constant, because the two pipelines built on this table want different answers:

``"best"``     score every rule the operator allows and keep the most balanced split.
               This is the real-data default: k-means alone routinely isolates a single
               taxon on real data (1/999, eta=999) -- a real pendant edge, so a validity
               gate passes it, but a split no sub-sample can recover. Both candidates'
               etas are recorded, so the choice is auditable.
``"kmeans"``   force k-means. What the generated operator-comparison benchmark sweeps.
``"sigma2"``   force the gap search. What that benchmark screens with.

Forcing a rule is not a fallback -- ``bpart_sweep_cache`` records that the generated
benchmark's accounting is like-for-like with a published figure, so its rules are fixed
and this table has to be able to express them.

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
    rules: Tuple[str, ...] = ("kmeans", "sign")        # candidates under policy "best"
    forced: Tuple[str, ...] = ("kmeans", "sign", "sigma2")   # what a policy may force

    def resolve_rule(self, policy: str) -> str:
        """The rule this operator uses under ``policy``.

        An operator with a single rule ignores a forced policy -- there is nothing to
        choose. That is not a fallback but the definition: on ``B`` the sign pattern IS
        the partition, so "cut B by the sigma2 gap search" names nothing.
        """
        if policy in self.forced:
            return policy
        if len(self.rules) == 1:
            return self.rules[0]
        raise ValueError(f"operator {self.key!r} cannot be cut by {policy!r}; "
                         f"it allows {self.forced}")

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
    "B":    Operator("B", "B", r"$B = H\mathcal{D}H$", "griffing", rules=("sign",),
                     forced=("sign",)),
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


RULES = ("kmeans", "sign", "sigma2")


def cut(vector: np.ndarray, rule: str, min_split: int = MIN_SPLIT,
        matrix: np.ndarray | None = None, num_gaps: int = 10) -> np.ndarray:
    """Threshold a vector into a boolean bipartition under one rule.

    ``sigma2`` scores candidate thresholds by the second singular value of the
    cross-partition block, so it needs the matrix the vector came from.
    """
    if rule == "sign":
        return np.asarray(vector >= 0).astype(bool)
    if rule == "kmeans":
        from src.runners.p_sweep_inner import _kmeans_bipartition
        return np.asarray(_kmeans_bipartition(vector, min_split)).astype(bool)
    if rule == "sigma2":
        if matrix is None:
            raise ValueError("the sigma2 rule needs the matrix the vector came from")
        from src.utils.metrics import compute_reference_partition_and_quality
        part, _, _ = compute_reference_partition_and_quality(
            vector, matrix, num_gaps=num_gaps, min_split=min_split)
        return np.asarray(part).astype(bool)
    raise ValueError(f"unknown cut rule {rule!r}; expected one of {RULES}")


def choose_cut(op: Operator, vector: np.ndarray, min_split: int = MIN_SPLIT,
               matrix: np.ndarray | None = None, num_gaps: int = 10,
               policy: str = "best"):
    """``(rule, partition, {rule: eta})`` under a cut policy.

    ``policy="best"`` scores every rule the operator allows and keeps the most balanced
    split; ties go to k-means, which was the historical default, so a tree whose rules
    agree keeps the verdict it already has on disk. Any other value names a rule and
    forces it, and then the returned eta map holds that rule alone.
    """
    from src.utils.partition_metrics import eta as _eta
    if policy != "best":
        rule = op.resolve_rule(policy)
        part = cut(vector, rule, min_split, matrix, num_gaps)
        return rule, part, {rule: _eta(part)}
    parts = {r: cut(vector, r, min_split, matrix, num_gaps) for r in op.rules}
    etas = {r: _eta(p) for r, p in parts.items()}
    rule = min(parts, key=lambda r: (etas[r], r != "kmeans"))
    return rule, parts[rule], etas


def best_cut(op: Operator, vector: np.ndarray, min_split: int = MIN_SPLIT):
    """The ``policy="best"`` case, kept as its own name because it is the common one."""
    return choose_cut(op, vector, min_split)
