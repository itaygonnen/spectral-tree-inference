"""Simulated trees as run sources.

The real-data side answers "which directories hold alignments" and hands the runner a
:class:`~src.runners.experiment_run.Source`. This is the same answer for simulated
trees: a model, a size and a sequence length define a population, ``build_ids`` names
individuals in it, and ``GeneratedLoader`` turns a name into matrices. Everything after
that -- screen, gate, sweep, export -- is the code the real data uses, unchanged.

One source per (model, size) cell, mirroring one source per real dataset: the cost model
and the per-row CSV output are both keyed on the taxon count, so mixing two sizes into
one source would make both meaningless.
"""
from __future__ import annotations

from typing import List, Sequence

from .benchmark_loaders import GeneratedLoader
from .experiment_run import Source


def source_name(model: str, n_taxa: int, seq_len: int, mutation_rate: float,
                etas: Sequence[int] | None = None) -> str:
    """The cache slug and CSV label for one simulated cell.

    Everything that changes the trees is in the name, because the name is what the
    resumable screen and sweep caches are keyed on -- two runs that differ in mutation
    rate must not inherit each other's verdicts.
    """
    name = f"{model} n{n_taxa} L{seq_len} mu{mutation_rate:g}"
    if etas:
        name += " eta" + "-".join(str(e) for e in etas)
    return name


def sources_from_plan(plan) -> List[Source]:
    """One :class:`Source` per (model, size) cell of a
    :class:`src.utils.generated_prompts.GeneratedPlan`."""
    from ..utils.tree_ids import parse_spec

    loader = GeneratedLoader(
        seq_len=plan.seq_len,
        n_taxa=plan.n_values[0] if len(plan.n_values) == 1 else None,
        mutation_rate=plan.mutation_rate,
        params=plan.params,
    )

    # group the plan's ids by the cell they belong to, so each source gets exactly the
    # trees that exist for it (the eta pool may be short in some cells)
    by_cell: dict = {}
    for tid in plan.tree_ids:
        try:
            spec = parse_spec(tid)
            key = (spec.model, spec.n_taxa or (plan.n_values[0]
                                               if len(plan.n_values) == 1 else 0))
        except ValueError:
            key = (plan.models[0], plan.n_values[0])
        by_cell.setdefault(key, []).append(tid)

    out: List[Source] = []
    for (model, n_taxa), ids in sorted(by_cell.items()):
        out.append(Source(
            name=source_name(model, n_taxa, plan.seq_len, plan.mutation_rate,
                             plan.etas),
            loader=loader,
            ids=sorted(ids),
            m=int(n_taxa),
        ))
    return out
