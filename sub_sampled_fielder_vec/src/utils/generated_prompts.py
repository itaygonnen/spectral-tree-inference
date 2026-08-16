"""Interactive setup for the benchmark's ``generated`` data source.

Imported by ``scripts/run_benchmark.py``.

Asks, in this order: which tree model(s), that model's properties, how many taxa
(a list is allowed), the alignment length and mutation rate, and finally whether
to draw from the eta pool and which bins. Every answer is Enter-for-default.

The total number of trees is *derived* from those answers — models x sizes x
(eta bins) x trees-per-cell — rather than asked for separately. That is
deliberate: one "how many trees" prompt applied to a flat list silently
truncates it, and since ids are laid out cell by cell, a truncated list can drop
a whole model or size class.

Eta pooling is served by the shared pool under ``cache/pool_sample``, the one
``scripts/build_eta_pool.py`` writes, so bins already on disk cost nothing and
anything built here is visible to the notebooks reading that pool.

"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .interactive_ui import (
    Colors, confirm, get_input, print_header, print_option, print_warning,
)
from .tree_ids import MODEL_OF_TAG, build_ids, make_tree_id

# tag -> one-line description. Tags are what appear in a tree id; the pool and
# ``src/models/tree_models.py`` call 'bd' 'birth_death' (hence MODEL_OF_TAG).
MODEL_MENU: List[Tuple[str, str]] = [
    ("kingman", "coalescent — random topology, realistic branch lengths"),
    ("bd", "birth-death — random topology"),
    ("lopsided", "caterpillar — fixed topology, only the alignment varies"),
    ("balanced_binary", "balanced — fixed topology, taxa must be a power of 2"),
]

# tag -> [(param, prompt, default)] — the knobs that model actually reads
MODEL_PROPS: Dict[str, List[Tuple[str, str, str]]] = {
    "kingman": [("pop_size", "population size", "1.0")],
    "bd": [("birth_rate", "birth rate", "0.5"),
           ("death_rate", "death rate", "0.0")],
    "lopsided": [("edge_length", "edge length", "1.0")],
    "balanced_binary": [("edge_length", "edge length", "1.0")],
}

# What every eta pool on disk was built with. Offered as the default so accepting
# defaults hits the existing cache instead of a fresh rejection-sampling run.
POOL_MU, POOL_SEQ_LEN, POOL_POP = 0.1, 10_000, 1.0

# Only these can be eta-pooled — it is what build_eta_pool.py --tree-model takes.
POOLABLE = ("kingman", "kingman_mean", "lopsided", "birth_death",
            "balanced_binary")


@dataclass
class GeneratedPlan:
    """Everything the run needs, plus the ids it intends to use."""

    models: List[str]
    n_values: List[int]
    seq_len: int
    mutation_rate: float
    params: Dict[str, float]
    etas: Optional[List[int]]
    per_cell: int
    max_attempts: int = 2000
    tree_ids: List[str] = field(default_factory=list)

    @property
    def pooled(self) -> bool:
        return bool(self.etas)

    @property
    def pop_size(self) -> float:
        return float(self.params.get("pop_size", POOL_POP))

    def pool_models(self) -> List[str]:
        return [MODEL_OF_TAG.get(t, t) for t in self.models]

    def data_key(self) -> dict:
        """The part of the source that changes the data itself."""
        return {"kind": "generated", "models": self.models,
                "n_values": self.n_values, "seq_len": self.seq_len,
                "mutation_rate": self.mutation_rate, "params": self.params,
                "etas": self.etas}

    def meta(self) -> dict:
        return {**self.data_key(), "data_key": self.data_key(),
                "per_cell": self.per_cell,
                "label": f"{'+'.join(self.models)} n={self.n_values}"
                         + (f" eta={self.etas}" if self.pooled else "")}


def _ask_models() -> List[str]:
    print_header("Step 2 — tree model")
    for i, (tag, desc) in enumerate(MODEL_MENU, 1):
        print_option(str(i), f"{tag:16s} {desc}")
    valid = {str(i): tag for i, (tag, _) in enumerate(MODEL_MENU, 1)}
    valid.update({tag: tag for tag, _ in MODEL_MENU})
    while True:
        raw = get_input("\nModels — comma-separated numbers or names",
                        default="1,2")
        picked = [valid.get(tok.strip().lower())
                  for tok in raw.split(",") if tok.strip()]
        if picked and all(picked):
            return list(dict.fromkeys(picked))  # dedupe, keep typed order
        print_warning(f"unknown model in {raw!r}; pick from "
                      f"{', '.join(t for t, _ in MODEL_MENU)}")


def _ask_props(models: Sequence[str]) -> Dict[str, float]:
    print_header("Step 3 — model properties")
    params: Dict[str, float] = {}
    for tag in models:
        for name, label, default in MODEL_PROPS.get(tag, []):
            if name not in params:
                params[name] = float(get_input(f"({tag}) {label}", default=default))
    if not params:
        print("  (nothing to set for the selected models)")
    return params


def _ask_sizes(models: Sequence[str]) -> List[int]:
    print_header("Step 4 — taxa")
    while True:
        raw = get_input("Number of taxa — comma-separated for several sizes",
                        default="1000")
        try:
            sizes = sorted({int(tok) for tok in raw.split(",") if tok.strip()})
        except ValueError:
            print_warning(f"not a list of integers: {raw!r}")
            continue
        if not sizes:
            continue
        bad = [n for n in sizes if n & (n - 1)] if "balanced_binary" in models else []
        if bad:
            print_warning(f"balanced_binary needs powers of two; {bad} are not")
            continue
        return sizes


def _ask_sequences() -> Tuple[int, float]:
    print_header("Step 5 — alignment")
    seq_len = int(get_input("Sequence length (sites per taxon)",
                            default=str(POOL_SEQ_LEN)))
    mu = float(get_input("JC mutation rate", default=str(POOL_MU)))
    return seq_len, mu


def _ask_eta(models: Sequence[str]) -> Tuple[Optional[List[int]], int, int]:
    print_header("Step 6 — eta pooling")
    print("  A tree's eta is the imbalance max(n1,n2)/min(n1,n2) of the split its\n"
          "  reference partition finds. Pooling keeps only trees whose eta lands\n"
          "  near a target (+/-1), drawn from the shared pool under\n"
          "  cache/pool_sample. Without it you get whatever imbalance the model\n"
          "  happens to produce.\n")
    unpoolable = [t for t in models if MODEL_OF_TAG.get(t, t) not in POOLABLE]
    if unpoolable:
        print_warning(f"cannot be eta-pooled: {', '.join(unpoolable)}")
        return None, int(get_input("Trees per (model, size) cell", default="300")), 0

    if not confirm("Use eta pooling?", default=False):
        return None, int(get_input("\nTrees per (model, size) cell",
                                   default="300")), 0

    while True:
        raw = get_input("Eta targets — comma-separated", default="1,5,10,15")
        try:
            etas = sorted({int(tok) for tok in raw.split(",") if tok.strip()})
            break
        except ValueError:
            print_warning(f"not a list of integers: {raw!r}")
    per_cell = int(get_input("Trees per (model, size, eta) cell", default="10"))
    max_attempts = int(get_input(
        "Max rejection attempts per size, when the pool is short", default="2000"))
    return etas, per_cell, max_attempts


def make_more_ids(plan: GeneratedPlan) -> Callable[[Dict[str, int]], List[str]]:
    """A supplier of replacement candidates, one draw at a time per short cell.

    The screen calls this when a cell ends up with too few *valid* trees. Pooled
    cells draw the next unused slots in their bin (bounded by what is on disk);
    plain cells just keep incrementing the index, which is unbounded because the
    tree is generated from the id.
    """
    from ..runners.eta_pool_bridge import available_samples

    issued: Dict[str, int] = {}

    def cell_key(model: str, n: int, eta: Optional[int]) -> str:
        parts = [model, f"n{n}"]
        if eta is not None:
            parts.append(f"eta{eta:02d}")
        return "|".join(parts)

    # cell -> (tag, n, eta, ceiling) where ceiling is None for unbounded supply
    supply: Dict[str, tuple] = {}
    kw = dict(seq_len=plan.seq_len, mu=plan.mutation_rate, pop_size=plan.pop_size)
    for tag, model in zip(plan.models, plan.pool_models()):
        for n in plan.n_values:
            for eta in (plan.etas or [None]):
                have = (len(available_samples(model, n, eta, **kw))
                        if eta is not None else None)
                supply[cell_key(tag, n, eta)] = (tag, n, eta, have)

    def more(short: Dict[str, int]) -> List[str]:
        out: List[str] = []
        for cell, need in short.items():
            spec = supply.get(cell)
            if spec is None:
                continue
            tag, n, eta, ceiling = spec
            start = issued.get(cell, plan.per_cell)
            stop = start + need if ceiling is None else min(start + need, ceiling)
            out += [make_tree_id(tag, i, n, eta) for i in range(start, stop)]
            issued[cell] = stop
        return out

    return more


def print_pool_table(plan: GeneratedPlan) -> Dict[Tuple[str, int, int], int]:
    """Show have/need per cell; returns the counts so callers can reuse them."""
    from ..runners.eta_pool_bridge import counts_by_cell

    models = plan.pool_models()
    counts = counts_by_cell(models, plan.n_values, plan.etas or [],
                            seq_len=plan.seq_len, mu=plan.mutation_rate,
                            pop_size=plan.pop_size)
    print(f"\n  pool availability (mu={plan.mutation_rate}, L={plan.seq_len}, "
          f"pop={plan.pop_size}):")
    for model in models:
        for n in plan.n_values:
            cells = "  ".join(f"eta{e:02d} {counts[(model, n, e)]}/{plan.per_cell}"
                              for e in (plan.etas or []))
            short = any(counts[(model, n, e)] < plan.per_cell
                        for e in (plan.etas or []))
            mark = (f"{Colors.YELLOW}build{Colors.RESET}" if short
                    else f"{Colors.GREEN}ready{Colors.RESET}")
            print(f"    {model:16s} n={n:<6d} {cells}   [{mark}]")
    if all(c == 0 for c in counts.values()):
        print_warning("nothing cached for these parameters — the pool on disk was "
                      f"built with mu={POOL_MU}, L={POOL_SEQ_LEN}, "
                      f"pop={POOL_POP}, kingman")
    return counts


def prompt_generated_plan() -> GeneratedPlan:
    """Run the generated-source dialogue and return the resolved plan."""
    models = _ask_models()
    params = _ask_props(models)
    n_values = _ask_sizes(models)
    seq_len, mu = _ask_sequences()
    etas, per_cell, max_attempts = _ask_eta(models)

    plan = GeneratedPlan(models=models, n_values=n_values, seq_len=seq_len,
                         mutation_rate=mu, params=params, etas=etas,
                         per_cell=per_cell, max_attempts=max_attempts)
    plan.tree_ids = build_ids(models=models, n_values=n_values, etas=etas,
                              per_cell=per_cell)

    cells = len(models) * len(n_values) * max(1, len(etas or []))
    print(f"\n  → {len(plan.tree_ids)} trees = {cells} cells x {per_cell} each")
    if plan.pooled:
        print_pool_table(plan)
    else:
        print(f"  {Colors.CYAN}note{Colors.RESET}: trees are deterministic in "
              "their id, so reruns rebuild the same data")
    return plan


def resolve_pool_ids(plan: GeneratedPlan) -> List[str]:
    """Top the pool up if needed, then return the ids that actually exist.

    Called after the run is confirmed, because building a bin is the expensive
    part. Cells still short after the build contribute what they have and are
    reported rather than silently dropped.
    """
    from ..runners.eta_pool_bridge import available_samples, ensure_cells

    models = plan.pool_models()
    kw = dict(seq_len=plan.seq_len, mu=plan.mutation_rate, pop_size=plan.pop_size)
    ensure_cells(models, plan.n_values, plan.etas or [], plan.per_cell,
                 max_attempts=plan.max_attempts, **kw)

    slots: Dict[Tuple[str, int, int], int] = {}
    short: List[str] = []
    for tag, model in zip(plan.models, models):
        for n in plan.n_values:
            for e in plan.etas or []:
                use = min(len(available_samples(model, n, e, **kw)), plan.per_cell)
                slots[(tag, n, e)] = use
                if use < plan.per_cell:
                    short.append(f"{model} n={n} eta{e:02d} {use}/{plan.per_cell}")
    if short:
        print_warning("cells still short after the build: " + "; ".join(short))

    # index-major, so any prefix of the list stays balanced across cells
    return [make_tree_id(tag, i, n, e)
            for i in range(plan.per_cell)
            for (tag, n, e), use in slots.items() if i < use]
