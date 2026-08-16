"""``results.json`` for the operator-comparison benchmark — the per-p metric table.

``summary.json`` answers "how did each operator do overall" (one p* per operator);
this answers "what was every metric at every p, for every size". One row per
(cell, operator, p), where a *cell* is one (model, n, eta) combination of the run.

The shape is the one the sweep experiments already use — ``{"columns": [...],
"rows": [...]}`` under ``n{n}_L{L}/`` — so ``src/utils/merge_results.py`` reads a
benchmark run with no changes, and so does anything else built against that
layout. The top-level ``results.json`` is the same table with every size in it,
for when you just want one file.

Metric columns are ``<metric>_<stat>`` over the trees in that cell: nmi, ari,
agreement, sign, dot x mean, median, std. Only metrics actually present in the
sweeps appear, so a run that swept just Griffing does not carry empty L columns.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from .benchmark_cache import SWEEP_METRICS, BenchmarkConfig, curve_key
from .benchmark_plots import P_STAR_THRESHOLD, p_star_per_tree

# Identity columns first, in this order, then the metric columns.
_ID_COLUMNS = ["p", "operator", "cell", "model", "num_taxa", "eta", "n_trees"]

__all__ = ["write_results_json", "build_rows"]


def _split_cell(cell: Optional[str]) -> Dict[str, object]:
    """``'kingman|n1000|eta05'`` -> model/eta fields; anything else -> just the label."""
    out: Dict[str, object] = {"cell": cell or "all", "model": None, "eta": None}
    if not cell:
        return out
    parts = cell.split("|")
    out["model"] = parts[0]
    for part in parts[1:]:
        if part.startswith("eta"):
            try:
                out["eta"] = int(part[3:])
            except ValueError:
                pass
    return out


def build_rows(cohorts: Dict[str, Sequence[str]], results: Dict[str, dict],
               cfg: BenchmarkConfig, cell_of: Dict[str, Optional[str]],
               n_taxa_of: Dict[str, Optional[int]]) -> List[dict]:
    """One row per (cell, operator, p): mean/median/std of each metric over trees."""
    P = np.asarray(cfg.p_values, float)
    rows: List[dict] = []

    for method in cfg.methods():
        trees = [t for t in cohorts.get(method, ())
                 if method in results.get(t, {})]
        by_cell: Dict[Optional[str], List[str]] = {}
        for t in trees:
            by_cell.setdefault(cell_of.get(t), []).append(t)

        for cell, cell_trees in sorted(by_cell.items(), key=lambda kv: str(kv[0])):
            # p* is a property of the (cell, operator) pair, not of one p; it is
            # repeated down the rows the way the sweep files repeat their
            # configuration columns.
            nmi = np.vstack([np.asarray(results[t][method], float)
                             for t in cell_trees])
            ps = p_star_per_tree(nmi, P)
            finite = np.isfinite(ps) & (ps > 0)
            p_star = float(np.median(ps[finite])) if finite.any() else None

            ns = {n_taxa_of.get(t) for t in cell_trees} - {None}
            stacks = {}
            for metric in SWEEP_METRICS:
                key = curve_key(method, metric)
                vals = [np.asarray(results[t][key], float)
                        for t in cell_trees if key in results[t]]
                if len(vals) == len(cell_trees):
                    stacks[metric] = np.vstack(vals)

            for i, p in enumerate(cfg.p_values):
                row: Dict[str, object] = {
                    "p": float(p), "operator": method,
                    "num_taxa": (int(next(iter(ns))) if len(ns) == 1 else None),
                    "n_trees": len(cell_trees),
                    "p_star_median": p_star,
                    "p_star_threshold": P_STAR_THRESHOLD,
                    **_split_cell(cell),
                }
                for metric, arr in stacks.items():
                    col = arr[:, i]
                    row[f"{metric}_mean"] = float(np.nanmean(col))
                    row[f"{metric}_median"] = float(np.nanmedian(col))
                    row[f"{metric}_std"] = float(np.nanstd(col))
                rows.append(row)
    return rows


def _table(rows: List[dict]) -> dict:
    columns = [c for c in _ID_COLUMNS if any(c in r for r in rows)]
    columns += sorted({c for r in rows for c in r} - set(columns))
    return {"columns": columns, "rows": rows}


def write_results_json(out_dir: Path, cohorts: Dict[str, Sequence[str]],
                       results: Dict[str, dict], cfg: BenchmarkConfig,
                       cell_of: Dict[str, Optional[str]],
                       n_taxa_of: Dict[str, Optional[int]],
                       seq_len: Optional[int] = None) -> List[Path]:
    """Write the run-level table, plus one ``n{n}_L{L}/results.json`` per size.

    The per-size files are what ``merge_results.py`` globs for; they are only
    written when ``seq_len`` is known (the generated path knows it, a FASTA
    directory does not) and when the rows carry a taxon count.
    """
    rows = build_rows(cohorts, results, cfg, cell_of, n_taxa_of)
    if not rows:
        return []

    written = [Path(out_dir) / "results.json"]
    written[0].write_text(json.dumps(_table(rows), indent=2, default=str))

    if seq_len:
        by_n: Dict[int, List[dict]] = {}
        for r in rows:
            if r.get("num_taxa"):
                by_n.setdefault(int(r["num_taxa"]), []).append(r)
        for n, sub in sorted(by_n.items()):
            d = Path(out_dir) / f"n{n}_L{int(seq_len)}"
            d.mkdir(parents=True, exist_ok=True)
            path = d / "results.json"
            path.write_text(json.dumps(_table(sub), indent=2, default=str))
            written.append(path)
    return written
