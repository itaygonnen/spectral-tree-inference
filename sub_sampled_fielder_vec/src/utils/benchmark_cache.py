"""Resumable on-disk state for the operator-comparison benchmark.

Two caches, both written incrementally so an interrupted run resumes:

``screen_table.csv``   one row per screened tree — per-operator validity flag and
                       partition imbalance eta. Human-readable on purpose; it is
                       both the resume shard and the record of *what each tree was
                       found valid for* (a tree can pass for S and fail for D).
``sweeps/<id>.npz``    one file per swept tree — the three NMI curves plus r(T).

``config.json`` holds the settings both caches were produced under, and mixing
incompatible settings into one directory is refused rather than silently merged.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

# Operator keys, in the order they appear in the screen table.
OPERATORS = ("S", "D", "Lsym")

# The same three operators carry two names: a *screen* key (the column in
# screen_table.csv) and a *method* key (the curve in sweeps/<id>.npz and in the
# figures). The mapping is not the identity — screen "S" is the unnormalized
# Laplacian of S, whose curve is called "L".
SCREEN_OF_METHOD = {"L": "S", "Lsym": "Lsym", "G": "D"}
METHOD_OF_SCREEN = {v: k for k, v in SCREEN_OF_METHOD.items()}
# Menu order — what the launcher numbers 1, 2, 3.
METHOD_KEYS = ("L", "Lsym", "G")
METHOD_LABELS = {"L": "L + k-means", "Lsym": "L_sym + k-means",
                 "G": "B = HDH (Griffing-on-D)"}
# Per-p metrics carried through the sweep. "nmi" is the primary curve (the one
# the figures plot and the cohort logic keys on); the rest ride along for
# results.json. The values are the keys ``bootstrap_p_sweep_simple`` returns.
SWEEP_METRIC_SOURCES = {
    "nmi": "partition_nmi_M",
    "ari": "partition_ari_M",
    "agreement": "partition_agreement_M",
    "sign": "sign_agreement",
    "dot": "dot_product",
}
SWEEP_METRICS = tuple(SWEEP_METRIC_SOURCES)


def curve_key(method: str, metric: str = "nmi") -> str:
    """``('L','nmi') -> 'L'``; ``('L','ari') -> 'L_ari'``.

    NMI keeps the bare method name so every sweep npz written before the other
    metrics existed still loads.
    """
    return method if metric == "nmi" else f"{method}_{metric}"


SCREEN_FIELDS = (["tree", "n_taxa"]
                 + [f for op in OPERATORS for f in (f"valid_{op}", f"eta_{op}")])


@dataclass
class BenchmarkConfig:
    """Knobs the professor is asked about, plus the notebook's fixed constants."""

    p_values: List[float] = field(
        default_factory=lambda: list(np.logspace(-2, 0, 20)))
    bootstrap_reps: int = 10
    n_compare: int = 100
    num_gaps: int = 10
    min_split: int = 5
    # Which operators to screen, sweep and plot. Method keys, menu order. All
    # three is the historical behaviour; a subset saves a full bootstrap sweep
    # per tree per operator dropped.
    operators: List[str] = field(default_factory=lambda: list(METHOD_KEYS))

    def methods(self) -> List[str]:
        """Selected method keys, always in menu order."""
        return [k for k in METHOD_KEYS if k in self.operators]

    def screen_ops(self) -> List[str]:
        """Selected operators as screen-table keys."""
        return [SCREEN_OF_METHOD[k] for k in self.methods()]

    def screen_meta(self) -> dict:
        return {"num_gaps": self.num_gaps, "min_split": self.min_split}

    def sweep_meta(self) -> dict:
        return {"p_values": [round(float(p), 9) for p in self.p_values],
                "bootstrap_reps": self.bootstrap_reps,
                "num_gaps": self.num_gaps, "min_split": self.min_split}


def check_meta(out_dir: Path, cfg: BenchmarkConfig, extra: dict) -> None:
    """Refuse to mix incompatible settings into one output directory."""
    path = out_dir / "config.json"
    meta = {"screen": cfg.screen_meta(), "sweep": cfg.sweep_meta(), **extra}
    if path.exists():
        old = json.loads(path.read_text())
        for section in ("screen", "sweep"):
            if old.get(section) != meta[section]:
                raise ValueError(
                    f"{path} was written with different {section} settings "
                    f"({old.get(section)} != {meta[section]}). Use a fresh output "
                    "directory, or delete this one to recompute.")
        # Generated runs describe their data entirely in source.data_key (model,
        # sizes, mu, alignment length, eta bins). Resuming a directory with a
        # different one would mix incomparable trees under the same ids, and the
        # screen/sweep sections above cannot see it. Absent on either side (older
        # run dirs) means no opinion.
        old_key = (old.get("source") or {}).get("data_key")
        new_key = (meta.get("source") or {}).get("data_key")
        if old_key and new_key and old_key != new_key:
            raise ValueError(
                f"{path} was written from a different data source "
                f"({old_key} != {new_key}). Use a fresh output directory.")
        meta = {**old, **meta}
    path.write_text(json.dumps(meta, indent=2, default=str))


def update_meta(out_dir: Path, **fields) -> None:
    path = out_dir / "config.json"
    meta = json.loads(path.read_text()) if path.exists() else {}
    meta.update(fields)
    path.write_text(json.dumps(meta, indent=2, default=str))


# --- screen table ----------------------------------------------------------

def read_screen_table(path: Path) -> List[dict]:
    """Rows already on disk, with the CSV strings parsed back to typed values."""
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["n_taxa"] = int(r["n_taxa"]) if r.get("n_taxa") else 0
        for op in OPERATORS:
            v, e = r.get(f"valid_{op}", ""), r.get(f"eta_{op}", "")
            # "" means unknown (no ground-truth tree), which is not the same as False.
            r[f"valid_{op}"] = {"True": True, "False": False}.get(v)
            r[f"eta_{op}"] = float(e) if e else None
    return rows


def append_screen_row(path: Path, row: dict) -> None:
    new = not path.exists()
    with path.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SCREEN_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def write_screen_npz(out_dir: Path, rows: List[dict],
                     tree_ids: Sequence[str]) -> None:
    """Also emit ``run_screen``-format caches so the notebooks can read them."""
    for op in OPERATORS:
        recs = [{"tree": r["tree"], "eta": r[f"eta_{op}"] or float("nan"),
                 "valid": bool(r[f"valid_{op}"])}
                for r in rows if r.get(f"valid_{op}") is not None]
        if recs:
            np.savez(out_dir / f"screen_{op}.npz",
                     screen_ids=np.array(list(tree_ids), dtype=object),
                     rows=np.array(recs, dtype=object))


# --- per-tree sweeps -------------------------------------------------------

def sweep_path(out_dir: Path, tree_id: str) -> Path:
    return out_dir / "sweeps" / f"{tree_id}.npz"


def load_sweep(path: Path, cfg: BenchmarkConfig,
               keys: Optional[Sequence[str]] = None
               ) -> Optional[Dict[str, object]]:
    """Cached curves for one tree, or ``None`` on a miss / settings mismatch.

    ``keys`` are the method curves the caller needs (default: everything the
    config selects). A file holding only some of them is a *miss*, which is why
    the operator selection is deliberately absent from ``sweep_meta``: an older
    all-three file still satisfies any subset, and a subset file simply misses
    when a later run asks for more.
    """
    if not path.exists():
        return None
    want = list(keys) if keys is not None else cfg.methods()
    try:
        z = np.load(path, allow_pickle=True)
        if dict(z["meta"].item()) != cfg.sweep_meta():
            return None
        have = set(z.files)
        if not set(want) <= have:
            return None
        extra = [curve_key(m, x) for m in want for x in SWEEP_METRICS
                 if curve_key(m, x) in have]
        return {k: z[k] for k in set(want) | set(extra)} | {"rT": float(z["rT"])}
    except Exception:  # noqa: BLE001 — a truncated npz from a kill is just a miss
        return None


def save_sweep(out_dir: Path, res: dict, cfg: BenchmarkConfig) -> None:
    """Persist the curves this result carries, merging with any already on disk.

    A tree swept for one operator today and another tomorrow ends up with both in
    the same file rather than losing the first.
    """
    path = sweep_path(out_dir, res["tree"])
    all_keys = [curve_key(m, x) for m in METHOD_KEYS for x in SWEEP_METRICS]
    curves = {k: np.asarray(res[k], float) for k in all_keys if k in res}
    if path.exists():
        try:
            z = np.load(path, allow_pickle=True)
            if dict(z["meta"].item()) == cfg.sweep_meta():
                for k in all_keys:
                    if k not in curves and k in z.files:
                        curves[k] = z[k]
        except Exception:  # noqa: BLE001 — unreadable old file, just overwrite
            pass
    np.savez(path, p_values=np.asarray(cfg.p_values, float),
             rT=res["rT"], meta=cfg.sweep_meta(), **curves)
