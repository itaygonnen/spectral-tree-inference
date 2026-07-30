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


def load_sweep(path: Path, cfg: BenchmarkConfig) -> Optional[Dict[str, object]]:
    """Cached curves for one tree, or ``None`` on a miss / settings mismatch."""
    if not path.exists():
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if dict(z["meta"].item()) != cfg.sweep_meta():
            return None
        return {k: z[k] for k in ("G", "L", "Lsym")} | {"rT": float(z["rT"])}
    except Exception:  # noqa: BLE001 — a truncated npz from a kill is just a miss
        return None


def save_sweep(out_dir: Path, res: dict, cfg: BenchmarkConfig) -> None:
    np.savez(sweep_path(out_dir, res["tree"]),
             p_values=np.asarray(cfg.p_values, float),
             G=res["G"], L=res["L"], Lsym=res["Lsym"], rT=res["rT"],
             meta=cfg.sweep_meta())
