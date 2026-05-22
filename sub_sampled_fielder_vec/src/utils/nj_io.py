"""Read NJ-sweep results in either layout:

  * legacy JSON-only:  ``nj_results.json`` with per-pair Q arrays inline
  * split (preferred): ``nj_meta.json`` + ``nj_arrays.npz``
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import numpy as np


_EXTRA_KEYS = (
    "linf_mean", "linf_std", "linf_per_rep",
    "rel_specnorm_mean", "rel_specnorm_std", "rel_specnorm_per_rep",
    "rel_frob_mean", "rel_frob_std", "rel_frob_per_rep",
)


def _merge_extra_metrics(meta: Dict, run_dir: Path) -> None:
    """Merge ``nj_meta_extra.json`` (written by
    ``scripts/nj_recompute_normalized_metrics.py``) into ``meta['per_p']``.
    Matches by ``p`` value. No-op if the file is missing — keeps the loader
    backward-compatible with pre-diagnostic sweeps."""
    extra_path = run_dir / "nj_meta_extra.json"
    if not extra_path.exists():
        return
    with open(extra_path) as f:
        extra = json.load(f)
    extras_by_p: Dict[float, Dict] = {float(e["p"]): e for e in extra.get("per_p", [])}
    for entry in meta["per_p"]:
        ex = extras_by_p.get(float(entry["p"]))
        if ex is None:
            continue
        for k in _EXTRA_KEYS:
            if k in ex:
                entry[k] = ex[k]
    for k in ("D_norm_2", "D_norm_F", "D_norm_inf"):
        if k in extra:
            meta[k] = extra[k]


def load_nj_results(run_dir) -> Dict:
    """Return a dict with the same shape the runner used to write as JSON.

    ``per_p[i]`` gains ``q_adj`` and ``q_nonadj`` as ``np.ndarray``s when
    loaded from the split layout (faster than rebuilding Python lists);
    legacy JSON returns them as lists.

    If ``nj_meta_extra.json`` is present in the same directory, its per-p
    fields (``linf_mean/std``, ``rel_specnorm_mean/std``,
    ``rel_frob_mean/std``) are merged into the matching ``per_p`` entries.
    """
    p = Path(run_dir)
    meta_path = p / "nj_meta.json"
    npz_path = p / "nj_arrays.npz"
    if meta_path.exists() and npz_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
        with np.load(npz_path) as arrs:
            adj = arrs["q_adj"]
            nonadj = arrs["q_nonadj"]
        for idx, entry in enumerate(meta["per_p"]):
            entry["q_adj"] = adj[idx]
            entry["q_nonadj"] = nonadj[idx]
        _merge_extra_metrics(meta, p)
        return meta

    legacy_path = p / "nj_results.json"
    if legacy_path.exists():
        with open(legacy_path) as f:
            meta = json.load(f)
        _merge_extra_metrics(meta, p)
        return meta

    raise FileNotFoundError(
        f"No NJ results found in {run_dir} (expected nj_meta.json + "
        f"nj_arrays.npz, or legacy nj_results.json)"
    )
