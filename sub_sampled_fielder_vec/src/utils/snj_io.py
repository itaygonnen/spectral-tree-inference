"""Read SNJ-sweep results in either layout:

  * legacy JSON-only:  ``snj_results.json`` with per-pair σ₂ arrays inline
  * split (preferred): ``snj_meta.json`` + ``snj_arrays.npz``
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import numpy as np


def load_snj_results(run_dir) -> Dict:
    """Return a dict with the same shape the runner used to write as JSON.

    ``per_p[i]`` gains ``sigma2_adj`` and ``sigma2_nonadj`` as ``np.ndarray``s
    when loaded from the split layout (faster than rebuilding Python lists);
    legacy JSON returns them as lists.
    """
    p = Path(run_dir)
    meta_path = p / "snj_meta.json"
    npz_path = p / "snj_arrays.npz"
    if meta_path.exists() and npz_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)
        with np.load(npz_path) as arrs:
            adj = arrs["sigma2_adj"]
            nonadj = arrs["sigma2_nonadj"]
        for idx, entry in enumerate(meta["per_p"]):
            entry["sigma2_adj"] = adj[idx]
            entry["sigma2_nonadj"] = nonadj[idx]
        return meta

    legacy_path = p / "snj_results.json"
    if legacy_path.exists():
        with open(legacy_path) as f:
            return json.load(f)

    raise FileNotFoundError(
        f"No SNJ results found in {run_dir} (expected snj_meta.json + "
        f"snj_arrays.npz, or legacy snj_results.json)"
    )
