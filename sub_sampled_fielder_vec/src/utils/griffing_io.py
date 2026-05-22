"""Read Griffing-sweep results.

Layout:
  * ``griffing_meta.json`` (always)        — per-p means/stds, cfg snapshot
  * ``griffing_arrays.npz`` (optional)     — ``v_ref`` (m,) and
                                              ``v_hat_stack`` (n_p, n_reps, m)
                                              when the runner was called with
                                              ``save_eigvecs=True``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import numpy as np


def load_griffing_results(run_dir) -> Dict:
    """Return the meta dict, optionally with eigenvector arrays attached.

    If ``griffing_arrays.npz`` exists alongside the meta JSON, the returned
    dict gains ``v_ref`` and ``v_hat_stack`` keys (np.ndarray). Otherwise
    those keys are absent — the recovery / spec-norm panels use only the
    aggregated per-p means in ``meta['per_p']``.
    """
    p = Path(run_dir)
    meta_path = p / "griffing_meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(
            f"No Griffing results found in {run_dir} (expected griffing_meta.json)"
        )
    with open(meta_path) as f:
        meta = json.load(f)

    npz_path = p / "griffing_arrays.npz"
    if npz_path.exists():
        with np.load(npz_path) as arrs:
            meta["v_ref"] = arrs["v_ref"]
            meta["v_hat_stack"] = arrs["v_hat_stack"]
    return meta
