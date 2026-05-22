"""Read first-layer Fiedler-recovery sweep results (``fiedler_meta.json``)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict


def load_fiedler_results(run_dir) -> Dict:
    p = Path(run_dir)
    meta_path = p / "fiedler_meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(
            f"No Fiedler results found in {run_dir} (expected fiedler_meta.json)"
        )
    with open(meta_path) as f:
        return json.load(f)
