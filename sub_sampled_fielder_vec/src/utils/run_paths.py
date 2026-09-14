"""Where a run's results and its resumable state live.

Two roots under ``<repo>/results/real_data/``, and the distinction matters:

    runs/<timestamp>-<name>/   ONE directory per run, holding every source that run
                               covered -- config, CSVs, plot, log. This is what you
                               download, mail or plot from, and it never changes once
                               the run finishes.
    _cache/<source>/           machine state so a killed run resumes: the screen rows
                               and one .npz per swept tree, per source, reused across
                               runs.

``$STR_RESULTS_DIR`` overrides the parent of both. The paths are exactly what the first
cluster runs produced, and must stay that way: changing one orphans a cache that cost
30 min per tree to fill.
"""
from __future__ import annotations

import os
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]      # spectral-tree-inference


def results_root() -> Path:
    env = os.environ.get("STR_RESULTS_DIR")
    return Path(env).expanduser() if env else _REPO / "results" / "real_data"


def slug(name: str) -> str:
    return "_".join(str(name).lower().split())


def cache_dir(source: str) -> Path:
    return results_root() / "_cache" / slug(source)


def screen_cache_path(source: str) -> Path:
    """Resumable screen state for one source: one row per tree."""
    return cache_dir(source) / "screen.npz"


def sweep_cache_dir(source: str) -> Path:
    """Resumable sweep state for one source: one .npz per tree."""
    return cache_dir(source) / "sweep"


def new_run_dir(name: str = "") -> Path:
    """``runs/<timestamp>[-<name>]/`` -- created empty, then filled by the run."""
    from datetime import datetime
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    d = results_root() / "runs" / (f"{stamp}-{slug(name)}" if name else stamp)
    d.mkdir(parents=True, exist_ok=True)
    return d
