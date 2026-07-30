#!/usr/bin/env python3
"""Sync/verify docs/overleafs/v9/figures/ against the notebooks that produce them.

Default mode copies any figure with a known alternate source (currently just
`coherence_vs_eta_by_n.png`, whose generator notebook is unowned this round
and not re-run) into v9/figures/, then prints a manifest of every entry in
PAPER_FIGURES: (filename, md5, mtime, producing notebook).

--check is read-only: reports any PAPER_FIGURES entry missing from
v9/figures/, whose notebook was modified long after the figure was written
(STALE), or whose sibling figures from the same notebook were written far
apart in time (MIXED-RUN — the signal that would have caught the v8 HBM
bug where two panels of one figure came from two different runs, see
open-items/00-R1.md [R1/23]). Exits non-zero if any are found.

Usage:
    python scripts/sync_paper_figures.py            # copy + manifest
    python scripts/sync_paper_figures.py --check     # lint only, no writes
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis" / "theoretical_interpretation"
FIGURES_DIR = ROOT / "docs" / "overleafs" / "v9" / "figures"

# filename -> producing notebook, relative to ANALYSIS.
PAPER_FIGURES = {
    # C3 (this agent) — synthesized/cbm_theory
    "pstar_synth_nmi_eta1.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_nmi_eta2.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_nmi_eta4.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_nmi_eta8.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_agr_eta1.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_agr_eta2.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_agr_eta4.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_synth_agr_eta8.png": "synthesized/cbm_theory/nonbalanced_flat_cbm.ipynb",
    "pstar_balanced_nmi.png": "synthesized/cbm_theory/balanced_binary_threshold.ipynb",
    "S_by_alpha.png": "synthesized/cbm_theory/hbm_spectral_gap_verification.ipynb",
    "spectral_gap_bound.png": "synthesized/cbm_theory/hbm_spectral_gap_verification.ipynb",
    # unowned this round — nobody re-runs this notebook; see open-items/16-C3.md
    "coherence_vs_eta_by_n.png": "generated/eta_pool_sweeps/fiedler_tree_partition_by_eta.ipynb",
    # C2 — owned by C2, listed here only so the manifest is complete
    "pstar_gen_kmeans_2panel.png": "generated/eta_pool_sweeps/eta_pool_sweep.ipynb",
    "recovery_grid_kmeans.png": "generated/eta_pool_sweeps/eta_pool_sweep.ipynb",
    "pstar_gen_3operators.png": "generated/eta_pool_sweeps/eta_pool_sweep.ipynb",
    "identity_scatter_kingman.png": "generated/distance_vs_similarity/simulation_distance_vs_similarity.ipynb",
}

# Figures with a known copy living outside v9/figures/ that default mode should pull in.
ALT_SOURCES = {
    "coherence_vs_eta_by_n.png": ANALYSIS / "figures" / "coherence_vs_eta_by_n.png",
}


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def _mtime(path: Path) -> str:
    return datetime.datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")


def sync() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    for name, src in ALT_SOURCES.items():
        dst = FIGURES_DIR / name
        if src.exists() and (not dst.exists() or _md5(src) != _md5(dst)):
            shutil.copy2(src, dst)
            print(f"copied {src} -> {dst}")

    print(f"{'filename':<32}{'md5':<12}{'mtime':<22}notebook")
    for name, notebook in PAPER_FIGURES.items():
        path = FIGURES_DIR / name
        if not path.exists():
            print(f"{name:<32}{'MISSING':<12}{'-':<22}{notebook}")
            continue
        print(f"{name:<32}{_md5(path)[:10]:<12}{_mtime(path):<22}{notebook}")


# Grace windows in seconds. A legitimate `nbconvert --execute` writes the .ipynb
# (with fresh outputs) moments after its last `fig.savefig`, and a trailing
# markdown-only edit can nudge the notebook mtime a bit further — neither is
# real staleness. STALE_GRACE_S bounds the notebook-vs-figure gap;
# PROVENANCE_GRACE_S bounds how far apart figures from the *same* notebook may
# land, which is what would have caught the v8 HBM two-different-runs bug.
STALE_GRACE_S = 600
PROVENANCE_GRACE_S = 120


def check() -> int:
    problems = []
    by_notebook: dict[str, list[tuple[str, float]]] = {}
    for name, notebook in PAPER_FIGURES.items():
        path = FIGURES_DIR / name
        nb_path = ANALYSIS / notebook
        if not path.exists():
            problems.append(f"MISSING  {name}  (generator: {notebook})")
            continue
        mtime = path.stat().st_mtime
        by_notebook.setdefault(notebook, []).append((name, mtime))
        if nb_path.exists() and nb_path.stat().st_mtime - mtime > STALE_GRACE_S:
            problems.append(f"STALE    {name}  (notebook {notebook} modified >{STALE_GRACE_S}s after figure)")
    for notebook, entries in by_notebook.items():
        mtimes = [m for _, m in entries]
        if len(entries) > 1 and max(mtimes) - min(mtimes) > PROVENANCE_GRACE_S:
            names = ", ".join(n for n, _ in entries)
            problems.append(f"MIXED-RUN  {notebook}: outputs span >{PROVENANCE_GRACE_S}s ({names})")
    if problems:
        print("sync_paper_figures --check: problems found")
        for p in problems:
            print(f"  {p}")
        return 1
    print(f"sync_paper_figures --check: all {len(PAPER_FIGURES)} figures present and not stale")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="read-only staleness/missing check")
    args = parser.parse_args()
    if args.check:
        sys.exit(check())
    sync()
