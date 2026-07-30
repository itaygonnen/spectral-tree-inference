#!/usr/bin/env python3
"""Validate the v9 paper's figure provenance, and regenerate PAPER_MAP.md.

The manifest lives in ``analysis/paper_figures.py``.
This script only reads it, so mapping / docs / checker cannot drift.

    python scripts/sync_paper_figures.py --check      # read-only; non-zero on problems
    python scripts/sync_paper_figures.py --write-map  # regenerate PAPER_MAP.md
    python scripts/sync_paper_figures.py              # manifest listing

Checks
------
MISSING     a mapped figure is absent from v9/figures/
UNMAPPED    an \\includegraphics in v9/sections/*.tex has no manifest entry
                -- the check that catches a wrong or stale mapping, which is how the
                   previous hand-maintained version drifted unnoticed
ORPHAN      a file in v9/figures/ is neither mapped nor listed as RETIRED
NO-NOTEBOOK the producing notebook does not exist at the recorded path
STALE       the notebook was modified well after the figure was written
MIXED-RUN   panels of ONE float were written far apart in time, i.e. that figure
                mixes two runs. This is the check that would have caught the v8 HBM
                bug where two panels of one figure came from different runs
                (open-items/00-R1.md [R1/23]).

                Grouped by float, NOT by notebook. One notebook legitimately produces
                several independent figures across different runs -- eta_pool_sweep
                emits Fig 3 (§5) plus Figs 8 and 9 (App G), and those were written 21
                minutes apart while each float stayed internally coherent. Grouping by
                notebook flagged that as a defect; it is not one.

STALE / MIXED-RUN are mtime heuristics, so they WARN rather than fail: a `git clone`
gives every file the same checkout time (hiding real staleness), and merely editing a
notebook's markdown bumps its mtime past its figures. MISSING / UNMAPPED /
UNUSED-MAP / ORPHAN / NO-NOTEBOOK are exact and gate the exit code.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis"
V9 = ROOT / "docs" / "overleafs" / "v9"
FIGURES_DIR = V9 / "figures"
SECTIONS_DIR = V9 / "sections"
MAP_PATH = ANALYSIS / "PAPER_MAP.md"

sys.path.insert(0, str(ANALYSIS))
import paper_figures as PF  # noqa: E402

# A legitimate `nbconvert --execute` writes the .ipynb moments after its last
# savefig, and a trailing markdown edit can nudge the mtime further -- neither is
# real staleness. PROVENANCE_GRACE_S bounds how far apart figures from the *same*
# notebook may land.
STALE_GRACE_S = 600
PROVENANCE_GRACE_S = 120


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def _mtime(path: Path) -> str:
    return datetime.datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")


def tex_images() -> dict[str, list[str]]:
    """Every \\includegraphics target in v9/sections/, keyed by filename."""
    found: dict[str, list[str]] = {}
    for tex in sorted(SECTIONS_DIR.glob("*.tex")):
        for img in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex.read_text()):
            found.setdefault(img, []).append(f"sections/{tex.name}")
    return found


def check() -> int:
    errors: list[str] = []
    warnings: list[str] = []
    mapped = PF.by_file()

    # --- exact checks: these gate the exit code -----------------------------
    for name, fig in mapped.items():
        if not (FIGURES_DIR / name).exists():
            errors.append(f"MISSING      {name}  (Fig {fig.number}, from {fig.notebook})")
        if not (ANALYSIS / fig.notebook).exists():
            errors.append(f"NO-NOTEBOOK  {fig.notebook}  (Fig {fig.number})")

    in_tex = tex_images()
    for img, where in in_tex.items():
        if img not in mapped:
            errors.append(f"UNMAPPED     {img}  used in {', '.join(where)} but absent from paper_figures.py")
    for name in mapped:
        if name not in in_tex:
            errors.append(f"UNUSED-MAP   {name}  mapped but no \\includegraphics references it")

    on_disk = {
        str(p.relative_to(FIGURES_DIR))
        for p in FIGURES_DIR.rglob("*.png")
        if "retired" not in p.parts
    }
    for name in sorted(on_disk - set(mapped) - set(PF.RETIRED)):
        errors.append(f"ORPHAN       {name}  in v9/figures/ but neither mapped nor RETIRED")

    # Every notebook must be accounted for: either it produces a figure or it is listed as
    # supporting. Without this, a new notebook is invisible to the map -- which is exactly
    # how the two distance-route notebooks went unlisted.
    listed = {f.notebook for f in PF.PAPER_FIGURES} | {s.notebook for s in PF.SUPPORTING}
    present = {
        str(p.relative_to(ANALYSIS))
        for p in ANALYSIS.rglob("*.ipynb")
        if "legacy" not in p.parts and ".ipynb_checkpoints" not in p.parts
    }
    for nb in sorted(present - listed):
        errors.append(f"UNLISTED-NB  {nb}  exists but is in neither PAPER_FIGURES nor SUPPORTING")
    for nb in sorted(listed - present):
        errors.append(f"GHOST-NB     {nb}  listed in the manifest but not on disk")

    # --- mtime heuristics: advisory only ------------------------------------
    # Grouped by FLOAT (fig.label), not by notebook: panels of one figure must come
    # from one run, but one notebook may legitimately produce several figures across
    # different runs.
    by_float: dict[str, list[tuple[str, float]]] = {}
    for name, fig in mapped.items():
        path = FIGURES_DIR / name
        if not path.exists():
            continue
        by_float.setdefault(fig.label, []).append((name, path.stat().st_mtime))
        nb = ANALYSIS / fig.notebook
        if nb.exists() and nb.stat().st_mtime - path.stat().st_mtime > STALE_GRACE_S:
            warnings.append(f"STALE        {name}  (notebook modified >{STALE_GRACE_S}s after figure)")
    for label, entries in by_float.items():
        times = [t for _, t in entries]
        if len(entries) > 1 and max(times) - min(times) > PROVENANCE_GRACE_S:
            names = ", ".join(n for n, _ in entries)
            span = int(max(times) - min(times))
            warnings.append(
                f"MIXED-RUN    {label}: panels span {span}s "
                f"(>{PROVENANCE_GRACE_S}s) -- one float, two runs? ({names})")

    for w in warnings:
        print(f"  warn  {w}")
    if warnings:
        print(f"\n{len(warnings)} advisory warning(s) -- mtime-based, do not gate\n")
    if errors:
        print("sync_paper_figures --check: FAILED\n")
        for e in errors:
            print(f"  {e}")
        print(f"\n{len(errors)} error(s)")
        return 1
    print(f"sync_paper_figures --check: OK -- {len(mapped)} figures across "
          f"{len(PF.PAPER_FIGURES)} floats, all present, mapped and cross-checked "
          f"against {len(in_tex)} \\includegraphics targets")
    return 0


def listing() -> None:
    print(f"{'figure':<52}{'md5':<12}{'mtime':<22}notebook")
    for fig in PF.PAPER_FIGURES:
        for name in fig.files:
            path = FIGURES_DIR / name
            if not path.exists():
                print(f"{name:<52}{'MISSING':<12}{'-':<22}{fig.notebook}")
            else:
                print(f"{name:<52}{_md5(path)[:10]:<12}{_mtime(path):<22}{fig.notebook}")


def write_map() -> None:
    L: list[str] = []
    A = L.append
    A("# PAPER_MAP — figure provenance for the v9 manuscript")
    A("")
    A("<!-- GENERATED FILE. Edit paper_figures.py, then run:")
    A("     python scripts/sync_paper_figures.py --write-map -->")
    A("")
    A("Where every figure in `docs/overleafs/v9/` comes from: which notebook produces it,")
    A("which claim it speaks to, which cache or results dir it consumes, and how to rebuild")
    A("it. Validated by `python scripts/sync_paper_figures.py --check`.")
    A("")
    A("## Paper figures")
    A("")
    A("| Fig | File(s) | Section | Label | Notebook | Claim(s) |")
    A("|---|---|---|---|---|---|")
    for f in PF.PAPER_FIGURES:
        files = "<br>".join(f"`{x}`" for x in f.files)
        claims = ", ".join(f"`{c}`" for c in f.claims) or "—"
        A(f"| **{f.number}** | {files} | `{f.section}` | `{f.label}` "
          f"| `{f.notebook}` | {claims} |")
    A("")
    A("## Inputs and rebuild")
    A("")
    for f in PF.PAPER_FIGURES:
        A(f"### Figure {f.number} — `{f.label}`")
        A("")
        A(f"- **Notebook**: `{f.notebook}`")
        A(f"- **Consumes**: {f.cache}")
        A(f"- **Rebuild**: `{f.rebuild}`")
        if f.notes:
            A(f"- **Note**: {f.notes}")
        A("")
    A("## Floats with no image file")
    A("")
    for label, what in PF.NON_IMAGE_FLOATS.items():
        A(f"- `{label}` — {what}")
    A("")
    A("## Retired figures")
    A("")
    A("Present in `v9/figures/` (or `figures/retired/`) but deliberately unreferenced.")
    A("Listed so a stray file reads as \"retired on purpose\", not \"someone forgot\".")
    A("")
    A("| File | Why retired |")
    A("|---|---|")
    for name, why in PF.RETIRED.items():
        A(f"| `{name}` | {why} |")
    A("")
    A("## Supporting notebooks (no paper figure)")
    A("")
    A("All 10 are kept. This table exists so nobody mistakes them for figure sources,")
    A("and so the evidence behind the distance route stays findable — App F /")
    A("`thm:main-dist` carries no figure of its own.")
    A("")
    A("| Notebook | Backs | Status |")
    A("|---|---|---|")
    for s in PF.SUPPORTING:
        A(f"| `{s.notebook}` | {s.backs} | {s.status} |")
    A("")
    MAP_PATH.write_text("\n".join(L) + "\n")
    print(f"wrote {MAP_PATH.relative_to(ROOT)}  "
          f"({len(PF.PAPER_FIGURES)} floats, {len(PF.all_files())} files, "
          f"{len(PF.SUPPORTING)} supporting notebooks)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="read-only validation")
    ap.add_argument("--write-map", action="store_true", help="regenerate PAPER_MAP.md")
    args = ap.parse_args()
    if args.check:
        sys.exit(check())
    if args.write_map:
        write_map()
    else:
        listing()
