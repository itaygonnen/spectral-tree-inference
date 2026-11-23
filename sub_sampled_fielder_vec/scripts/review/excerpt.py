"""Crop the real page image out of the built PDF, for each review item.

A paraphrase of a result is not the result. synctex maps a .tex line to a page and a
vertical position in the typeset PDF, so every item can show the passage as it is
printed -- numbering, equation, cross-reference colours and all.

Needs `latexmk -synctex=1`; run `--build` once if `thesis_v11.synctex.gz` is missing.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

from . import register

V11 = register.V11
PDF = V11 / "thesis_v11.pdf"
OUT = V11 / "review" / "app" / "excerpts"
SIDECAR = V11 / "review" / "app" / "excerpts.json"

DPI = 150
SCALE = DPI / 72.0
PAGE_H = 792.0            # letter
LEFT, WIDTH = 61.0, 490.0  # 1in margin less a hair, text width plus a hair
TOP_LIMIT, BOTTOM_LIMIT = 64.0, 730.0
PAD, MAX_SPAN = 13.0, 330.0

ENVS = ("theorem", "lemma", "proposition", "corollary", "definition",
        "assumption", "condition", "remark", "figure", "table", "algorithm")


def synctex(path: Path, line: int) -> dict | None:
    try:
        out = subprocess.run(
            ["synctex", "view", "-i", f"{line}:1:{path}", "-o", str(PDF)],
            capture_output=True, text=True, cwd=V11, timeout=30).stdout
    except Exception:
        return None
    got = {}
    for key in ("Page", "v", "H", "h", "W"):
        m = re.search(rf"^{key}:([0-9.]+)$", out, re.M)
        if m:
            got[key] = float(m.group(1))
    return {"page": int(got["Page"]), "v": got["v"], "h": got.get("h", LEFT),
            "H": got.get("H", 11.0)} if "Page" in got and "v" in got else None


def live_line(lines: list[str], start: int, direction: int) -> int | None:
    """Nearest line that actually prints. A commented line has no page position."""
    i = start
    while 1 <= i <= len(lines):
        text = lines[i - 1].strip()
        if text and not text.startswith("%") and not text.startswith("\\label"):
            return i
        i += direction
    return None


def span_for(path: Path, line: int) -> tuple[int, int]:
    """The lines to show: the whole environment if the anchor opens one, else its paragraph."""
    lines = path.read_text(encoding="utf-8").splitlines()
    idx = min(max(line, 1), len(lines))
    m = re.search(r"\\begin\{(" + "|".join(ENVS) + r")\}", lines[idx - 1])
    if m:
        env = m.group(1)
        for j in range(idx, min(len(lines), idx + 60)):
            if f"\\end{{{env}}}" in lines[j]:
                return idx, j + 1
        return idx, min(idx + 12, len(lines))
    first = idx
    while first > 1 and lines[first - 2].strip():
        first -= 1
    last = idx
    while last < len(lines) and lines[last].strip():
        last += 1
    return max(first, idx - 12), min(last, idx + 12)


def crop(page: int, top: float, height: float, dest: Path) -> bool:
    top = max(TOP_LIMIT, min(top, BOTTOM_LIMIT - 20))
    height = max(28.0, min(height, BOTTOM_LIMIT - top))
    args = ["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(DPI), "-png",
            "-x", str(int(LEFT * SCALE)), "-y", str(int(top * SCALE)),
            "-W", str(int(WIDTH * SCALE)), "-H", str(int(height * SCALE)),
            str(PDF), str(dest)]
    if subprocess.run(args, capture_output=True, cwd=V11).returncode != 0:
        return False
    made = dest.parent / f"{dest.name}-{page}.png"
    if not made.is_file():
        made = dest.parent / f"{dest.name}-{page:02d}.png"
    if made.is_file():
        made.rename(dest.with_suffix(".png"))
        return True
    return False


def one_item(item: dict) -> dict | None:
    f = item["fields"]
    m = register.ANCHOR_RE.match(f.get("Anchor", ""))
    if not m:
        return None
    rel, line = m.group(1), int(m.group(2))
    path = (V11 / rel).resolve()
    if not path.is_file() or path.suffix != ".tex" or path.parent.name != "sections":
        # The .bib and the master's preamble print nothing of their own. Falling back to
        # the nearest line that does print lands on the title page, which is worse than
        # showing no excerpt at all.
        return None

    lines = path.read_text(encoding="utf-8").splitlines()
    start, end = span_for(path, line)
    a = synctex(path, start) or (lambda n: synctex(path, n) if n else None)(
        live_line(lines, start, 1))
    if not a:
        return None
    b = synctex(path, end) or (lambda n: synctex(path, n) if n else None)(
        live_line(lines, end, -1)) or a

    slug = item["id"].replace("/", "-")
    images = []
    if b["page"] == a["page"] and b["v"] >= a["v"]:
        top = a["v"] - a["H"] - PAD
        height = min((b["v"] - a["v"]) + b["H"] + 2 * PAD, MAX_SPAN)
        if crop(a["page"], top, height, OUT / slug):
            images.append(f"excerpts/{slug}.png")
    else:
        # the passage runs over a page break: show the foot of one page and the head
        # of the next, which is also the honest picture of how it reads
        if crop(a["page"], a["v"] - a["H"] - PAD, MAX_SPAN, OUT / f"{slug}-a"):
            images.append(f"excerpts/{slug}-a.png")
        if b["page"] != a["page"] and crop(b["page"], TOP_LIMIT,
                                           min(b["v"] - TOP_LIMIT + b["H"] + PAD, MAX_SPAN),
                                           OUT / f"{slug}-b"):
            images.append(f"excerpts/{slug}-b.png")
    if not images:
        return None
    return {"page": a["page"], "pages": sorted({a["page"], b["page"]}), "images": images}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--build", action="store_true", help="rebuild the PDF with synctex first")
    args = ap.parse_args(argv)

    if args.build or not (V11 / "thesis_v11.synctex.gz").is_file():
        subprocess.run(["latexmk", "-pdf", "-synctex=1", "-g",
                        "-interaction=nonstopmode", "thesis_v11.tex"],
                       cwd=V11, capture_output=True)
    if not PDF.is_file():
        print("no built PDF", file=sys.stderr)
        return 2

    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.png"):
        old.unlink()
    found = {}
    for item in register.load_items():
        got = one_item(item)
        if got:
            found[item["id"]] = got
        print(f"  {item['id']:8} {'page ' + str(got['page']) if got else 'no page image'}")
    SIDECAR.write_text(json.dumps(found, indent=2) + "\n", encoding="utf-8")
    print(f"{len(found)} excerpts -> {OUT.relative_to(register.REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
