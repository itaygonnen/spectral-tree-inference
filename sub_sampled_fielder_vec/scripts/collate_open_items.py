#!/usr/bin/env python3
"""Assemble v9/open-items/*.md fragments into the single register v9/OPEN_ITEMS.md.

Each restructure agent writes only its own fragment; concurrent appends to one file
would clobber. Fragments are concatenated in filename order (hence the NN- prefixes).

Also lints the mandatory entry format and reports violations, so a malformed entry is
caught here rather than by the reader:
    ### [ID/NN] title
    - **Anchor:**   required
    - **Type:**     question | doubt | dilemma | flaw | direction
    - **Item:**     required, <= 3 sentences
    - **Options:**  required iff Type is dilemma, forbidden otherwise
    - **Blocking:** required, must start yes/no

Usage:
    python scripts/collate_open_items.py            # write the register
    python scripts/collate_open_items.py --check     # lint only, no write
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs" / "overleafs" / "v9"
FRAGMENT_DIR = DOCS / "open-items"
REGISTER = DOCS / "OPEN_ITEMS.md"

VALID_TYPES = {"question", "doubt", "dilemma", "flaw", "direction"}
ENTRY_RE = re.compile(r"^### \[([A-Za-z0-9]+)/(\d+)\]\s+(.+)$")
FIELD_RE = re.compile(r"^-\s+\*\*(\w+):\*\*\s*(.*)$")

HEADER = """# v9 — open items register

Single register for everything the restructure could not settle: unresolved questions, doubts,
dilemmas, suspected flaws, and directions worth pursuing. One entry per item.

`Blocking: no` entries record the default that was taken, so every one is visible and reversible.

**This file is generated** — edit the fragments in `open-items/` and re-run
`python scripts/collate_open_items.py`.

"""


def split_sentences(text: str) -> list[str]:
    """Rough sentence count; good enough to catch entries that are really three items."""
    text = re.sub(r"`[^`]*`", "X", text)          # code spans may contain periods
    text = re.sub(r"\b(e\.g|i\.e|cf|vs|resp)\.", r"\1", text)
    return [s for s in re.split(r"(?<=[.?!])\s+", text.strip()) if s]


def lint(path: Path) -> tuple[int, list[str]]:
    """Return (entry count, problems) for one fragment."""
    problems: list[str] = []
    lines = path.read_text(encoding="utf-8").splitlines()

    entries: list[tuple[int, str]] = [
        (n, m.group(0)) for n, line in enumerate(lines, 1) if (m := ENTRY_RE.match(line))
    ]
    if not entries:
        return 0, problems

    bounds = [n for n, _ in entries] + [len(lines) + 1]
    for idx, (start, heading) in enumerate(entries):
        body = lines[start : bounds[idx + 1] - 1]
        fields: dict[str, str] = {}
        for line in body:
            if m := FIELD_RE.match(line.strip()):
                fields[m.group(1).lower()] = m.group(2)
            elif fields and line.startswith("  ") and line.strip():
                last = list(fields)[-1]                  # continuation line
                fields[last] += " " + line.strip()

        tag = f"{path.name}:{start} {heading[:60]}"
        for required in ("anchor", "type", "item", "blocking"):
            if required not in fields or not fields[required].strip():
                problems.append(f"{tag} — missing **{required.title()}:**")

        kind = fields.get("type", "").strip().lower()
        if kind and kind not in VALID_TYPES:
            problems.append(f"{tag} — Type '{kind}' not in {sorted(VALID_TYPES)}")
        if kind == "dilemma" and "options" not in fields:
            problems.append(f"{tag} — Type is dilemma but **Options:** is absent")
        if kind and kind != "dilemma" and "options" in fields:
            problems.append(f"{tag} — **Options:** present but Type is '{kind}', not dilemma")

        blocking = fields.get("blocking", "").strip().lower()
        if blocking and not blocking.startswith(("yes", "no")):
            problems.append(f"{tag} — Blocking must start with yes/no, got '{blocking[:20]}'")

        n_sent = len(split_sentences(fields.get("item", "")))
        if n_sent > 3:
            problems.append(f"{tag} — Item has ~{n_sent} sentences (max 3)")

    return len(entries), problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="lint only; do not write the register")
    args = ap.parse_args()

    if not FRAGMENT_DIR.is_dir():
        print(f"no fragment directory at {FRAGMENT_DIR}", file=sys.stderr)
        return 1

    fragments = sorted(FRAGMENT_DIR.glob("*.md"))
    if not fragments:
        print(f"no fragments in {FRAGMENT_DIR}", file=sys.stderr)
        return 1

    parts, total, all_problems = [], 0, []
    for frag in fragments:
        count, problems = lint(frag)
        total += count
        all_problems += problems
        parts.append(frag.read_text(encoding="utf-8").rstrip() + "\n")
        print(f"  {frag.name:24s} {count:3d} entries" + ("" if not problems else f"  ({len(problems)} problems)"))

    for p in all_problems:
        print(f"  LINT {p}")

    if not args.check:
        REGISTER.write_text(HEADER + "\n---\n\n".join(parts), encoding="utf-8")
        print(f"\nwrote {REGISTER.relative_to(DOCS.parent.parent)} — {total} entries from {len(fragments)} fragments")
    else:
        print(f"\n{total} entries from {len(fragments)} fragments; {len(all_problems)} lint problems")

    return 1 if all_problems else 0


if __name__ == "__main__":
    sys.exit(main())
