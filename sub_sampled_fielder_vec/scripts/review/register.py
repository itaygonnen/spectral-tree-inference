"""Parse, lint, collate and export the v11 review items.

The field grammar is v9's, proven by `scripts/collate_open_items.py`: a `### [ID]`
heading, `- **Field:** value` lines, two-space continuations. What is new is one
item per FILE (the app needs per-item identity) and four fields the skills
themselves supply: Category, Source, Patch and Verified.

`Verified:` exists because a tool's finding is a claim, not a fact. The first
latex-paper-en run on v11 reported six P0 "undefined reference" criticals; all six
were its own failure to split `\\Cref{a,b}`, and the build has zero undefined refs.
An unverified finding never becomes a defect in front of the author.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from . import deps

REPO_ROOT = Path(__file__).resolve().parents[2]
V11 = REPO_ROOT / "docs" / "overleafs" / "v11"
ITEMS_DIR = V11 / "review" / "items"
PATCH_DIR = V11 / "review" / "patches"
COLLATED = V11 / "review" / "REVIEW_ITEMS.md"
APP_JSON = V11 / "review" / "app" / "items.json"

ENTRY_RE = re.compile(r"^### \[([A-Za-z0-9]+)/(\d+)\]\s+(.+)$")
FIELD_RE = re.compile(r"^-\s+\*\*(\w+):\*\*\s*(.*)$")
# Anchors are relative to v11/: a sections file, the master, or references.bib.
ANCHOR_RE = re.compile(r"^([\w./-]+\.(?:tex|bib)):(\d+)(?:\s+\(([^)]+)\))?$")

REQUIRED = ["Anchor", "Source", "Category", "Type", "Severity", "Context", "Item",
            "Affects", "Patch", "Status", "Verified"]
OPTIONAL = ["Context", "Coverage", "Applied", "Decision"]
# Context says what the thing at the anchor actually SAYS, in plain words. Without it
# an item names a label and makes the author go find it -- which is how a register
# becomes unreadable.  # Applied/Decision written by `cli.py apply`
# Coverage answers "does this patch catch every occurrence?" -- for a notation change,
# fixing one site and leaving the rest is worse than leaving all of them.
TYPES = {"phrasing", "notation", "citation", "claim-discipline", "structure",
         "proof-gap", "statement", "assumption", "figure"}
SEVERITIES = {"blocking", "material", "editorial"}
STATUSES = {"open", "approved", "rejected", "needs-work", "applied"}
# proof-checker's derived labels -> this register's three
SEVERITY_MAP = {"FATAL": "blocking", "CRITICAL": "blocking",
                "MAJOR": "material", "MINOR": "editorial"}


def parse_item(path: Path) -> dict:
    """One item file -> a dict. Raises ValueError with the file name on a malformed file."""
    item: dict = {"path": path, "fields": {}}
    field = None
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.rstrip()
        if not line:
            continue
        head = ENTRY_RE.match(line)
        if head:
            if item.get("id"):
                raise ValueError(f"{path.name}: a second `### [ID]` heading at line {lineno}")
            item["tag"], item["num"], item["title"] = head.group(1), head.group(2), head.group(3)
            item["id"] = f"{head.group(1)}/{head.group(2)}"
            continue
        fld = FIELD_RE.match(line)
        if fld:
            field = fld.group(1)
            item["fields"][field] = fld.group(2).strip()
        elif field and raw.startswith("  "):
            item["fields"][field] += " " + line.strip()
        else:
            raise ValueError(f"{path.name}:{lineno}: not a heading, field or continuation")
    if not item.get("id"):
        raise ValueError(f"{path.name}: no `### [ID/NN] title` heading")
    return item


def load_items() -> list[dict]:
    if not ITEMS_DIR.is_dir():
        return []
    return [parse_item(p) for p in sorted(ITEMS_DIR.glob("*.md"))]


def _sentences(text: str) -> int:
    return len([s for s in re.split(r"(?<=[.!?])\s+", text.strip()) if s])


def lint(items: list[dict], graph: deps.Graph) -> list[str]:
    problems: list[str] = []
    seen: dict[str, str] = {}

    def bad(item, msg):
        problems.append(f"{item['path'].name} [{item.get('id','?')}] {msg}")

    for item in items:
        f = item["fields"]
        if item["id"] in seen:
            bad(item, f"duplicate id, also in {seen[item['id']]}")
        seen[item["id"]] = item["path"].name

        for key in REQUIRED:
            if key not in f:
                bad(item, f"missing required field {key}")
        if f.get("Type") and f["Type"] not in TYPES:
            bad(item, f"Type {f['Type']!r} not one of {sorted(TYPES)}")
        if f.get("Severity") and f["Severity"] not in SEVERITIES:
            bad(item, f"Severity {f['Severity']!r} not one of {sorted(SEVERITIES)}")
        if f.get("Status") and f["Status"] not in STATUSES:
            bad(item, f"Status {f['Status']!r} not one of {sorted(STATUSES)}")
        if f.get("Item") and _sentences(f["Item"]) > 3:
            bad(item, f"Item is {_sentences(f['Item'])} sentences, max 3")
        if f.get("Verified", "").lower().startswith("no"):
            bad(item, "Verified is 'no' -- check the claim against the source "
                      "before this item reaches the author, or delete it")

        m = ANCHOR_RE.match(f.get("Anchor", ""))
        if not m:
            bad(item, f"Anchor {f.get('Anchor')!r} is not 'sections/<file>.tex:<line> (label)'")
        else:
            rel, line, label = m.group(1), int(m.group(2)), m.group(3)
            target = V11 / rel
            if not target.is_file():
                bad(item, f"Anchor file {rel} does not exist")
            elif (line > len(target.read_text(encoding="utf-8").splitlines())
                  and f.get("Status") != "applied"):
                # Applying one item shifts the line numbers of every other item in the
                # same file, so a stale anchor on an OPEN item is a real defect to fix.
                bad(item, f"Anchor line {line} is past the end of {rel} -- re-anchor it")
            if label and label not in graph.labels:
                bad(item, f"Anchor label {label!r} is not a label in v11")
            elif label:
                computed = {p[-1] for p in graph.affects(label)}
                claimed = {t.strip() for chain in f.get("Affects", "").split(";")
                           for t in chain.split("->")} - {"", "none"}
                stray = claimed - computed - {label}
                if stray:
                    bad(item, f"Affects names {sorted(stray)}, which the graph does not derive")

        patch = f.get("Patch", "none")
        # An applied patch's anchor is gone by definition -- that is what applied means.
        # Checking it would turn every successful apply into a lint failure.
        if f.get("Status") == "applied":
            patch = "none"
        if patch not in ("none", ""):
            pp = V11 / "review" / patch if not patch.startswith("review/") else V11 / patch
            if not pp.is_file():
                bad(item, f"Patch {patch} does not exist")
            else:
                try:
                    spec = json.loads(pp.read_text(encoding="utf-8"))
                except json.JSONDecodeError as exc:
                    bad(item, f"Patch {patch} is not valid JSON: {exc}")
                    continue
                # A hunk may name its own file, so one item can carry a change that
                # runs across sections -- a notation or punctuation fix usually does.
                for i, hunk in enumerate(spec.get("hunks", []), 1):
                    rel = hunk.get("file", spec.get("file", ""))
                    src = V11 / rel
                    if not src.is_file():
                        bad(item, f"Patch hunk {i} target {rel!r} does not exist")
                        continue
                    text = src.read_text(encoding="utf-8")
                    n = text.count(hunk.get("anchor_old", "\0"))
                    if n == 0:
                        bad(item, f"Patch hunk {i}: anchor_old not found in {rel}")
                    elif n > 1:
                        bad(item, f"Patch hunk {i}: anchor_old occurs {n} times -- not unique")
    return problems


def collate(items: list[dict]) -> str:
    by_sev = {s: [i for i in items if i["fields"].get("Severity") == s] for s in
              ("blocking", "material", "editorial")}
    out = ["# v11 review items", "",
           "> Generated by `scripts/review/register.py collate`. Do not hand-edit:",
           "> edit the files in `review/items/` and regenerate.", "",
           f"{len(items)} items: "
           + ", ".join(f"{len(v)} {k}" for k, v in by_sev.items()), ""]
    for sev in ("blocking", "material", "editorial"):
        if not by_sev[sev]:
            continue
        out += [f"## {sev}", ""]
        for item in sorted(by_sev[sev], key=lambda i: i["id"]):
            f = item["fields"]
            out.append(f"### [{item['id']}] {item['title']}")
            for key in REQUIRED + OPTIONAL:
                if key in f:
                    out.append(f"- **{key}:** {f[key]}")
            out.append("")
    return "\n".join(out) + "\n"


def pathlib_name(p: str | None) -> str | None:
    """Anchors are written `sections/x.tex`; the map keys on the bare file name."""
    return p.rsplit("/", 1)[-1] if p else None


def related_map(rows: list[dict]) -> dict[str, list[dict]]:
    """Which items should be judged together, and WHY.

    The reason is the useful half: "same statement" means one of the two edits may
    undo the other, while "same section" only means you are already reading there.
    Strongest relation wins, so a pair is never listed under the weaker one.
    """
    rel: dict[str, list[dict]] = {r["id"]: [] for r in rows}
    for a in rows:
        for b in rows:
            if a["id"] == b["id"]:
                continue
            if a["label"] and a["label"] == b["label"]:
                why = "same statement"
            elif b["label"] and any(b["label"] in chain.split(" -> ")
                                    for chain in a["affects"]):
                why = "depends on this"
            elif a["label"] and any(a["label"] in chain.split(" -> ")
                                    for chain in b["affects"]):
                why = "this depends on it"
            elif a["file"] and a["file"] == b["file"]:
                why = "same section"
            else:
                continue
            rel[a["id"]].append({"id": b["id"], "why": why, "file": pathlib_name(b["file"]),
                                 "line": b["line"], "severity": b["severity"]})
    order = {"same statement": 0, "depends on this": 1, "this depends on it": 2, "same section": 3}
    for k in rel:
        rel[k].sort(key=lambda r: (order[r["why"]], r["id"]))
    return rel


def export(items: list[dict], graph: deps.Graph) -> dict:
    # Page images are produced separately by excerpt.py, which needs the built PDF.
    sidecar = V11 / "review" / "app" / "excerpts.json"
    shots = json.loads(sidecar.read_text(encoding="utf-8")) if sidecar.is_file() else {}
    rows = []
    for item in sorted(items, key=lambda i: i["id"]):
        f = item["fields"]
        m = ANCHOR_RE.match(f.get("Anchor", ""))
        label = m.group(3) if m and m.group(3) else None
        patch = f.get("Patch", "none")
        hunks = []
        if patch not in ("none", ""):
            pp = V11 / patch if patch.startswith("review/") else V11 / "review" / patch
            if pp.is_file():
                spec = json.loads(pp.read_text(encoding="utf-8"))
                hunks = [{"file": h.get("file", spec.get("file")),
                          "before": h.get("anchor_old", ""),
                          "after": h.get("replacement_new", "")}
                         for h in spec.get("hunks", [])]
        rows.append({
            "id": item["id"],
            "title": item["title"],
            "file": m.group(1) if m else None,
            "line": int(m.group(2)) if m else None,
            "label": label,
            "source": f.get("Source"),
            "category": f.get("Category"),
            "type": f.get("Type"),
            "severity": f.get("Severity"),
            "context": f.get("Context"),
            "excerpt": shots.get(item["id"]),
            "note": f.get("Item"),
            "verified": f.get("Verified"),
            "status": f.get("Status"),
            "coverage": f.get("Coverage"),
            "applied_on": f.get("Applied"),
            "decision_note": f.get("Decision"),
            "affects": [" -> ".join(p) for p in graph.affects(label)] if label else [],
            # Where the thing is SAID, and everywhere it is USED -- the two marks the
            # document map draws. Without the second, a citation item cannot show the
            # author what reading the change would mean elsewhere.
            "stated_at": ({"file": pathlib_name(graph.labels[label]["file"]),
                           "line": graph.labels[label]["line"]} if label in graph.labels
                          else {"file": pathlib_name(m.group(1)) if m else None,
                                "line": int(m.group(2)) if m else None}),
            "cited_at": graph.cite_sites.get(label, []) if label else [],
            "prose_sites": graph.prose_sites(label) if label else [],
            "proof_at": graph.proofs.get(label) if label else None,
            "statement": graph.statement_tex(label) if label else None,
            "hunks": hunks,
        })
    for row in rows:
        row["related"] = []
    rel = related_map(rows)
    for row in rows:
        row["related"] = rel[row["id"]]
    return {"generated_from": "docs/overleafs/v11/review/items",
            "document": deps.document(V11 / "sections", V11 / "thesis_v11.tex"),
            "items": rows,
            "counts": {s: sum(1 for r in rows if r["severity"] == s)
                       for s in ("blocking", "material", "editorial")}}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("command", choices=["lint", "collate", "export"])
    args = ap.parse_args(argv)

    graph = deps.Graph().scan(V11 / "sections")
    try:
        items = load_items()
    except ValueError as exc:
        print(f"malformed item: {exc}", file=sys.stderr)
        return 2

    if args.command == "lint":
        problems = lint(items, graph)
        for p in problems:
            print(f"  {p}")
        print(f"{len(items)} items, {len(problems)} problem(s)")
        return 1 if problems else 0

    if args.command == "collate":
        COLLATED.parent.mkdir(parents=True, exist_ok=True)
        COLLATED.write_text(collate(items), encoding="utf-8")
        print(f"wrote {COLLATED.relative_to(REPO_ROOT)} ({len(items)} items)")
        return 0

    APP_JSON.parent.mkdir(parents=True, exist_ok=True)
    APP_JSON.write_text(json.dumps(export(items, graph), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {APP_JSON.relative_to(REPO_ROOT)} ({len(items)} items)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
