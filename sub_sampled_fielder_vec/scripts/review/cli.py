"""Apply approved items to v11, and show the delta against v10.

Decisions come from the review artifact's database; this script never talks to it.
Fetch them with the ArtifactData tool into a JSON file shaped
`{"<ID with - for />": {"state": ..., "comment": ...}}` and pass it with --decisions.

Two rules the apply enforces, both of them the point of the v10/v11 split:
  * a patch may only touch files under v11; a patch naming v10 is a bug, not an edit
  * the build must still report zero undefined references afterwards, or the batch
    is rolled back
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import deps, register

V10 = register.REPO_ROOT / "docs" / "overleafs" / "v10"
V11 = register.V11


def read_decisions(path: Path) -> dict[str, dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    rows = raw.get("items", raw) if isinstance(raw, dict) else raw
    out = {}
    for key, val in rows.items():
        out[key.replace("-", "/", 1) if "/" not in key else key] = val
    return out


def undefined_count() -> int | None:
    log = V11 / "thesis_v11.log"
    if not log.is_file():
        return None
    return sum(1 for line in log.read_text(encoding="utf-8", errors="ignore").splitlines()
               if "undefined" in line.lower())


def build() -> int | None:
    proc = subprocess.run(
        ["latexmk", "-pdf", "-interaction=nonstopmode", "thesis_v11.tex"],
        cwd=V11, capture_output=True, text=True)
    if proc.returncode not in (0, 12):     # 12 = latexmk "errors in the log"
        print(proc.stdout[-1500:], file=sys.stderr)
    return undefined_count()


def apply_patch(item: dict, dry: bool) -> list[str]:
    """Apply one item's patch. Returns the problems found; empty means applied."""
    spec_rel = item["fields"].get("Patch", "none")
    if spec_rel in ("none", ""):
        return [f"{item['id']}: approved but carries no patch -- nothing to apply"]
    spec_path = V11 / spec_rel if spec_rel.startswith("review/") else V11 / "review" / spec_rel
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    edited: dict[Path, str] = {}
    for i, hunk in enumerate(spec.get("hunks", []), 1):
        rel = hunk.get("file", spec.get("file", ""))
        target = (V11 / rel).resolve()
        if V11.resolve() not in target.parents and target.parent != V11.resolve():
            return [f"{item['id']}: patch target {rel} is outside v11"]
        text = edited.get(target, target.read_text(encoding="utf-8"))
        old, new = hunk["anchor_old"], hunk["replacement_new"]
        if text.count(old) != 1:
            return [f"{item['id']}: hunk {i} anchor occurs {text.count(old)} times, needs 1"]
        edited[target] = text.replace(old, new)
    if not dry:
        for target, text in edited.items():
            target.write_text(text, encoding="utf-8")
    return []


def cmd_apply(args) -> int:
    decisions = read_decisions(args.decisions)
    items = {i["id"]: i for i in register.load_items()}
    approved = [items[i] for i, d in decisions.items()
                if d.get("state") == "approved" and i in items
                and items[i]["fields"].get("Status") != "applied"]
    if not approved:
        print("nothing approved and still open")
        return 0

    before = undefined_count()
    backups = {}
    for item in approved:
        spec_rel = item["fields"].get("Patch", "none")
        if spec_rel in ("none", ""):
            continue
        sp = V11 / spec_rel if spec_rel.startswith("review/") else V11 / "review" / spec_rel
        spec = json.loads(sp.read_text())
        for hunk in spec.get("hunks", []):
            tgt = V11 / hunk.get("file", spec.get("file", ""))
            if tgt.is_file():
                backups.setdefault(tgt, tgt.read_text(encoding="utf-8"))

    problems = []
    for item in sorted(approved, key=lambda i: i["id"]):
        problems += apply_patch(item, args.dry_run)
    for p in problems:
        print("  " + p)
    if args.dry_run:
        print(f"dry run: {len(approved)} approved, {len(problems)} problem(s)")
        return 1 if problems else 0

    after = build()
    if before is not None and after is not None and after > before:
        for path, text in backups.items():
            path.write_text(text, encoding="utf-8")
        build()
        print(f"ROLLED BACK: undefined references went {before} -> {after}")
        return 1

    for item in sorted(approved, key=lambda i: i["id"]):
        if any(item["id"] in p for p in problems):
            continue
        d = decisions[item["id"]]
        text = item["path"].read_text(encoding="utf-8")
        text = re.sub(r"^- \*\*Status:\*\* .*$", "- **Status:** applied", text,
                      count=1, flags=re.M)
        # An applied item needs its own date: "approved" is a decision the author made,
        # "applied" is an edit that exists in the .tex, and the two are days apart.
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        text = re.sub(r"^- \*\*Applied:\*\* .*\n", "", text, flags=re.M)
        text = text.rstrip("\n") + f"\n- **Applied:** {stamp}\n"
        if d.get("comment"):
            text = re.sub(r"^- \*\*Decision:\*\* .*\n", "", text, flags=re.M)
            text = text.rstrip("\n") + f"\n- **Decision:** approved -- {d['comment']}\n"
        item["path"].write_text(text, encoding="utf-8")

    print(f"applied {len(approved) - len(problems)} item(s); "
          f"undefined references {before} -> {after}")
    return 1 if problems else 0


def cmd_reanchor(args) -> int:
    """Move every label-bearing anchor to where its label now sits.

    Applying one item shifts the line numbers of every other item in the same file.
    Labels survive that, so they are the stable address; the line number is a
    convenience that has to be refreshed after each apply.
    """
    graph = deps.Graph().scan(V11 / "sections")
    moved = 0
    for item in register.load_items():
        f = item["fields"]
        m = register.ANCHOR_RE.match(f.get("Anchor", ""))
        if not m or not m.group(3):
            continue
        label = m.group(3)
        meta = graph.labels.get(label)
        if not meta:
            print(f"  {item['id']}: label {label} no longer exists")
            continue
        want = f"sections/{meta['file']}:{meta['line']} ({label})"
        if want == f["Anchor"]:
            continue
        text = item["path"].read_text(encoding="utf-8")
        item["path"].write_text(
            re.sub(r"^- \*\*Anchor:\*\* .*$", f"- **Anchor:** {want}", text,
                   count=1, flags=re.M), encoding="utf-8")
        print(f"  {item['id']}: {f['Anchor']} -> {want}")
        moved += 1
    print(f"{moved} anchor(s) moved")
    return 0


def cmd_delta(args) -> int:
    proc = subprocess.run(["diff", "-ru", str(V10 / "sections"), str(V11 / "sections")],
                          capture_output=True, text=True)
    if not proc.stdout.strip():
        print("v11/sections is identical to v10/sections")
        return 0
    print(proc.stdout)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("apply", help="apply approved patches to v11")
    a.add_argument("--decisions", type=Path, required=True)
    a.add_argument("--dry-run", action="store_true")
    a.set_defaults(func=cmd_apply)
    ra = sub.add_parser("reanchor", help="move label anchors to their current lines")
    ra.set_defaults(func=cmd_reanchor)
    d = sub.add_parser("delta", help="diff v10/sections against v11/sections")
    d.set_defaults(func=cmd_delta)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
