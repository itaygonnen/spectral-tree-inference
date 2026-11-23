"""Label dependency graph over the v11 manuscript.

Why this exists: an item that changes a lemma must name the results that now need
re-checking. That list is computable -- every numbered object in the manuscript
carries a \\label, and every citation of it is a \\Cref/\\cref/\\ref/\\eqref.

Ownership heuristic (documented because it is a heuristic, not a parse):
a citation belongs to the innermost numbered statement whose \\label was last
seen in that file, and a `proof` environment is attributed to the statement it
follows. Outside any statement, a citation belongs to its section.

cleveref takes comma lists -- \\Cref{asm:gtr,asm:clock} cites TWO labels. Splitting
on the comma is not optional: not splitting it reports both labels as undefined,
which is how a clean manuscript acquires six phantom P0 criticals.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SECTIONS = REPO_ROOT / "docs" / "overleafs" / "v11" / "sections"

STATEMENT_ENVS = {
    "theorem", "lemma", "proposition", "corollary",
    "definition", "assumption", "condition", "remark",
}
NUMBERED_ENVS = STATEMENT_ENVS | {
    "equation", "align", "gather", "multline", "figure", "table", "algorithm",
}

BEGIN_RE = re.compile(r"\\begin\{([A-Za-z*]+)\}")
END_RE = re.compile(r"\\end\{([A-Za-z*]+)\}")
LABEL_RE = re.compile(r"\\label\{([^}]*)\}")
REF_RE = re.compile(r"\\(?:[Cc]ref|crefrange|ref|eqref|autoref)\*?\{([^}]*)\}")
SECTION_RE = re.compile(r"\\(?:sub)*section\*?\{")
PROOF_OF_RE = re.compile(r"\\begin\{proof\}\[([^\]]*)\]")
LABEL_IN_TITLE_RE = re.compile(r"\\(?:[Cc]ref|ref|eqref|autoref)\*?\{([^}]*)\}")


def _strip_comments(line: str) -> str:
    """Drop everything after an unescaped %. Commented LaTeX is not a citation."""
    out = []
    for i, ch in enumerate(line):
        if ch == "%" and (i == 0 or line[i - 1] != "\\"):
            break
        out.append(ch)
    return "".join(out)


class Graph:
    def __init__(self) -> None:
        self.labels: dict[str, dict] = {}          # label -> {env, file, line}
        self.cites: dict[str, set[str]] = defaultdict(set)   # citer -> cited
        self.dangling: list[tuple[str, str, int]] = []       # label, file, line
        self.proofs: dict[str, str] = {}            # label -> "file:line" of its proof
        self.cite_sites: dict[str, list[dict]] = defaultdict(list)  # label -> where it is cited

    # -- construction -------------------------------------------------------
    def scan(self, sections: Path) -> "Graph":
        for path in sorted(sections.glob("*.tex")):
            self._scan_file(path)
        self._resolve()
        return self

    def _scan_file(self, path: Path) -> None:
        """Attribute every citation to its owning object.

        owner = the statement whose environment (or whose following `proof`) the
        citation sits in; otherwise the enclosing section. A citation in the prose
        AFTER a theorem belongs to the section, not to the theorem -- forgetting to
        reset on \\end{theorem} is what turns a section into a false dependency hub.
        """
        env_stack: list[str] = []
        current_section = f"file:{path.stem}"
        in_stmt: str | None = None
        last_stmt: str | None = None
        proof_owner: str | None = None

        for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            line = _strip_comments(raw)

            for env in BEGIN_RE.findall(line):
                env_stack.append(env.rstrip("*"))

            for label in LABEL_RE.findall(line):
                env = next((e for e in reversed(env_stack) if e in NUMBERED_ENVS), "section")
                self.labels[label] = {"env": env, "file": path.name, "line": lineno}
                if env in STATEMENT_ENVS:
                    in_stmt = label
                elif env == "section":
                    current_section = label

            titled = PROOF_OF_RE.search(line)
            if titled:
                # `\\begin{proof}[Proof of \\Cref{X}]` -- the proof belongs to X. The
                # label in the title is a heading, not a dependency, so strip the
                # bracket before scanning refs or the preceding lemma inherits it.
                named = LABEL_IN_TITLE_RE.search(titled.group(1))
                proof_of = named.group(1).split(",")[0].strip() if named else None
                line = line[:titled.start()] + line[titled.end():]
            else:
                proof_of = None
            if proof_of:
                proof_owner = proof_of
                self.proofs[proof_of] = f"{path.name}:{lineno}"

            if "proof" in env_stack:
                owner = proof_owner or in_stmt or last_stmt or current_section
            else:
                owner = in_stmt or current_section

            for group in REF_RE.findall(line):
                for label in (t.strip() for t in group.split(",")):
                    if label:
                        self.cites[owner].add(label)
                        self.cite_sites[label].append(
                            {"file": path.name, "line": lineno, "by": owner})
                        if label not in self.labels:
                            self.dangling.append((label, path.name, lineno))

            for env in END_RE.findall(line):
                env = env.rstrip("*")
                if env in env_stack:
                    for i in range(len(env_stack) - 1, -1, -1):
                        if env_stack[i] == env:
                            del env_stack[i]
                            break
                if env == "proof":
                    proof_owner = None
                if env in STATEMENT_ENVS and in_stmt:
                    last_stmt, in_stmt = in_stmt, None

    def _resolve(self) -> None:
        self.dangling = [d for d in self.dangling if d[0] not in self.labels]
        self.affected_by: dict[str, set[str]] = defaultdict(set)
        for citer, cited_set in self.cites.items():
            for cited in cited_set:
                if cited in self.labels:
                    self.affected_by[cited].add(citer)

    # -- queries ------------------------------------------------------------
    def affects(self, label: str, depth: int = 2) -> list[list[str]]:
        """Paths of NUMBERED STATEMENTS that would need re-checking if `label` changed.

        Statements only. Section and appendix labels are containers, not results:
        letting a path hop through one makes every statement in an appendix look
        like a dependant of every other.
        """
        def stmt_deps(node: str) -> list[str]:
            return sorted(n for n in self.affected_by.get(node, ())
                          if self.labels.get(n, {}).get("env") in STATEMENT_ENVS)

        paths: list[list[str]] = []
        frontier = [[label]]
        for _ in range(depth):
            nxt = []
            for path in frontier:
                for nb in stmt_deps(path[-1]):
                    if nb in path:
                        continue
                    paths.append(path + [nb])
                    nxt.append(path + [nb])
            frontier = nxt
        longest: dict[str, list[str]] = {}
        for p in paths:
            if len(p) > len(longest.get(p[-1], [])):
                longest[p[-1]] = p
        return [longest[k] for k in sorted(longest, key=lambda k: (-len(longest[k]), k))]

    def prose_sites(self, label: str) -> list[str]:
        """Sections and appendices whose prose cites `label` -- locations, not impact."""
        return sorted(n for n in self.affected_by.get(label, ())
                      if self.labels.get(n, {}).get("env") not in STATEMENT_ENVS)

    def statement_tex(self, label: str, sections: Path | None = None) -> str | None:
        """The LaTeX of the environment carrying `label`, for rendering in the app."""
        meta = self.labels.get(label)
        if not meta:
            return None
        base = sections or SECTIONS
        path = base / meta["file"]
        if not path.is_file():
            return None
        lines = path.read_text(encoding="utf-8").splitlines()
        env = meta["env"]
        start = next((i for i in range(meta["line"] - 1, -1, -1)
                      if f"\\begin{{{env}}}" in lines[i]), None)
        if start is None:
            return None
        end = next((i for i in range(start, len(lines))
                    if f"\\end{{{env}}}" in lines[i]), len(lines) - 1)
        return "\n".join(lines[start:end + 1])

    def statements(self) -> dict[str, dict]:
        return {k: v for k, v in self.labels.items() if v["env"] in STATEMENT_ENVS}


INPUT_RE = re.compile(r"\\input\{sections/([\w-]+)\}")
HEADING_RE = re.compile(r"\\section\*?\{(.+?)\}")


def document(sections: Path | None = None, master: Path | None = None) -> list[dict]:
    """The manuscript as an ordered list of parts -- the spine the map is drawn from.

    Order comes from the master's \\input sequence, not from sorting filenames: the
    appendices sort after the sections by luck, not by rule.
    """
    base = sections or SECTIONS
    mst = master or (base.parent / "thesis_v11.tex")
    order = INPUT_RE.findall(mst.read_text(encoding="utf-8")) if mst.is_file() else []
    if not order:
        order = [p.stem for p in sorted(base.glob("*.tex"))]
    parts = []
    in_appendix = False
    n_sec = 0
    for stem in order:
        path = base / f"{stem}.tex"
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        head = HEADING_RE.search(text)
        title = head.group(1).replace("---", "\u2014").strip() if head else stem.split("-", 1)[-1]
        if stem[0].isalpha() and stem[0].isupper():
            in_appendix = True
            mark = stem[0]
        else:
            mark = "" if stem.startswith("0") else str(stem[0])
            n_sec += 1
        parts.append({"file": path.name, "stem": stem, "title": title.title()
                      if title.isupper() else title,
                      "lines": len(text.splitlines()),
                      "mark": mark, "appendix": in_appendix})
    return parts


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sections", type=Path, default=SECTIONS)
    ap.add_argument("--check", action="store_true", help="print the integrity report")
    ap.add_argument("--affects", metavar="LABEL", help="show what one label affects")
    ap.add_argument("--json", action="store_true", help="dump the whole graph as JSON")
    ap.add_argument("--document", action="store_true", help="print the section spine")
    args = ap.parse_args(argv)

    if not args.sections.is_dir():
        print(f"no sections directory at {args.sections}", file=sys.stderr)
        return 2
    g = Graph().scan(args.sections)

    if args.affects:
        for path in g.affects(args.affects):
            print(" -> ".join(path))
        return 0
    if args.document:
        for part in document(args.sections):
            print(f"  {part['mark'] or '-':>2}  {part['lines']:>4} lines  {part['title']}")
        return 0
    if args.json:
        print(json.dumps({
            "labels": g.labels,
            "affects": {k: g.affects(k) for k in g.statements()},
            "dangling": g.dangling,
        }, indent=2, sort_keys=True))
        return 0

    stmts = g.statements()
    # Two different facts, and conflating them is how a cited result reads as dead:
    # no statement depends on it (so changing it breaks no proof) vs nothing mentions
    # it at all (so it may be deletable).
    no_dependants = [k for k in stmts if not g.affects(k, depth=1)]
    prose_only = sorted(k for k in no_dependants if g.prose_sites(k))
    rest = [k for k in no_dependants if not g.prose_sites(k)]
    proved_unused = sorted(k for k in rest if k in g.proofs)
    unreferenced = sorted(k for k in rest if k not in g.proofs)
    print(f"labels: {len(g.labels)}   statements: {len(stmts)}   "
          f"citing objects: {len(g.cites)}")
    print(f"dangling references: {len(g.dangling)}")
    for label, fname, line in g.dangling:
        print(f"  DANGLING {label}  {fname}:{line}")
    print(f"no statement depends on them, but prose cites them "
          f"({len(prose_only)}): {', '.join(prose_only) or 'none'}")
    print(f"proved but no statement or prose uses them ({len(proved_unused)}): "
          f"{', '.join(proved_unused) or 'none'}")
    print(f"nothing references them and no proof found ({len(unreferenced)}): "
          f"{', '.join(unreferenced) or 'none'}")
    if args.check:
        deepest = max(((len(p), p) for k in stmts for p in g.affects(k)),
                      default=(0, []))
        print("deepest impact chain:", " -> ".join(deepest[1]) or "none")
    return 1 if g.dangling else 0


if __name__ == "__main__":
    raise SystemExit(main())
