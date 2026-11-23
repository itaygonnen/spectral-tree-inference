# What the review tools get wrong on this manuscript

Findings that were checked against the source and turned out to be the tool's
defect, not the manuscript's. They are recorded here instead of in the register so
the register stays a list of real defects -- and so nobody re-reports them.

## `latex-paper-en` / `check_references.py` -- `\Cref{a,b}` reads as one label

Reports `[Severity: Critical] [Priority: P0] Undefined reference` for every
cleveref comma list. Six of them on v11:

| reported at | the list |
|---|---|
| `6-topological-scaling.tex:15` | `asm:gtr,asm:clock` |
| `7-proof-outline.tex:11` | `app:comp1,app:comp2,app:comp3` |
| `7-proof-outline.tex:92` | `lem:A-decomp,lem:gap` |
| `C-entrywise-perturbation.tex:101` | `lem:C-inc,prop:C-within` |
| `C-entrywise-perturbation.tex:126` | `lem:pernode,prop:C-bind` |
| `E-split-choice.tex:6` | `asm:gtr,asm:clock` |

All eleven labels exist, and `thesis_v11.log` reports **0** undefined references.
`scripts/review/deps.py` splits on the comma, which is why it agrees with the build.

## `latex-paper-en` / `verify_bib.py` -- some braced surnames read as a missing author

Reports `[ERROR] Missing required field 'author'` for `blum2006which` and
`erdos1999logs`. Both have an `author` field. What those two share and the other
five braced-author entries do not is a brace nested inside a brace group within a
surname: `Fran{\c{c}}ois`, `Erd{\H{o}}s`. The exact parser rule was not pinned
down; the finding is false either way.

## The lesson, which is now a register field

`Verified:` is required on every item, and `register.py lint` fails an item whose
`Verified:` starts with "no". Two of the three highest-severity findings the tools
produced on the first run were their own bugs. An unchecked finding reaching the
author costs more than the finding is worth.

## `latex-paper-en` / `deai_check.py` -- appendices are never scanned

Run over the whole document it reports traces only for the abstract, the introduction and
the empirical section. Run directly on any of the five appendix files it reports zero and
prints an empty priority ranking. It buckets text by recognised section names, and an
appendix matches none, so the content is never examined. Anything the appendices get has
to come from reading them.

## `deai_check.py` -- four of six flags on this manuscript are false

| flagged | verdict |
|---|---|
| `0-abstract.tex:8` binary-contrast shell | **real**, and the whole abstract is worse than the one line (ABS/01) |
| `1-introduction.tex:39` vague referent | **real** (STY/03) |
| `1-introduction.tex:25` "a sequence-pair estimate ... rather than a lookup" | false: a real technical contrast naming a real alternative |
| `1-introduction.tex:74` "coordinate-wise rather than in Euclidean norm" | false: same |
| `1-introduction.tex:84` over-confident "always a split" | false: `thm:stdr-split` holds unconditionally, so "always" is accurate |
| `8-empirical-study.tex:23` binary-contrast shell | false: an orienting sentence describing a shared figure layout, with no contrast in it |

## `analyze_sentences.py` -- its rewrites destroy citations

It detects long sentences correctly, then "fixes" them by replacing commas with full
stops, which turns `\cite{a, b, c}` into `\cite{a. b. c}`. Detection is usable; the
suggested replacements are not, and none was taken.
