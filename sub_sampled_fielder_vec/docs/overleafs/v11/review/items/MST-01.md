### [MST/01] Two deleted shorthand commands came back from Overleaf and are unused
- **Anchor:** thesis_v11.tex:49
- **Source:** manual
- **Category:** dead-macro
- **Type:** notation
- **Severity:** material
- **Context:** The list of shorthand commands at the top of the master file, which the
  sections use to keep the notation consistent.
- **Item:** The Overleaf copy put back two commands that had been deleted, together with
  comments pointing at files that no longer exist. Neither command is used anywhere in the
  paper.
- **Affects:** none
- **Patch:** review/patches/MST-01.json
- **Status:** applied
- **Coverage:** `\Sout` and `\Dmat` have zero uses in any section file, so removing the two
  definitions changes no rendered output
- **Verified:** yes -- defined at thesis_v11.tex:49 and :74; `grep -rc` over sections/
  returns zero uses
- **Applied:** 2026-11-23
- **Decision:** approved -- approved in chat
