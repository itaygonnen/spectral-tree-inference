### [STY/03] The introduction says "two things" and names neither
- **Anchor:** sections/1-introduction.tex:39
- **Source:** latex-paper-en:deai
- **Category:** ai-trace
- **Type:** phrasing
- **Severity:** editorial
- **Context:** The sentence opening the introduction's account of why the problem is hard,
  just before the two obstacles are described.
- **Item:** It opens "Two things keep this from being routine", where things and this point
  at nothing in particular, and the following clauses describe the obstacles anyway. Naming
  them directly costs no words.
- **Affects:** none
- **Patch:** review/patches/STY-03.json
- **Status:** applied
- **Coverage:** the only vague-referent hit in the document; four of the checker's other
  five flags were checked and found false, and are recorded in TOOL_NOTES.md
- **Verified:** yes -- read at 1-introduction.tex:39-42; the replacement keeps the citation
  and both clauses intact
- **Applied:** 2026-11-23
- **Decision:** approved -- approved in chat
