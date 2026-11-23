### [STY/01] Nineteen em-dash asides, replaced with commas, colons or parentheses
- **Anchor:** sections/1-introduction.tex:21
- **Source:** latex-paper-en:deai
- **Category:** punctuation
- **Type:** phrasing
- **Severity:** editorial
- **Context:** Em-dashes used to set off asides, spread across the introduction, the proof
  outline and Appendix A.
- **Item:** The manuscript uses 22 em-dashes in live text where the de-AI check caps a
  document at five, and the project style note bans the aside-in-dashes. The patch rewrites
  19 of them, choosing per site: a comma for an appositive, a colon where a restatement or
  list follows, parentheses for a true aside.
- **Affects:** none
- **Patch:** review/patches/STY-01.json
- **Status:** applied
- **Coverage:** 19 of the 22; the other three sit in AD/01 and AE/01, whose patches remove
  them, so no two patches touch the same text
- **Verified:** yes -- every site read in place before choosing its replacement, and no
  replacement reintroduces a dash
- **Applied:** 2026-11-23
- **Decision:** approved -- approved in chat
