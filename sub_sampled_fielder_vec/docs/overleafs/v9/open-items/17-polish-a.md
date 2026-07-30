## POLISHA — Proof Outline / distance sibling, cut to author-approved size

Fragment owner: POLISHA. Anchors are `v9` labels in `sections/proof-outline.tex` and
`sections/main-result-dist.tex` unless a prior register ID (`R1/xx`, `B1/xx`) is given
explicitly. Scope: the author's review verdict — keep §Proof Outline but cut it to
~1/4 size, delete the flowchart, trim the distance-route sibling's exposition — not a
math change.

### [POLISHA/01] `lem:neumann`/`lem:pernode` kept in `proof-outline.tex`, not moved to the appendix
- **Anchor:** `\label{lem:neumann}`, `\label{lem:pernode}` (`sections/proof-outline.tex`)
- **Type:** direction
- **Item:** Per the brief's specified implementation choice, both boxed lemma statements stay physically in `proof-outline.tex` rather than relocating into `appendix-C.tex`. This is zero-risk: `appendix-C.tex` and `appendix-G.tex` already cite these labels by name via `\Cref`, not by file, and neither file was touched.
- **Blocking:** no — implementation choice already decided by the brief; recorded here for the label map's sake.

### [POLISHA/02] `rem:gap-scope`/`rem:dist-tradeoff` merged into one boxed remark
- **Anchor:** `\label{rem:gap-scope}\label{rem:dist-tradeoff}` (`sections/main-result-dist.tex:59`)
- **Type:** direction
- **Item:** The two donor remarks now share a single `\begin{remark}` environment carrying both labels, since `appendix-G.tex` cites `rem:dist-tradeoff` and the theorem's own proof-pointer sentence cites `rem:gap-scope`. Both claims' substance is kept (balanced-only closed form; worse constants vs.\ better edge choice are different axes), but the donor's "see the register entry for the alternative and its cost" pointer and its "neither substitutes for the other" closing line were cut as redundant with the adjacent `[B1/03]` citation and the preceding sentence.
- **Blocking:** no — took brevity over the redundant closing clauses; no claim was dropped.

### [POLISHA/03] Both files land above their line-count targets
- **Anchor:** `sections/proof-outline.tex` (138 lines, target ~100-110), `sections/main-result-dist.tex` (74 lines, target ~35-45)
- **Type:** flaw
- **Item:** In both files the content that must stay unchanged already exceeds the low end of the target before any prose is counted: `proof-outline.tex`'s two full boxed lemmas plus the untouched §6.4/§6.5 total ~66 lines, and `main-result-dist.tex`'s unchanged theorem, proof-pointer sentence, and macro block total ~27 lines. Compressing the remaining prose further than 3-6 (resp.\ 2-3) sentences per paragraph would breach the brief's own register floor.
- **Blocking:** no — compressed every paragraph to the mandated sentence floor and stopped there; flagging the arithmetic for the author rather than cutting substance to force the count.
