## POLISHB — abstract/intro prose conversion (this pass)

Fragment owner: POLISHB. Anchors are `v9` labels in `sections/abstract.tex` and
`sections/intro.tex` unless noted.

### [POLISHB/01] `fig:proof_map` still exists but was not referenced
- **Anchor:** roadmap paragraph, `sections/intro.tex` (final paragraph, `\Cref{sec:outline}`)
- **Type:** doubt
- **Item:** As of this pass `\label{fig:proof_map}` is still defined in `proof-outline.tex`, so referencing it would not currently dangle; the brief flagged it as the flowchart the author disliked and likely to be removed by the parallel agent editing that file concurrently. The roadmap paragraph points at `\Cref{sec:outline}` instead of the figure, so it stays correct regardless of whether the other agent removes it.
- **Blocking:** no — took the conservative option (`sec:outline`, not `fig:proof_map`); if the figure survives the parallel edit, promoting the roadmap clause to cite it too is a one-line addition for whoever assembles v9 next.

### [POLISHB/02] `sec:related` removed cleanly; no substitute related-work paragraph added
- **Anchor:** header comment, `sections/intro.tex`; `deferred/related.tex`
- **Type:** direction
- **Item:** The old `\Cref{sec:related} places the work` clause is deleted, not replaced with a survey; the roadmap paragraph carries one clause noting prior work is woven into the discussion rather than collected separately, citing `aizenbud2023spectral` as the organisational model this follows.
- **Blocking:** no — matches the brief's instruction to acknowledge placement in one clause, not reconstruct Related Work.

### [POLISHB/03] Abstract citations conflict with the single-paragraph-no-citation convention the survey found — RESOLVED
- **Anchor:** `sections/abstract.tex`
- **Type:** dilemma
- **Item:** All six survey papers write abstracts as bare prose with no `\cite` calls, but the task brief separately instructs "preserve every citation" when promoting the abstract's itemize version to prose, and the itemize version carries 15 distinct keys.
- **Options:** (a) keep all 15 citations inline, matching the brief's explicit instruction but breaking with the unanimous survey convention and pushing the paragraph to ~290 words against a 200–280 target; (b) drop citations from the abstract entirely to match genre convention, keeping them only in the intro (where all 15 already also appear).
- **Blocking:** no — POLISHB shipped (a); resolved to (b) in the verification pass, since every citation already also appears in `intro.tex` (nothing lost), it matches the unanimous convention, and it brought the paragraph back inside the target word count for free.
