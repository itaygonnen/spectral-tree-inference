# v9 refine pass — editorial restructure notes

Narrative/structure pass over `thesis_v9.tex`. **No mathematics was evaluated or
changed**: every `equation`, `theorem`, `lemma`, `proposition`, `corollary`,
`assumption`, `remark` and `proof` body is byte-identical to the pre-pass text
unless listed below as *unwrapped* (packaging removed, contents preserved).
Math-audit material stays in `open-items/`.

The same notes appear inline in the `.tex` sources as `% [Reviewer Note: ...]`,
invisible in the compiled PDF.

---

## Abstract — `sections/abstract.tex`

Rewritten as **problem → approach → bound → takeaway**. The previous version
spent five of eleven sentences on method exposition (the Neumann-vs-Davis–Kahan
contrast, the recovered $\sqrt m$, the two failure modes, the coherence identity
as a separate beat) — that is the body of the paper, not its advertisement. Each
survives as a subordinate clause. The closing "we prove no matching lower bound"
sentence is kept verbatim: it is a claim-discipline commitment
(`open-items/12-A3.md`).

## §1 Introduction — `sections/intro.tex`

**Dissolved:** the five numbered difficulties `(i)`–`(v)` and the eight-bullet
contributions list (≈160 lines) → two paragraphs plus a compressed roadmap.

Rationale: both lists landed exactly where the reader had just been promised a
question ("how small can $p$ be") and was ready for an answer. A list also
flattens dependency — difficulties (i), (ii) and (iv) are one compound obstacle
(coherence under imbalance) but read as three unrelated items. Now:

- **¶1** carries that compound obstacle as one escalating argument, ending on the
  two-mechanisms distinction (gap collapse → uninformative; coherence inflation →
  unaffordable).
- **¶2** carries the contributions in the order the paper delivers them.
  Difficulty (iii) moved here: it is not an obstacle to the paper but the reason
  for its method. Difficulty (v) is a scope limitation and is now stated as one,
  at the end.
- **Roadmap** compressed from a sentence-per-section walkthrough to four
  sentences; it retains the "Related Work is woven in, not sectioned" note, which
  does real work now that Related Work is deferred out of the document.

All 15 citation keys from the deleted lists are preserved.

## §2 Problem Definition — `sections/problem.tex`

**Dissolved:** subsections 2.1–2.5 (`sec:topology`, `sec:generative`,
`sec:blockmodel`, `sec:params`, `sec:problemdef`) → one continuous narrative,
former headings carried by transition sentences. Five headings over the setup
announced structure the reader could already see, and each restarted the prose.

**Unwrapped:** the three `\begin{definition}` environments — topological
imbalance, structural margin, subspace coherence. Displayed equations
`eq:imbalance_ratio`, `eq:margin_definition`, `eq:coherence_definition` are
unchanged and keep their labels; only the declarative packaging is gone, so
$\eta$, $\rho$ and $\mu(U)$ now arrive where they are first used.

**Received:** the construction of the centred distance operator
$\mathcal B = H\mathcal D H$, moved here from the old §4.2. It is setup, not
result — introducing a second operator inside the results section forced that
section to re-open the modelling discussion just after the main theorem had
closed it.

**Not touched:** all four `\begin{assumption}` environments and the assumption
chain table. `asm:margin` and its siblings are cited as hypotheses *inside*
`thm:main-sim` and must remain referenceable.

## §3 — `sections/bridging.tex`

**Retitled** "Bridging Tree Topology to Linear Algebra" →
**"Topological Scaling and Spectral Gap"**: the old title named the activity
rather than the content.

**Dissolved:** the three `\subsection*` headings ("Imbalance raises the subspace
coherence", "Imbalance erodes the spectral gap", "Two failure modes") and the two
`\paragraph` headings under the last of them. Each announced in a sentence what
the proposition or lemma immediately below stated formally, so the reader was
told everything twice. Proposition 1 and Lemma 1 are now stated directly and the
two failure modes run as continuous prose.

## §4 Main Result — `sections/main-result.tex`, `sections/main-result-dist.tex`

**Reordered to put Theorem 1 and Theorem 2 back to back.** Previously the section
stated Theorem 1, spent a corollary and a remark on the *reach* of that theorem,
and only then — in a subsection of its own — introduced a second operator from
scratch before stating Theorem 2; the paper's second main result arrived after
the section had begun winding down. Now:

1. lead-in (entry-wise criterion → a condition on $p$),
2. `thm:main-sim`,
3. one-sentence sibling-not-corollary bridge,
4. `thm:main-dist`,
5. shared four-step proof pointer to §5,
6. coda on reach: `cor:infeasible`, `rem:infeasible-mech`, `rem:gap-scope`.

`\subsection{The Distance Route}` and the `\paragraph{Two operators, one
criterion.}` heading are gone; the `\providecommand` macro block is promoted to
the preamble of `thesis_v9.tex` (kept as a fallback in the section file).

## §5 Proof Outline — `sections/proof-outline.tex`

**Dissolved:** subsections 5.1–5.5. A proof sketch is where momentum matters
most — the reader is being carried through an argument, not consulting a manual —
and five headings across two pages announced each step twice, once in the roadmap
sentence and once as a title. `lem:neumann`, `lem:pernode` and the proof
environment (with its internal Step 1–4 structure) are untouched.

## §6 Empirical Results — `sections/empirical.tex`

**Integrated:** §6.1 "Scoring the recovered split" into the experimental setup. A
metric is apparatus, not a result: a numbered subsection put a definition between
the reader and the first experiment and forced a forward reference ("§6.1 fixes
the metric before either sweep is reported") that the merge makes unnecessary.
The NMI definition, `cor:nmi` with its proof, and the measured-vs-population
caveat (three paragraphs → one) now sit in the setup; the section then flows
straight into the synthesized and generated regimes, which keep their headings
because they are results.

---

## Cross-reference repairs (mechanical, forced by the above)

| Was | Now | Where |
|---|---|---|
| `\Cref{def:imbalance}` | `\eqref{eq:imbalance_ratio}` | `intro`, `bridging`, `main-result` |
| `\Cref{def:margin}` | `\eqref{eq:margin_definition}` | `intro`, `main-result`, `deferred/related` |
| `\Cref{def:coherence}` | `\eqref{eq:coherence_definition}` | `bridging`, `appendix-A`, `appendix-B` |
| `\Cref{sec:blockmodel}` | `\Cref{sec:problem}` | `problem`, `appendix-G` |
| `\Cref{sec:gap-erodes}` | `\Cref{sec:bridging}` | `problem` (×2) |
| `\Cref{sec:problemdef}` | `\Cref{sec:problem}` | `main-result` |
| `\Cref{sec:failure-modes}` | `\Cref{sec:bridging}` | `intro`, `main-result-dist` |
| `\Cref{sec:outline-synth}` | `\Cref{sec:outline}` | `appendix-C` (×2) |
| `\Cref{sec:outline-gap/-neumann/-bern/-dist}` | prose ("Step 1", "above") | `proof-outline`, `main-result-dist` |

Two near-duplications introduced by the moves were also removed: the "recovery is
an entry-wise event" sentence no longer opens §4 (it closes §2), and the
"hypotheses are not nested" clause appears once, in §4 where Theorem 2 sits.

## Build state after the pass

```
latexmk -pdf thesis_v9.tex   →  exit 0
undefined references/citations: 0   (baseline: 0)
bibitems: 22, all still cited        (bibtex warnings: 0)
pages: 31                            (baseline: 33)
```
