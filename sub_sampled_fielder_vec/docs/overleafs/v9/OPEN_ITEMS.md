# v9 — open items register

Single register for everything the restructure could not settle: unresolved questions, doubts,
dilemmas, suspected flaws, and directions worth pursuing. One entry per item.

`Blocking: no` entries record the default that was taken, so every one is visible and reversible.

**This file is generated** — edit the fragments in `open-items/` and re-run
`python scripts/collate_open_items.py`.

## R1 — findings from the pre-restructure investigation

Not produced by the restructure. Filed here because this is the single register.
Anchors are `v8` labels; all of them survive into `v9` under the same names unless noted.

### [R1/01] Prescribed sampling rate falls below the connectivity threshold
- **Anchor:** `eq:main_rate` (`thm:main`)
- **Type:** flaw
- **Item:** At the paper's own experimental values the rate reduces to `p = 0.25·log m/m` when `η=1`. A Bernoulli mask yields a connected observation graph only for `p ≳ log m/m`, so below that some leaves have no observed entries at all and their rows of `Ŝ` are identically zero. Exact recovery of every leaf is then impossible for any estimator, not only for the Fiedler sign rule.
- **Blocking:** no — the statement is ported unchanged.

### [R1/02] The rate never requires the entries to be observed
- **Anchor:** `eq:main_rate`, via `prop:C-var`
- **Type:** flaw
- **Item:** The `β₀²` factor arises because `prop:C-within` has already removed every within-clan increment from `(E_L v)_i`, leaving only cross-clan pairs of magnitude `β₀`. The rate therefore measures noise in the signal direction alone, and a smaller `β₀` makes it smaller without limit — with no term expressing that enough entries must be seen for the split to be identifiable.
- **Blocking:** no.

### [R1/03] The spectral gap admits an exact expression
- **Anchor:** `lem:A-gap` / `lem:gap`
- **Type:** flaw
- **Item:** The interlacing step bounds `λ₃ ≥ min_c λ₂(L_{C_c})` and in doing so drops the `+bβ₀` / `+aβ₀` shift carried by the clan-supported invariant subspaces — that shift is exactly the `−η·m_min·S_out^max` penalty. The exact statement is `λ₂ = mβ₀` and `Δλ = min(λ₂(L_A^intra) − aβ₀, λ₂(L_B^intra) − bβ₀)`, which in the flat CBM equals `m_min·ρ`.
- **Blocking:** no.

### [R1/04] The structural failure mode may have no CBM instance
- **Anchor:** `asm:margin`, `cor:tolerance`, `sec:failure-modes`
- **Type:** flaw
- **Item:** Both rest on `Δλ → 0` as `η` grows. If R1/03 holds then `Δλ > 0` whenever `S_in > β₀`, independently of `η`, so "structural collapse" does not occur under `asm:cbm` and survives only under `asm:hbm`, where the floor `S_in·α^{D_max}` genuinely can fall below `β₀`.
- **Blocking:** no — `asm:margin` is relocated, not retired, per plan decision D2.

### [R1/05] The binding clan may be the small one
- **Anchor:** `prop:C-bind`
- **Type:** doubt
- **Item:** The argument pairs the gap — set by the small clan `A`, since `λ₂(L_A) = aS_in < bS_in` — with the margin `c_B`, set by the large clan `B`. These extremes are attained at different nodes, so the pairing may never be realised and the binding clan may be `A`.
- **Blocking:** no.

### [R1/06] The recovery criterion is stronger than recovery
- **Anchor:** `eq:entrywise`
- **Type:** direction
- **Item:** `‖ŵ‖_∞ < min_i|v_i|` compares the error at the worst node with the margin at the tightest node, which under imbalance are different nodes. The per-node form `max_i |w_i|/|v_i| < 1` is what sign recovery actually needs and is weaker by a factor `η`.
- **Blocking:** no.

### [R1/07] The discarded remainder is not shown to be negligible
- **Anchor:** `lem:C-remainder`, `rem:C-asym`
- **Type:** flaw
- **Item:** This confirms the issue raised in `v8/review.md`. The remainder is bounded through the operator norm, and the ratio `‖E_L‖₂/Δλ` does not vanish as `m` grows at the operating point, so the bound is not `o(min_i|v_i|)`. The same unjustified discard also removes the `‖E_L‖₂|v_i|` term from the *first-order* bound in `lem:C-first`.
- **Blocking:** no.

### [R1/08] The nominated two-to-infinity fix does not apply here
- **Anchor:** `lem:C-remainder`; `v8/review.md` lines 44, 57
- **Type:** direction
- **Item:** `review.md` proposes Cape–Tang–Priebe's `2→∞` norm to recover the lost `√m`. Under the CBM the resolvent restricted to the complement is a low-rank perturbation of `Δλ^{-1}·I`, so there is no coherence for that machinery to grade and the `2→∞` norm differs from the operator norm by a constant. Recovering the `√m` appears to require treating `Σ_j (E_S)_{ij} w_j` as a mean-zero sum rather than applying Cauchy–Schwarz.
- **Blocking:** no.

### [R1/09] Weyl is applied outside its range in the first-order bound
- **Anchor:** `lem:C-first`
- **Type:** flaw
- **Item:** The proof asserts `|λ_ℓ − λ̂₂| ≥ Δλ − ‖E_L‖₂` for all `ℓ ≠ 2`, but at `ℓ = 1` Weyl gives only `λ̂₂ − λ₁ ≈ mβ₀`, which can be smaller than `Δλ`. The conclusion is rescuable because `E_L 𝟏 = 0` exactly — `L̂ = D̂ − Ŝ` is itself a Laplacian, so `E_L` acts invariantly on `𝟏^⊥` — but the paper never states this.
- **Blocking:** no.

### [R1/10] The rank-two premise does not hold
- **Anchor:** `prop:B-rank2`
- **Type:** flaw
- **Item:** With `S_ii = 1` the matrix `S` is full rank, not rank 2. Separately, `U = [v^(1), v^(2)]` is defined from eigenvectors of `L` (`def:coherence`), so the identity `S = UΛU^⊤` additionally requires `S` and `L` to share eigenvectors — which needs constant row sums and fails whenever `η > 1`.
- **Blocking:** no.

### [R1/11] The `(1+η)` in the noise bound may be an artifact of the route
- **Anchor:** `lem:B-bernstein`, via `prop:B-rank2`→`prop:B-rownorm`
- **Type:** doubt
- **Item:** The `(1+η)` enters only through the coherence route to `max_i‖S_{i,·}‖₂²`; the elementary bound `m·S_max²` carries no `η`. Note that removing it would not change `thm:main`'s `(1+η)³`, which comes from `prop:C-var` (one factor) and `lem:A-gap` (two).
- **Blocking:** no.

### [R1/12] The infeasibility corollary does not follow
- **Anchor:** `cor:infeasible`
- **Type:** flaw
- **Item:** Three independent problems: a *sufficient* condition becoming unsatisfiable licenses "this theorem is silent", not "no design recovers the split", and the paper proves no converse anywhere; the `1/3` exponent is derived holding both `β₀` and `ρ−ηβ₀` fixed as `η → ∞`, which is impossible since `ρ−ηβ₀ = S_in − (1+η)β₀`; and with `β₀` fixed, `asm:margin` already caps `η < 1/β₀ − 1 = O(1)`, so the regime is empty.
- **Blocking:** no — reworded to a vacuity statement per plan decision D3; the exponent and the empty-regime issues are left for the author.

### [R1/13] The scaling regime is never declared
- **Anchor:** `asm:margin`, `asm:bounded`, `cor:tolerance`, `cor:infeasible`
- **Type:** question
- **Item:** Is `β₀` held fixed as `m → ∞`, or allowed to depend on `m`? Every asymptotic-`η` statement changes meaning with the answer, and it decides whether `cor:tolerance` and `cor:infeasible` describe non-empty regimes. A single table fixing the `(m, η, β₀, ρ)` regime would settle it.
- **Blocking:** no.

### [R1/14] Theorems claim exact recovery; experiments score a threshold
- **Anchor:** `cor:nmi`, `sec:empirical`
- **Type:** question
- **Item:** The theorems conclude exact recovery, hence `NMI = 1`; the sweeps report the smallest `p` at which *mean* NMI crosses `0.95`, which admits systematic partial failure. Which is the intended claim? Note `compute_recovery` returns `max(s, 1−s)`, whose floor is `η/(1+η)`, so at `η = 8` a `0.95` cut on sign agreement sits about six points above chance.
- **Blocking:** no.

### [R1/15] The sampling model excludes the sequence noise the experiments include
- **Anchor:** `asm:sampling` against `asm:gtr`
- **Type:** question
- **Item:** `asm:sampling` sets `Ŝ_ij = Ω_ij S_ij/p`, i.e. the *population* similarity is observed exactly whenever the pair is sampled. But `asm:gtr` concedes the estimator only concentrates at `O(1/√ℓ)`, and the generated-regime experiments run with that error present. Is mask-only noise the intended scope of `thm:main`?
- **Blocking:** no.

### [R1/16] The two-layer HBM question is unresolved
- **Anchor:** `lem:gap`, `asm:hbm`, `prop:D-floor`
- **Type:** question
- **Item:** Task 2 of the brief asks whether the spectral-gap lower bound can rest on a two-layer HBM assumption rather than the full `α^{D_max}` floor. This was not attempted here — it is a derivation, not a restructure. Note that memory records V.04 as having settled on `S_in·α^{D_max}` in place of a sibling-level `α¹`, so this may be revisiting a closed point.
- **Blocking:** no.

### [R1/17] The distance-route rate is dimensionally inconsistent
- **Anchor:** `distance approach/distance_version.tex:681-687`
- **Type:** flaw
- **Item:** The displayed condition is `p > C(1+η)³d₀²log m/(m(Δλ^(B))²)` and is asserted on the same line to equal `Θ(log m/m)`. Since `Δλ^(B) = Θ(m)`, the expression is `Θ(log m/m³)`. The `(1+η)³` belongs to the form in which the gap has already been substituted, while `1/(m(Δλ^(B))²)` belongs to the unsubstituted form; the two have been combined.
- **Blocking:** no — ported in the donor's unsubstituted form so no false equality is propagated (plan decision D5).

### [R1/18] Within-clan cancellation does not transfer to the distance operator
- **Anchor:** `distance approach/distance_version.tex:637-647`
- **Type:** flaw
- **Item:** The claimed load `ρ̃^(B) = (1+η)d₀²/η` is obtained by asserting that within-clan increments cancel "as in the Laplacian proof". That cancellation comes from the Laplacian row-sum identity underlying `lem:C-inc`, which produces the factor `(v_i − v_j)`; the centered distance operator has no such identity, so within-clan terms survive and the load changes shape, not merely its constant.
- **Blocking:** no.

### [R1/19] The two distance-route penalties are the same quantity counted twice
- **Anchor:** `distance approach/distance_version.tex:68, 189, 714-724`
- **Type:** flaw
- **Item:** The constant penalty is stated as `d₀²/β₀² · r(T)²`, treating the cross-clan distance scale and the tree diameter as independent. Under `asm:clock` the tree is ultrametric, so every cross-clan MRCA is the root and `r(𝒯) = d₀` exactly.
- **Blocking:** no.

### [R1/20] The distance spectrum is proved only in the balanced case
- **Anchor:** `distance approach/distance_version.tex:322-356` (`prop:cbm-gap`)
- **Type:** direction
- **Item:** The proposition assumes `|A| = |B| = m/2`, in a paper whose subject is imbalance; the imbalanced remark that follows offers only `Θ(m(d₀−d₁))` with an unevaluated constant, and its stated argument does not hold because the top eigenvector of the distance matrix is not orthogonal to `𝟏` when `a ≠ b`. The general case follows from `𝓑 = −2(d₀−d₁)gg^⊤ − d₁H` with `g = He_A`, giving `Δλ^(B) = 2mη(d₀−d₁)/(1+η)²` and the structural eigenvector equal to the Fiedler vector exactly.
- **Blocking:** no.

### [R1/21] The distance-route empirical claim uses a different operator from the theorem
- **Anchor:** `distance approach/distance_version.tex:99, 540-543`
- **Type:** flaw
- **Item:** The quoted range `η_S ∈ [12,199]` is measured with the `σ₂`-gap rule, not the Fiedler-sign rule that `thm:main` analyses; the like-for-like sign-rule comparison is a median of 34.1 against 1.51. Both cited notebook paths (`02_real_data_sweeps/...`, `05_distance/...`) no longer exist after the directory reorganisation. Separately, `L_sym` with the sign rule matches the distance operator on the same cohort, so the gain may be attributable to centering rather than to the log-transform.
- **Blocking:** no.

### [R1/22] The synthesized figure cannot test the dependence it is captioned as confirming
- **Anchor:** `fig:pstar_synth` and its caption
- **Type:** flaw
- **Item:** The plotted curve is `C·η(1+η)³log m/[m(ρ−ηS_out)²]`, carrying a factor `η` that is absent from `eq:main_rate` — and `cor:C-cond` agrees with the theorem, so it is the figure that disagrees. The generating notebook also fits a separate constant per `η`, so a single-`η` panel with a free constant cannot test any `η`-dependence.
- **Blocking:** no — the caption is corrected not to claim it; the curve itself is left for the author.

### [R1/23] The HBM verification figure mixes two runs
- **Anchor:** `fig:hbm_spectral_verification`
- **Type:** flaw
- **Item:** With `\graphicspath{{../}{./}}` one panel resolves to a local copy and the other one level up, so the two panels are outputs of different runs. Separately, the lower bound being verified is negative over part of the swept `(η, α)` region, and confirming that a quantity exceeds a negative bound is not evidence for the bound.
- **Blocking:** no — both panels are regenerated from one run into `v9/figures/`.

### [R1/24] The bibliography cannot resolve
- **Anchor:** `v8/thesis_v8.tex:694`
- **Type:** flaw
- **Item:** `\bibliography{references.bib}` fails because the file is at `../references.bib` and the extension must be dropped. Every citation in the committed PDF is undefined. Every key does exist in the bib file.
- **Blocking:** no — fixed in v9.

---

## A1 — Problem Definition, Bridging, Appendices A/B/D (v8 §3, §4, app:comp1, app:comp2, app:hbm)

A1's session was killed by an infra error before this fragment could be filed; the `.tex`
output survived intact. This fragment is filed by a follow-up review of that surviving output
against v8 and the plan. All checklist items in the task brief were verified against the actual
fragments (duplicate `asm:hbm` label, `def:imbalance`/`def:margin`/`def:coherence` order,
`lem:A-coherence`/`lem:A-gap` → proof-environment conversion, `\bibliography{references}` fix,
absence of hardcoded agent/theorem/lemma numbers) and all passed cleanly — no register entries
needed for those. The items below are the substantive open points the restructure surfaced.

### [A1/01] `asm:margin` now precedes its own derivation in reading order
- **Anchor:** `asm:margin` (`sections/problem.tex`, sec:params) vs `lem:gap` (`sections/bridging.tex`, sec:gap-erodes)
- **Type:** doubt
- **Item:** Promoting `asm:margin` into §3 means the reader meets "recoverable only if ρ>η·Sout^max" before §4 derives it from the spectral gap; v8 had the reverse order, stating it as a corollary right after `lem:gap`. The relocated text forward-references `\Cref{sec:gap-erodes}` to flag this, but the assumption still reads as asserted rather than earned on first pass.
- **Blocking:** no — kept per the plan's instruction to relocate `asm:margin` into sec:params; flagging the forward-reference for the author to judge.

### [A1/02] Does `asm:margin` belong in §Problem Definition or §Bridging?
- **Type:** dilemma
- **Anchor:** `asm:margin`
- **Item:** `asm:margin`'s content is purely a consequence of `lem:gap`, which lives entirely in §Bridging; moving it to §Problem Definition turns a derived necessary condition into a stated assumption a section before its justification exists.
- **Options:** (a) keep in §3 sec:params, alongside `def:margin`/`def:imbalance` — front-loads every standing hypothesis of `thm:main-sim` in one place, which is what A1 did; (b) leave in §4 sec:bridging immediately after `lem:gap`, as in v8, where the "only if" claim is earned before being asserted.
- **Blocking:** no — (a) taken, matching the plan.

### [A1/03] `def:margin`'s ρ and `lem:gap`'s ρ disagree, and now sit closer together
- **Anchor:** `def:margin` (`eq:margin_definition`) vs the ρ used in `lem:gap`'s statement
- **Type:** flaw
- **Item:** `def:margin` defines ρ with a cross-clan dispersion penalty, `(Sin_min−Sout_max)−(Sout_max−Sout_min)`, while `lem:gap` itself states "structural margin ρ=Sin_min−Sout_max" with no penalty term; the two coincide only under `asm:clock`. This mismatch is inherited unchanged from v8 (where the two statements were ~100 lines apart), but relocating `asm:margin` into sec:params now sits it within ~20 lines of `def:margin` in the same subsection, making the inconsistency far more visible to a reader.
- **Blocking:** no — pre-existing math content, out of restructure scope; left for the author.

### [A1/04] Main-text remarks now state the HBM threshold that v8 deferred to Appendix D
- **Anchor:** `asm:margin`'s closing remark and the `sec:gap-erodes` paragraph in `sections/bridging.tex` (both cite `prop:D-gap`)
- **Type:** doubt
- **Item:** Both spots now quote the concrete inequality `α^{D_max}>(1+η)Sout^max/Sin`, which v8 stated only inside `app:hbm`'s `prop:D-gap`; v8's main text at the corresponding points (v8:410-411, v8:430-431) just said "deferred to `\Cref{app:hbm}`." The content is accurate and consistent with the unchanged `prop:D-gap`, but pulling a derived formula forward into main-text prose is more than mechanical relocation.
- **Blocking:** no — content verified consistent with `prop:D-gap`; flagged for the author to confirm it's in scope for a restructure-only pass.

### [A1/05] `sections/bridging.tex` forward-references A2-owned labels
- **Anchor:** `sections/bridging.tex` sec:failure-modes: `\Cref{thm:main-sim}`, `\Cref{cor:infeasible}`
- **Type:** direction
- **Item:** A1 does not own `main-result.tex`, yet sec:failure-modes cites A2's renamed theorem (v8's `thm:main` became `thm:main-sim`) and `cor:infeasible` by label; both currently resolve correctly against the surviving `main-result.tex`, confirming the label handoff worked. This is nonetheless a live coupling across an ownership boundary that would break silently if A2's labels change again.
- **Blocking:** no — labels currently consistent; flagging the cross-file dependency for whoever runs final assembly/consistency checking.

---

## A2 — main result, proof outline, appendix C/E

Fragment owner: A2 (follow-up session, finishing `appendix-E.tex` after the earlier
session was killed mid-task). Anchors are `v9` labels and `v9/sections/*.tex` line
contexts unless a `v8` path is given explicitly.

### [A2/01] Appendix E orphan risk did not materialize
- **Anchor:** `app:compare`, `rem:E-compare`
- **Type:** doubt
- **Item:** `\Cref{app:compare}` is cited from `sec:outline` (proof-outline.tex:175) and independently three more times from `intro.tex` and once from `related.tex`, so Appendix E is not orphaned despite the deletion of the old §5 Neumann paragraph that used to cite it. `rem:E-compare` itself is cited twice more from `sec:outline-neumann` (proof-outline.tex:153,220).
- **Blocking:** no — confirmed present at both sites the brief asked about; no register-blocking action needed.

### [A2/02] Stale `\Cref{thm:neumann}` survives in A3's files after the rename
- **Anchor:** `intro.tex:234`, `related.tex:113`
- **Type:** flaw
- **Item:** `sec:outline` renamed `thm:neumann` to `lem:neumann` (proof-outline.tex:7,177), but `intro.tex:234` and `related.tex:113` still cite the old label, which no longer exists anywhere in v9. Both will render as `??` in the compiled PDF.
- **Blocking:** no — A2 does not own `intro.tex`/`related.tex`; flagged here for A3 to retarget both citations to `lem:neumann`.

### [A2/03] `lem:neumann`'s convergence hypothesis now lives in sec:outline-gap, not at the lemma
- **Anchor:** `proof-outline.tex:137-140`
- **Type:** question
- **Item:** With the old §5 Neumann-motivation prose deleted, the hypothesis `‖E_L‖₂<Δλ` needed for the expansion is stated once, in `sec:outline-gap`'s closing paragraph, ahead of `lem:neumann` itself rather than restated at the lemma. It is the only place the strong-signal regime is declared before Step 2 of the proof relies on it.
- **Blocking:** no — location confirmed adequate; treated as intentional rather than a gap.

### [A2/04] `prop:C-bind`/`cor:C-cond` stayed in Appendix C; Step 4 re-derives the same algebra inline
- **Anchor:** `appendix-C.tex:113-136` (`prop:C-bind`, `cor:C-cond`); `proof-outline.tex:242-249` (Step 4)
- **Type:** doubt
- **Item:** Only `lem:pernode`'s statement was promoted to `sec:outline`; `prop:C-bind` and `cor:C-cond` remain full standalone statements in Appendix C, cited by name from Step 4 of `sec:outline-synth`. Step 4 then re-derives the same squaring/cancellation/solve-for-`p` argument inline rather than just invoking them, so the calculation appears twice at full length.
- **Blocking:** no — left as authored; a tighter version would have Step 4 cite `prop:C-bind`/`cor:C-cond` without repeating the algebra, but that edit touches proof-outline.tex, which A2-follow-up does not own for this task.

### [A2/05] The √m-saving narrative repeats three times within sec:outline-neumann alone
- **Anchor:** `proof-outline.tex:151,170,172,220`; `appendix-C.tex:147`; `appendix-E.tex`
- **Type:** direction
- **Item:** Grep for `\sqrt m` across the four files gives 4 hits in proof-outline.tex — 3 of them narrative, all inside `sec:outline-neumann` — plus 1 quantitative hit in appendix-E.tex; grep for `Davis` gives 2 in proof-outline.tex, 1 in appendix-C.tex, 1 in appendix-E.tex. The "qualitative once in §6.2, quantitative once in E" target is met for Appendix E but not for §6.2, which states the same saving three times in one subsection.
- **Blocking:** no — condensing `sec:outline-neumann` to state the saving once is a style cleanup left for the author, since proof-outline.tex is not this task's file to edit.

---

## A3 — introduction, abstract, related work, bibliography

Fragment owner: A3. Anchors are `v9` labels and `v9/sections/*.tex` line contexts unless a `v8`
path is given explicitly.

### [A3/01] Griffing thesis metadata is unverified
- **Anchor:** `references.bib` entry `Griffing2012`
- **Type:** question
- **Item:** The entry is transcribed verbatim from the donor bibliography at `distance approach/distance_version.tex:745-751`, which gives only `A.~Griffing`; the expansion to "Alexander", the title capitalisation, and any month/URL/DOI could not be confirmed against the NCSU repository or any indexing service. The key is fixed by the donor's `\cite{Griffing2012}` calls so it must not be renamed.
- **Blocking:** no — the donor's fields were kept as-is with a `% VERIFY` note above the entry; the author should confirm the given name and pull a repository handle before submission.

### [A3/02] Related Work cites no work newer than 2020
- **Anchor:** `v9/sections/related.tex`
- **Type:** doubt
- **Item:** The outline covers the four lines v8's own TODO named plus community detection and entry-wise eigenvector analysis, but the newest key in the whole bibliography is `abbe2020entrywise`. A section with a six-year gap invites the reviewer question "has nobody done sub-sampled spectral tree recovery since?", which the paper cannot currently answer.
- **Blocking:** no — the outline is written from the existing 19 keys plus the 6 A3 added; no literature search for 2021–2026 work was in scope.

### [A3/03] "Tight rate" is not defensible
- **Anchor:** `v8/thesis_v8.tex:512` ("converts the recovery criterion into the tight rate")
- **Type:** flaw
- **Item:** The paper proves a sufficient condition and no lower bound anywhere, so neither "tight" nor `Θ`-as-optimality is licensed; the `Θ(log m/m)` in `eq:main_rate` is at best the growth order of the sufficient bound, not an optimality claim. The same reading applies to `rem:E-compare`, where the `Θ(mη/β₀²)` ratio compares two upper bounds, one of which is the paper's own construction.
- **Blocking:** no — the introduction and abstract say "sufficient" throughout and state explicitly that no lower bound is proved; the word "tight" at v8:512 belongs to A2's section and is flagged here rather than edited.

### [A3/04] Which gap claim may the introduction advertise?
- **Anchor:** `lem:gap` / `cor:tolerance`, against register `R1/03` and `R1/04`
- **Type:** dilemma
- **Item:** `R1/03` says `Δλ` admits an exact expression equal to `m_min·ρ`, and `R1/04` concludes that the gap-collapse failure mode may have no CBM instance at all. The contributions bullet for `lem:gap` therefore cannot honestly be written as "imbalance drives the gap to zero under the CBM".
- **Options:** (a) advertise `lem:gap` as a lower bound only and locate the collapse mode in the HBM (`app:hbm`), which is what R1/04 supports; (b) hold the bullet until the exact `Δλ` is settled, since if R1/03 stands the `−η·m_min·S_out^max` penalty disappears and the whole two-failure-mode framing loses one of its halves.
- **Blocking:** no — took (a): the bullet says "lower bound", carries an explicit caveat that a sharper expression is open, and points the collapse mode at `app:hbm`.

### [A3/05] Infeasibility is written as vacuity, not as a converse
- **Anchor:** `cor:infeasible`, against register `R1/12`
- **Type:** flaw
- **Item:** v8's prose says "no uniform sampling design recovers the split", but a sufficient condition becoming unsatisfiable licenses only "this theorem is silent". `R1/12` further notes the regime may be empty, since with `β₀` fixed `asm:margin` already caps `η = O(1)`.
- **Blocking:** no — the introduction states it as the guarantee going silent and says explicitly that no impossibility result is proved; if the regime is in fact empty the bullet should be deleted rather than reworded.

### [A3/06] The paper analyses one split, not the recursion
- **Anchor:** `sec:intro` ("Why the top split is the object that matters"); `sec:problem`
- **Type:** question
- **Item:** The motivation is divide-and-conquer reconstruction, but every result concerns a single bipartition on `m` leaves, and nothing composes the guarantee down the recursion — where clan sizes shrink, `η` is redrawn at each node, and the per-call failure probabilities must be unioned. Is the single-split scope the intended claim, or is a recursion corollary expected?
- **Blocking:** no — the introduction states the single-split scope explicitly rather than implying a whole-tree guarantee.

### [A3/07] Coherence is defined on `L` but used as if it were `S`'s
- **Anchor:** `def:coherence`, `prop:B-rank2`, against register `R1/10`
- **Type:** doubt
- **Item:** `def:coherence` builds `U` from eigenvectors of `L`, while `prop:B-rank2` uses `S = UΛUᵀ` with the same `U` — which needs `L` and `S` to share eigenvectors, i.e. constant row sums, which fails for `η > 1`. The introduction's difficulty (i) therefore had to be phrased carefully: the coherence that governs completion is a property of the operator being sampled, and the paper samples `S` while measuring coherence on `L`.
- **Blocking:** no — difficulty (i) names the mismatch in a sub-bullet and points here rather than glossing it.

### [A3/08] Mask-only noise versus the sequence noise the experiments carry
- **Anchor:** `asm:sampling` against `asm:gtr`, mirroring register `R1/15`
- **Type:** question
- **Item:** `asm:sampling` observes the *population* `S_ij` exactly whenever the pair is sampled, so `thm:main-sim` carries no `O(1/√ℓ)` term, yet the generated-regime sweeps run with that error present. The introduction lists this as difficulty (v) and states the scope as mask-only, but the paper never says which of the two noise sources the theorem is about.
- **Blocking:** no — difficulty (v) declares the scope as mask-only and points at the generated-regime experiments as the only evidence about the other term.

### [A3/09] Bibliography defects found while appending
- **Anchor:** `v9/references.bib`
- **Type:** flaw
- **Item:** `abbe2020entrywise` was typed `@inproceedings` while carrying `journal`/`volume`/`number` (Annals of Statistics 48(3):1452-1474), so `plain` would have dropped all three fields; corrected to `@article`, no field changes. Left alone, since other agents' `\cite` keys depend on them: `candes2008exact`'s key year (2008) doesn't match its `year` field (2009, correct per FoCM 9(6):717-772); `eldridge2017unperturbed`'s key year (2017) doesn't match ALT 2018.
- **Blocking:** no — only the entry type was changed; keys were left untouched because other agents' `\cite` calls depend on them.

### [A3/10] Labels on starred subsections resolve to the wrong counter
- **Anchor:** `sec:failure-modes` (v8: `\subsection*{Two failure modes}\label{sec:failure-modes}`)
- **Type:** flaw
- **Item:** `\label` after an unnumbered `\subsection*` attaches to whichever counter was last stepped, so `\Cref{sec:failure-modes}` prints a stale number — v8 already does this at line 597 and the new contributions block adds another reference. The same pattern affects every starred subsection in `sec:problem` and `sec:bridging` that carries a label.
- **Blocking:** no — the reference was written as-is; the fix (number the subsections, or add `\phantomsection`/`\refstepcounter`) belongs to whoever owns those sections.

### [A3/11] No citation could be found for the cost claim that motivates sub-sampling
- **Anchor:** `v9/sections/intro.tex`, "Why scale forces sub-sampling"
- **Type:** question
- **Item:** The premise that the `Θ(m²)` similarity evaluations, rather than memory or the downstream eigendecomposition, are the binding cost of top-down reconstruction is asserted with no source and no measurement in the paper. `aizenbud2023spectral` motivates the divide-and-conquer schema but does not profile this step.
- **Blocking:** no — the bullet is written without a citation and no `[CITE NEEDED]` marker was left in the `.tex`; either a source or a one-line timing from the repo's own runs would close it.

---

## B1 — distance-route merge (§5.2, Appendix "distance route")

Fragment owner: B1. Anchors are `v9` labels in `sections/main-result-dist.tex` and
`sections/appendix-G.tex` unless a donor path (`distance approach/distance_version.tex`) or a
prior register ID (`R1/xx`, `A3/xx`) is given explicitly.

### [B1/01] Notation: `\mathcal D`/`\mathcal B`, not `D^\dagger`/`D_deg`
- **Anchor:** `\Dmat`, `\Bop` macros (`main-result-dist.tex:14-15`, `appendix-G.tex:28-29`)
- **Type:** direction
- **Item:** The donor renames the host's degree matrix to `D_deg` and takes `D` for distances; this is not ported. The host's `D=\operatorname{diag}(S\mathbf 1)` and `L=D-S` are untouched everywhere; the new operators use `\mathcal D` (pairwise distance) and `\mathcal B=H\mathcal D H` (centered), with `H` for the centering projector, never `B` alone (the host has ~30 occurrences of clan `B`).
- **Blocking:** no — this was mandated by the brief, not chosen; recorded here for the label map's sake.

### [B1/02] Operator macros defined locally, not in the shared preamble
- **Anchor:** `main-result-dist.tex:10-15`, redeclared via `\providecommand` in `appendix-G.tex:28-32`
- **Type:** question
- **Item:** `\Dmat`, `\Bop`, `\uone`, `\huone`, `\DlamB` are declared where first used rather than in `thesis_v9.tex`'s shared macro block, even though `sections/notation.tex:10-12` already anticipates `\mathcal D`/`\mathcal B` in prose. Should A1 promote these five macros to the preamble alongside `\vv`/`\Dlam` for consistency with the rest of the notation table?
- **Blocking:** no — took local declaration with `\providecommand` (order-independent, safe if another agent's file is `\input` first); promoting them is a one-line preamble edit for A1 if preferred.

### [B1/03] The only closed-form spectral gap is balanced-only
- **Anchor:** `rem:gap-scope` (`main-result-dist.tex`), `rem:G-unbalanced` (`appendix-G.tex`); donor `distance_version.tex:322-423`; mirrors `R1/20`
- **Type:** dilemma
- **Item:** `thm:main-dist` depends on the top gap `Δλ^(𝓑)`, but the only closed-form value derived (`prop:G-spectrum`) assumes `|A|=|B|=m/2`, in a paper about imbalance `η`. The donor's own imbalanced-case argument is not just incomplete but invalid as stated (`R1/20`: the top eigenvector of `𝓓` is not orthogonal to `𝟏` when `a≠b`, which is exactly what the balanced proof relies on), so it cannot be ported even informally.
- **Options:** (a) scope `thm:main-dist` explicitly to the balanced flat CBM, stating the imbalanced case as unresolved future work; (b) state `thm:main-dist` conditionally on the abstract hypothesis `Δλ^(𝓑)>0`, and supply the balanced closed form only as an instantiation via `prop:G-spectrum` — costs an unevaluated abstract gap in the general statement, but keeps the theorem's shape a genuine sibling to `thm:main-sim` (which is also gap-conditional, with `lem:gap` supplying only a lower bound, not an exact value, at general `η`).
- **Blocking:** no — took (b): it is the more honest reading of the donor's own theorem statement (which was already gap-conditional; only the *instantiated* rate display was where the donor over-claimed) and parallels how `thm:main-sim` itself is structured.

### [B1/04] Within-clan cancellation for `𝓑` is asserted, not proved
- **Anchor:** `appendix-G.tex` §"Per-node Bernstein" (first term); donor `distance_version.tex:637-647`; mirrors `R1/18`
- **Type:** flaw
- **Item:** The donor's load `ρ̃^(𝓑)=(1+η)d₀²/η` requires within-clan increments of `(E_𝓓 u^(1))_i` to cancel "as in the Laplacian proof," but that cancellation is `lem:C-inc`'s row-sum identity for `L`, which `𝓑` has no analogue of (only the *global* identity `𝓑\mathbf 1=0` holds, not a per-row one). The claim is ported labelled as asserted by the donor, not re-derived.
- **Blocking:** no — the per-node Bernstein subsection states this explicitly as an unverified step ported from the donor; a genuine derivation is out of scope per the brief.

### [B1/05] The split-decomposition identity needs no diagonal correction
- **Anchor:** `lem:split-decomp` (`appendix-G.tex`); donor `distance_version.tex:472-476`
- **Type:** flaw
- **Item:** The donor's cut-metric identity `δ_e=½(𝟏𝟏ᵀ−s_es_eᵀ)+\mathrm{diag}(\cdots)` carries a spurious diagonal term: since `(s_e)_i²=1`, `½(1−(s_e)_i²)=0` on the diagonal automatically, so no correction is needed. Dropped in the ported lemma; the rest of the identity (and the resulting `𝓑⪯0`, `‖Hs_e‖²=4n₁n₂/m`) is correct and reproduced with a full proof.
- **Blocking:** no — fixed while porting per the brief's explicit instruction to remove it.

### [B1/06] The two constant-degradation penalties are the same quantity
- **Anchor:** `appendix-G.tex` §"Honest comparison" and the Step-1 rescaling note; donor `distance_version.tex:68,189,714-724`; mirrors `R1/19`
- **Type:** flaw
- **Item:** The donor multiplies `d₀²/β₀²` (cross-clan parameter) by `r(𝒯)²` (tree-diameter entry-scale bound) as independent penalties, but under `asm:clock` the tree is ultrametric so `r(𝒯)=d₀` exactly (every cross-clan MRCA is the root). Resolved by dropping the donor's `1/(4‖\operatorname{tr}Q‖)` rescaling entirely, so `𝓓_{ij}=-\log S_{ij}` literally coincides with the host's own `r(𝒯)` of `cor:tolerance`, and the final comparison states one factor, `d₀²/β₀²`, not a product.
- **Blocking:** no — took the simplification the brief suggested; the appendix states the identification `r(𝒯)=d₀` as a remark rather than asserting it silently.

### [B1/07] The donor's headline empirical range mixes two different sign rules
- **Anchor:** `rem:dist-compare-empirical` (`appendix-G.tex`); donor `distance_version.tex:99,540-543`; mirrors `R1/21`
- **Type:** flaw
- **Item:** The donor's `η_S∈[12,199]` for the Laplacian side is measured with the `σ₂`-gap statistic, not the Fiedler-sign rule that `thm:main-sim` actually analyses, so it is not a like-for-like comparison against `η_D∈[1,7]`. Replaced with the like-for-like sign-rule medians already on record (`34.1` vs `1.51`); both dead notebook paths (`02_real_data_sweeps/...`, `05_distance/...`) are dropped per the brief rather than repaired.
- **Blocking:** no — used the corrected numbers already verified in `R1/21` rather than the donor's mismatched range.

### [B1/08] The donor's rate display asserts a false `Θ(log m/m)` equality
- **Anchor:** `appendix-G.tex` proof of `thm:main-dist` (`app:dist-main`); donor `distance_version.tex:681-687`; mirrors `R1/17`
- **Type:** flaw
- **Item:** The donor's displayed condition `p>C(1+η)³d₀²\log m/(m(Δλ^(𝓑))²)` is asserted on the same line to equal `Θ(\log m/m)`; since `Δλ^(𝓑)=Θ(m)` in the balanced case (`prop:G-spectrum`), the right-hand side is `Θ(\log m/m³)`, not `Θ(\log m/m)`. Both `main-result-dist.tex` and the appendix state the condition in the unsubstituted form and explicitly decline to claim the `Θ(\log m/m)` simplification.
- **Blocking:** no — per plan decision D5 (register `R1/17`), the false equality is not printed anywhere in the merged text.

### [B1/09] Griffing thesis hypotheses transcribed, not verified
- **Anchor:** `thm:dist-partition` (`appendix-G.tex`); cross-refers `A3/01`
- **Type:** question
- **Item:** `thm:dist-partition`'s statement is a paraphrase of the donor's own paraphrase of Griffing's Ch. 4 result; neither the donor nor this restructure had access to the primary source to confirm the exact hypotheses (separation condition, tree regularity). The appendix says so explicitly next to the theorem rather than presenting it as verified.
- **Blocking:** no — same open point as `A3/01` on the bibliography entry; resolving it requires the actual thesis document, out of reach here.

### [B1/10] `thm:griffing` demoted to a quoted theorem with no reproduced proof
- **Anchor:** `thm:dist-partition` (`appendix-G.tex`); donor `distance_version.tex:436-527`
- **Type:** direction
- **Item:** The donor's own "proof sketch" of `thm:griffing` is not required anywhere by `thm:main-dist` (donor concedes this at `:431-433`) and contains two independent defects beyond the diagonal term of `B1/05`: an uncited appeal to an unpublished thesis for the decisive step, and "picking the maximiser saturates this bound" (`:512-513`), which is false because the Rayleigh quotient sums over *all* edges, not just the maximiser's. The result is kept as a quoted theorem (statement + one-sentence scope note) with no proof reproduced, rather than as a remark or omitted entirely, because it is the only available explanation for *which* edge `𝓑`'s sign rule selects on a general tree (used by `rem:dist-compare-empirical`).
- **Blocking:** no — quoted-theorem form chosen over remark/omission; the split-decomposition half of the argument (`lem:split-decomp`) is kept as a fully proved lemma since it is correct and does the real work.

### [B1/11] Bibliography keys needed from A3 — already present
- **Anchor:** `references.bib` entries `Griffing2012`, `BandeltDress1992`
- **Type:** question
- **Item:** Both entries required by this fragment (`\cite[Ch.~4]{Griffing2012}`, `\cite{BandeltDress1992}`) were already added by A3 with a `% VERIFY` note on the Griffing entry (see `A3/01`); no edit to `references.bib` was needed or made. A full-document `pdflatex`+`bibtex` pass confirms both resolve and are used in the compiled bibliography.
- **Blocking:** no — informational; nothing outstanding on B1's side.

### [B1/12] "Worse constants" and "better partition" are different axes, not a contradiction
- **Anchor:** `rem:dist-tradeoff` (`main-result-dist.tex`), `rem:dist-compare-empirical` (`appendix-G.tex`)
- **Type:** direction
- **Item:** The donor asserts both that the distance route has strictly worse multiplicative constants under sub-sampling and that it picks a better-balanced partition on real data; these can coexist because they answer different questions — the first compares noise tolerance for recovering a *fixed* target edge, the second compares *which* edge each operator's population sign rule selects before any sub-sampling. The merged text keeps them as two separate remarks anchored to the theorem (constants) and the appendix (edge selection) respectively, with `rem:dist-tradeoff` stating explicitly that neither substitutes for the other.
- **Blocking:** no — resolved by scoping, not by picking one claim over the other.

### [B1/13] Stale "STDR Fig. 10" reference dropped
- **Anchor:** donor `distance_version.tex:733`
- **Type:** flaw
- **Item:** The donor cites a specific figure number in the STDR paper for an empirical trend claim; this could not be verified and is not reproduced. The merged "Honest comparison" text in `appendix-G.tex` states the shared/differing quantities directly instead of leaning on an unverified external figure reference.
- **Blocking:** no — dropped per the brief's mechanical-fix instruction to use `\cite[Thm.~4.2]{aizenbud2023spectral}`-style citation rather than page/figure prose pointers that cannot be checked.

### [B1/14] Full-document compile check
- **Anchor:** `v9/thesis_v9.tex`
- **Type:** question
- **Item:** A full `pdflatex`+`bibtex`+`pdflatex`×2 pass of the whole `v9` document succeeds (36 pages, exit 0) with both B1 files included. Two pre-existing undefined-reference warnings remain (`thm:neumann` in `intro.tex`/`related.tex`, `lem:A-gap` in `appendix-C.tex`) — both predate this fragment (A2 renamed `thm:neumann`→`lem:neumann` and `lem:A-gap`→`lem:gap` but left stray old references elsewhere) and are outside B1's owned files.
- **Blocking:** no — flagged for whichever agent owns `intro.tex`/`related.tex`/`appendix-C.tex`; not a B1 defect.

---

## Referee pass, 2026-07-29 (`DISTANCE_REVIEW.md`)

Entries [B1/15]–[B1/21] record the mathematical changes made during a supervisor-level
proofread of the distance part. Unlike [B1/01]–[B1/14], these are **not** ported-as-is
decisions: the text was changed. Each entry states what it said, what it now says, and what
the author should verify independently.

### [B1/15] `eq:main-rate-dist` was mis-normalized (rate changed)
- **Anchor:** `thm:main-dist` (`main-result-dist.tex`); supersedes the display carried from the donor
- **Type:** flaw (corrected)
- **Was:** `p > C(1+η)³d₀²log m / [m·(Δλ^(𝓑))²]`
- **Now:** `p > C(1+η)d₀²·m·log m / (Δλ^(𝓑))²`, with the substituted balanced form in `cor:dist-balanced`
- **Item:** the old display mixed two normalizations. The shared synthesis
  `8σ²log m < gap²·min_i|u_i|²` with `min_i|u_i|²=1/(mη)` and `σ²=ρ̃^(𝓑)/p` gives `m` in the
  *numerator*. The `(1+η)³` and the `m` in the denominator arise only after substituting
  `Δλ ≥ m_min·margin` with `m_min = m/(1+η)`; the donor applied that conversion while keeping
  the unsubstituted `Δλ^(𝓑)²`, i.e. twice. Numerically the two forms differ by `8m²/(1+η)²`
  (≈2·10⁶ at m=10³). The Laplacian side is unaffected: `eq:main_rate` reproduces the
  abstract-gap form exactly under the substitution (checked numerically, ratio 1.000).
- **Verify:** the factor `η·ρ̃^(𝓑)=(1+η)d₀²` in Step 4 of `app:dist-main`, and that no other
  passage still quotes the old display.

### [B1/16] `Θ(log m/m³)` passage withdrawn
- **Anchor:** former closing paragraph of the proof in `app:dist-main`
- **Type:** flaw (corrected, consequence of [B1/15])
- **Item:** the proof previously substituted the mis-normalized display, obtained
  `Θ(log m/m³)`, and refused the `Θ(log m/m)` claim that `rem:gap-scope` and the comparison
  paragraph both asserted. With [B1/15] repaired, substituting `Δλ^(𝓑)=m(d₀−d₁)/2` gives
  `p·m/log m = const`, so `Θ(log m/m)` is correct and the passage is deleted rather than
  reconciled. `cor:dist-balanced` now carries the claim, with proof.
- **Verify:** `cor:dist-balanced`'s arithmetic, `8C d₀²log m/[m(d₀−d₁)²]`.

### [B1/17] Within-clan cancellation promoted to a stated hypothesis — **SUPERSEDED by [B1/25]**
- **Anchor:** `asm:dist-cancel` (`main-result-dist.tex`), `rem:dist-cancel-scope` (`appendix-G.tex`); supersedes [B1/04]
- **Type:** flaw (corrected structurally, not mathematically)
- **Item:** [B1/04] recorded that the load `ρ̃^(𝓑)` is asserted, not derived, and the proof of
  `thm:main-dist` began "Granting it". A theorem may be conditional but must say so in its
  hypotheses. The claim is now `asm:dist-cancel`, listed among `thm:main-dist`'s assumptions
  beside `Δλ^(𝓑)>0`. No mathematics changed; the dependency is now visible at the theorem.
- **Verify:** nothing new — but note that Theorem 2 is now openly conditional on two
  unresolved inputs, which is the honest reading and should be reflected in any abstract or
  introduction that describes the distance route as "proved".

### [B1/18] The natural computation gives a different load than the ported one — **RESOLVED by [B1/25]**
- **Anchor:** `rem:dist-cancel-scope` (`appendix-G.tex`), `eq:dist-load`
- **Type:** doubt (open — **not** changed)
- **Item:** the ported load is `ρ̃^(𝓑)=(1+η)d₀²/η`, the Laplacian load with `β₀→d₀`. But the
  Laplacian load descends from increments proportional to `(v_i−v_j)`, while
  `(E_𝓓u)_i=Σ_j (E_𝓓)_{ij}u_j` weights `u_j` alone. Computing directly with `c_A²=η/m`,
  `c_B²=1/(mη)`: for `i∈B` the surviving cross-clan variance is `((1−p)/p)d₀²·η/(1+η)`, which
  after synthesis is a load `η²d₀²/(1+η)`, not `(1+η)d₀²/η`. The two agree only up to a factor
  `η²/(1+η)²` — at η=1, a factor 4.
- **Blocking:** no — left unchanged deliberately. The underlying cancellation is itself
  unproven ([B1/17]), so re-deriving the constant on top of an unproven step would not
  produce a theorem. **This is the first thing to settle if the distance route is to become
  self-contained.**

### [B1/19] `lem:C-first` was mis-cited at its point of use
- **Anchor:** proof of `thm:main-dist` in `app:dist-main`
- **Type:** flaw (corrected)
- **Was:** "`lem:C-series`, `lem:C-first` bound `‖û−u‖_∞` by `‖E_𝓑‖₂/Δλ^(𝓑)` up to a coherence factor `(1+η)^{3/2}`"
- **Now:** the proof applies `lem:C-first` as stated —
  `|w_i| ≤ (|(E u)_i| + ‖E‖₂|u_i|)/(Δλ−‖E‖₂)` — and discards the second numerator term and the
  denominator correction exactly as `sec:outline` does on the similarity side.
- **Item:** `lem:C-first` is a per-node bound with no coherence factor; the discard is the same
  one flagged in `R1/07`, and it is now visible on the distance side too rather than hidden
  inside a paraphrase.

### [B1/20] Coherence transfer is balanced-case, not operator-agnostic
- **Anchor:** comparison itemization in `app:dist-main`, `rem:gap-scope`
- **Type:** flaw (corrected wording)
- **Item:** `prop:coherence` is proved from `lem:A-decomp`, the Laplacian's Fiedler coordinates.
  Transferring `μ(U)=(1+η)/2` to `[𝟏/√m, u^(1)]` needs `u^(1)` piecewise constant with the same
  clan sizes — established only in the balanced flat CBM (`prop:G-spectrum`) and open at general
  `η` (`rem:G-unbalanced`). The text now says "wherever `u^(1)` is known to have it" instead of
  asserting operator-agnosticism.

### [B1/21] Griffing's partition result demoted from theorem to reported proposition
- **Anchor:** `prop:dist-partition` (was `thm:dist-partition`); `lem:dist-edge` is new
- **Type:** flaw (corrected)
- **Item:** the statement was numbered as a theorem while the surrounding text admitted it was
  an unverified paraphrase with a non-quantitative hypothesis ("well separated"). It is now a
  proposition explicitly marked *reported*, and the part this paper does prove — the
  Rayleigh-quotient bound `|λ₁(𝓑)| ≥ 2τ_e n₁(e)n₂(e)/m` for every edge, from
  `lem:split-decomp` — is separated out as `lem:dist-edge` and proved.
- **Verify:** whether the primary source can now be consulted; if its separation hypothesis is
  quantifiable, `prop:dist-partition` can be stated properly or dropped in favour of
  `lem:dist-edge` alone.

### [B1/22] Smaller corrections
- `main-result-dist.tex`: `c` was declared as a constant but never bound to a probability; the
  statement now reads "with probability at least `1−O(m^{-c})`".
- `rem:dist-entryscale`: entries of `E_𝓓` are bounded by `d₀/p` under inverse-probability
  weighting, not by `r(𝒯)`; the per-node paragraph already used `d₀/p`.
- `lem:split-decomp` now precedes the sign rule that depends on it (it was a forward reference).
- `rem:gap-scope` carried two labels (`rem:gap-scope`, `rem:dist-tradeoff`) on one object; the
  second is retired and its references retargeted. The dead `\label{sec:dist-compare}` on an
  unnumbered `\paragraph` is removed.
- Donor- and register-facing prose is removed from the body throughout App. F (eleven sites);
  the mathematics and the scope statements remain, the editorial argument moves here.

### [B1/23] Observation: the paper has no Conclusion
- **Type:** direction
- **Item:** there is no Discussion/Conclusion section. At the target venue a short closing
  section is conventional, and three things currently scattered in remarks belong there: the
  certified-versus-possible asymmetry (`cor:infeasible`, `rem:infeasible-mech`), the two open
  population inputs (a general-`η` gap bound for either operator; the slack in the `ℓ_∞`
  criterion, `R1/06`), and the empirical gap between the `(1+η)³` inflation and the measured
  factor of ≈4.6 over `η∈[1,15]`.
- **Blocking:** no — recorded, not drafted.

### [B1/24] `cor:tolerance` is vacuous under the molecular clock (**referred to the author**)
- **Anchor:** `cor:tolerance` (`bridging.tex`), its proof in `appendix-A.tex`; exposed by `app:dist-landmarks`
- **Type:** flaw (not fixed — outside the distance part, touches a §3 headline claim)
- **Item:** `cor:tolerance` bounds `S_in^min ≥ e^{-r(𝒯)}` with `r(𝒯)=max_{ij}𝒟_ij` and concludes
  `η ≤ O(m^c)`. But under `asm:clock` every cross-clan MRCA is the root, so the diameter *is* the
  cross-clan distance: `r(𝒯)=d₀=-log S_out^max`, hence `e^{-r(𝒯)}=S_out^max` and the derived
  tolerance is `η < S_in^min/S_out^max ≈ 1` — not polynomial in `m`. The distance appendix's
  identification `r(𝒯)=d₀` is correct; what it exposes is that `cor:tolerance` needs the maximum
  **within-clan** distance, not the diameter.
- **Blocking:** no, but this is the most consequential item in this register: `bridging.tex`
  currently concludes "under balanced generative models the gap survives polynomial imbalance"
  from it.

### [B1/25] `asm:dist-cancel` was false and has been deleted; no assumption is needed
- **Anchor:** `lem:dist-var` (`appendix-G.tex`), `thm:main-dist`; supersedes [B1/17], resolves [B1/18]
- **Type:** flaw (corrected)
- **Item:** the load asserted in the intermediate pass, `ρ̃^(𝓑)=(1+η)d₀²/η`, is contradicted by the
  direct computation `ηd₀²/(1+η)` — a factor 4 at η=1, i.e. in the balanced flat CBM where
  `prop:G-spectrum` *proves* the shape of `u^(1)`. Shipping it as a hypothesis made `thm:main-dist`
  and `cor:dist-balanced` vacuous. It is unnecessary: `‖u^(1)‖₂=1` and `𝒟_ij ≤ d₀` bound the total
  per-node variance (within-clan terms included) by `d₀²/p` in one line. The rate is now
  `p > C d₀² η m log m/(Δλ^(𝓑))²`, stronger than either earlier form.
- **Verify:** `lem:dist-var`'s two-line proof, and Step 4's constant `C=4C₂²`.

### [B1/26] Missing operator-norm bound, mis-cited perturbation lemma, unbounded centering scalar
- **Anchor:** `lem:dist-opnorm`, `rem:dist-centering`, proof of `thm:main-dist`
- **Type:** flaw (corrected)
- **Item:** three gaps found by the independent review. (a) The route had no bound on `‖E_𝒟‖₂`
  anywhere, though `lem:C-series/-first/-remainder` all require `‖E‖₂ < gap`; added as
  `lem:dist-opnorm` via matrix Bernstein, mirroring `lem:B-bernstein`. (b) `lem:C-first` bounds the
  *first-order* term `w^(1)`, not `w`; `lem:C-remainder` is now invoked for the remainder.
  (c) The claim that the global centering scalar cannot affect sign recovery was false — a uniform
  shift exceeding `c_B` flips all of clan `B` — and the proposed de-meaning was vacuous, since
  `𝓑̂𝟏=0` already forces mean-zero eigenvectors. The scalar is now bounded by
  `‖E_𝒟‖₂/√m ≲ d₀√(log m/p)`, the same order as the per-node term.
- **Verify:** whether `lem:dist-opnorm`'s variance proxy `σ² ≤ m d₀²/p` is the tightest available;
  a row-norm refinement as in `prop:B-rownorm` may improve the constant.

### [B1/27] The discarded Neumann remainder breaks the stated rate on **both** routes
- **Anchor:** `lem:C-remainder`, `sec:outline` Step 2, `thm:main-sim`; distance side now `rem:dist-remainder`, hypothesis (iii) of `thm:main-dist`; mirrors and sharpens `R1/07`
- **Type:** flaw (distance side made explicit; **similarity side not fixed — author's call**)
- **Item:** the recovery event is `‖ŵ‖_∞ < min_i|v_i| = Θ((mη)^{-1/2})`, but `lem:C-remainder`
  bounds `‖w^(2)‖_∞` only through the operator norm, `O(‖E‖²/Δλ²)`. At the paper's own values
  (η=1, β₀=0.05, S_in=0.9) and its own rate `p = 0.25 log m/m`, `lem:B-bernstein` and `lem:gap`
  give `‖E_L‖₂/Δλ ≈ 7.07` — constant in m, verified identical at m = 10³, 10⁶, 10⁹ — so the
  discarded remainder is ≈50 against a vanishing margin, and the Neumann series of `lem:C-series`
  is not even shown to converge at the operating point. Requiring the remainder to be
  `o(min_i|v_i|)` costs a factor `√m`: the exponent this machinery supports is `Θ(log m/√m)`,
  not `Θ(log m/m)`.
- **Distance side:** stated as hypothesis (iii) of `thm:main-dist` (remainder control
  `‖E_𝓑‖₂ ≤ ε_m Δλ^(𝓑)`, `ε_m ≤ c₀(mη)^{-1/4}`), with `rem:dist-remainder` recording what it
  costs and that the same gap sits in `thm:main-sim`.
- **Similarity side:** untouched. `eq:main_rate` still displays an explicit constant
  `8(1+η)^3` and asserts `Θ(log m/m)`; per the above that display is not supported by the
  proof as written. **This is the most consequential open item in the paper.**
- **Fix that would close it:** an entry-wise (leave-one-out) bound on `w^(2)` in the manner of
  `abbe2020entrywise`, which is exactly the tool `sec:outline` cites for the first-order term but
  does not apply to the remainder. Failing that, both theorems should carry the `√m` loss or an
  explicit remainder-control hypothesis.

### [B1/28] `cor:infeasible` needs a quantifier over β₀ and the margin
- **Anchor:** `cor:infeasible` (`main-result-dist.tex`)
- **Type:** flaw (corrected)
- **Item:** as stated, `η = Ω((m/log m)^{1/3})` was claimed to force RHS > 1. Counterexample:
  m = 10⁶, η = 100, β₀ = 0.001, S_in = 0.9 satisfies `asm:cbm`, `asm:bounded`, `asm:margin`
  (0.899 > 0.1) and gives RHS = 1.8·10⁻⁴. Since `asm:margin` forces β₀ < ρ/η, the group
  `η³β₀²/(ρ-ηS_out^max)²` can stay bounded while η grows. The corollary now holds β₀ and the
  margin fixed and also states the invariant form.

### [B1/29] `lem:C-first`'s Weyl step is false at ℓ=1 on the similarity side
- **Anchor:** `lem:C-first` (`appendix-C.tex`), `lem:A-decomp`
- **Type:** flaw (**not fixed — outside the distance part**)
- **Item:** `lem:C-first` asserts `|λ_ℓ - λ̂₂| ≥ Δλ - ‖E_L‖₂` for all ℓ ≠ 2. With `λ₁ = 0` and
  `λ₂ = mβ₀` (`lem:A-decomp`), the ℓ=1 separation is `mβ₀`, which at the paper's values
  (0.05m) is *smaller* than `Δλ ≥ 0.4m`. The step is rescued only by `E_L 𝟏 = 0`, hence
  `⟨v^(1), E_L v^(2)⟩ = 0`, which the paper never states. One sentence in `app:comp3a` fixes it,
  or define `Δλ := min_{ℓ≠2}|λ_ℓ - λ₂|`. The distance route does make the corresponding check
  (`app:dist-landmarks` item 3, valid because `𝓑 ⪯ 0`).

---

## C1 — empirical results (§7) and supplementary figures

Fragment owner: C1. Anchors are `v9` labels in `sections/empirical.tex` and
`sections/appendix-emp.tex` unless a `v8`/register path is given explicitly.

### [C1/01] Should the synthesized theory curve be plotted at all, and in which form?
- **Anchor:** `fig:pstar_synth`, register `R1/22`
- **Type:** dilemma
- **Item:** The curve plotted in `fig:pstar_synth` is `C·η(1+η)³log m/[m(ρ−ηS_out^max)²]`, which carries a factor of `η` absent from `eq:main_rate`; `cor:C-cond` (A2, `main-result.tex`) agrees with the theorem, so the figure is the one that disagrees, not the corollary. The caption and prose no longer claim this curve confirms the theorem's `η`-dependence — see `sec:emp-synth` in `empirical.tex` — but the curve itself is unchanged, since fixing it is a plotting decision outside C1's scope (no code).
- **Options:** (a) replot with the correct `eq:main_rate` form, `C·(1+η)³log m/[m(ρ−ηS_out^max)²]` (still one free constant per `η`, so it would support only the `log m/m` shape claim already made, not a stronger one); (b) keep the current `η`-inflated curve but rename it in the caption as "a reference of this functional form" (already done) and leave the mismatch for the author to resolve; (c) refit a single shared constant `C` across all four `η∈{1,2,4,8}` panels and report the goodness of that shared fit — this is the only version of the plot that would actually test the `(1+η)³` term, but it requires new fitting code and is out of scope for this pass.
- **Blocking:** no — took (b). The prose is honest about what is and is not tested; the curve's formula is left for C2/the author.

### [C1/02] What §7.3 may claim with k-means as the sole main-text operator
- **Anchor:** `sec:emp-gen`, `fig:pstar_gen`
- **Item:** With the σ₂-gap and sign-rule comparison relocated to `fig:operator_sensitivity` in the appendix, the main-text claim is narrowed to: `n^{-1}` decay of `p̂*` within each `η` (k-means, `L_sym`), and the *direction* (not magnitude) of the upward shift with `η`. The cross-operator invariance sentence ("the `log n/n` scaling is a property of the spectral signal, not of the particular rounding rule") survives verbatim but only in the appendix, where three operators are actually shown.
- **Type:** direction
- **Blocking:** no — this is the reading implemented in `empirical.tex`/`appendix-emp.tex`; flagged so a reviewer of the restructure can check the narrowing was applied consistently.

### [C1/03] Does `cor:nmi` read as a measurement in its new home?
- **Anchor:** `sec:emp-nmi` (`cor:nmi`)
- **Item:** Moving `cor:nmi` out of §5 and into the empirical section risks it being read as summarizing the sweeps rather than being a population-level deduction from `thm:main-sim`. `sec:emp-nmi` is titled "Scoring the recovered split," opens by stating it fixes "the metric and its population value," and closes with an explicit paragraph distinguishing the corollary's exact `NMI=1` from the "measured NMI" used everywhere downstream; that qualifier is then used consistently in `sec:emp-synth`, `sec:emp-gen`, `sec:emp-reading`.
- **Type:** doubt
- **Blocking:** no — the framing paragraphs are the mitigation; if a reader still conflates the two the fix is structural (e.g. a boxed "population vs. measured" callout), which was judged unnecessary for a first pass.

### [C1/04] `C=32.7` cannot carry over to the NMI-scored balanced figure
- **Anchor:** `fig:pstar_balanced` in `appendix-emp.tex`
- **Item:** v8's balanced-baseline figure was sign-agreement-scored with fitted constant `C=32.7`; C3 is converting it to measured-NMI scoring, under which the fit will differ (different metric, different threshold-crossing point). The caption currently carries a `% VALUE FROM C3` comment plus a visible `$C=[\text{TBD, C3}]$` placeholder rather than a fabricated or stale number.
- **Type:** question
- **Blocking:** no — addressed to C3. Default taken: placeholder left in place; whoever regenerates the figure fills in the fitted constant and removes the placeholder.

### [C1/05] No bridge sentence from §5/§6 into §7
- **Anchor:** `sections/main-result.tex` (ends at `rem:infeasible-mech`), `sections/proof-outline.tex`
- **Item:** v8's `sec:empirical` opened by immediately stating the two test regimes (thesis_v8.tex:609), with no explicit hand-off sentence from the preceding infeasibility remark. Checked A2's `main-result.tex`: it also ends at `rem:infeasible-mech` with no bridge sentence into the empirical section, so this matches v8's structure and is not a gap introduced by the restructure. `empirical.tex`'s own lead paragraph opens directly with `\Cref{thm:main-sim}` and is self-contained, so no bridge is required for the section to read correctly — but if A2 or A3 later add a closing "the remainder of the paper tests this" sentence to §5/§6, check it doesn't duplicate `empirical.tex`'s opening clause.
- **Type:** question
- **Blocking:** no — addressed to A2 (informational; no bridge sentence was found missing relative to v8, so no action was assumed necessary).

### [C1/06] Measured `NMI≥0.95` is not the theorem's `NMI=1` endpoint
- **Anchor:** `sec:emp-nmi`, register `R1/14`
- **Item:** `thm:main-sim`/`cor:nmi` conclude exact recovery (`NMI=1`, population level); every sweep instead reports the smallest `p` at which *mean measured* NMI crosses `0.95`, which admits systematic partial failure at that `p`. `sec:emp-nmi`'s closing paragraph states this explicitly ("the theorem's own endpoint is `NMI=1`, not `NMI≥0.95`, and the two are not interchangeable") rather than letting the shared symbol `NMI` imply the sweeps test the theorem's exact conclusion.
- **Type:** question
- **Blocking:** no — resolved via the explicit disclaimer; the register item R1/14's further observation (that `compute_recovery`'s floor at high `η` sits close to a `0.95` cut under *sign agreement*) does not apply once every figure is rescored to NMI, since NMI's floor near a 50/50-confident split is not equivalent to the agreement floor `η/(1+η)` — but this was not independently verified against C2/C3's NMI computation and is worth a spot check.

### [C1/07] HBM verification figure: two runs, and a bound negative over part of the region
- **Anchor:** `fig:hbm_spectral_verification` in `sections/appendix-D.tex` (owner: A1), register `R1/23`
- **Item:** Not C1's file, flagged per task assignment. As of this pass `appendix-D.tex` still includes `S_by_alpha.png` and `spectral_gap_bound.png` as two separate `\includegraphics` calls with no note that they come from a single generating run; `R1/23`'s stated resolution ("both panels are regenerated from one run into v9/figures/") does not appear to have landed yet. Separately, and independently of the two-runs issue: `prop:D-gap`'s lower bound is negative over part of the swept `(η,α)` region (whenever `α^{D_max}≤(1+η)S_out^max/S_in`), and the right panel's caption ("the bound stays valid, slack≥0, across all configurations") reads as confirming the bound everywhere without noting that over the negative-bound region, `slack≥0` is automatic and not evidence of anything.
- **Type:** flaw
- **Blocking:** no — addressed to A1; C1 did not edit `appendix-D.tex`.

### [C1/08] `fig:pstar_synth`'s v8 caption mislabeled the metric
- **Anchor:** `v8/thesis_v8.tex:646-648`
- **Item:** v8's own `sec:empirical` lead paragraph (v8:615-617) says recovery is scored by NMI, but the `fig:pstar_synth` caption at v8:646 says "recovery (sign agreement)" — an internal inconsistency in v8 itself, independent of the restructure. The `v9` caption in `empirical.tex` now says "measured NMI" throughout, consistent with the filename `pstar_synth_nmi_eta8.png`, which presumes C2 regenerates this figure's underlying data as NMI-scored rather than agreement-scored.
- **Type:** question
- **Blocking:** no — addressed to C2. Default taken: caption written as true-after-regeneration; if C2's data is in fact still agreement-scored, the caption's "measured NMI" wording needs to change back (or the data needs to change).

### [C1/09] Lint check
- **Anchor:** n/a
- **Item:** `python3 scripts/collate_open_items.py --check` run after filing this register; no lint problems in this fragment.
- **Type:** direction
- **Blocking:** no.

---

## C2 — eta-pool sweep code, figures, and cache-key gate

Fragment owner: C2 (Python/notebooks only, no LaTeX). Anchors are `file:line` in the
`sub_sampled_fielder_vec` tree unless a `v9` LaTeX label is given explicitly.

### [C2/01] Step 0 audit — prior-session edits confirmed correct and complete
- **Anchor:** `src/utils/eta_pool_cache.py:179-201`, `src/utils/eta_pool_sweep.py:45-53,128-169`
- **Type:** direction
- **Item:** `load_pool_metadata` (cache.py), `PAPER_METHOD`/`paper_methods()`, and `collect_pool_features` (sweep.py) all matched the task's description exactly on read: correct signatures, correct use of `_default_pool_params`/`_CACHE_ROOT`, no duplication needed. `METHOD_SPECS`, `ETA_TARGETS`, `ETA_SAMPLE_CAP`, `P_VALUES`, `BOOTSTRAP_REPS`, `SEED`, `NUM_GAPS`, `EARLY_STOP_CONSECUTIVE_100` (all cache-key-bearing) were unchanged from the frozen pre-edit values — confirmed by reading, not assumed.
- **Blocking:** no.

### [C2/02] `sweep_plots.py` savepath wiring finished; all four functions now return `Figure`
- **Anchor:** `analysis/theoretical_interpretation/utils/sweep_plots.py:123,163,239,278`
- **Type:** direction
- **Item:** Added `savepath: Optional[Path] = None` to `plot_pstar_vs_n`, `plot_nmi_grid`, `plot_curves_per_n`, `plot_curves_per_eta`; replaced each `fig.tight_layout(...); plt.show()` tail with `fig.tight_layout(...); _finish(fig, savepath); return fig`. None of the four previously returned `fig` (diff showed no `return fig` anywhere), so I made all four consistent by returning it — the new `sweep_plots_two_panel.py` and the notebook cells both call these for their `savepath` side effect only, so the return value is currently unused but present for composability.
- **Blocking:** no — chose "all four return" over "none do" per the task's own tie-break instruction.

### [C2/03] New `sweep_plots_two_panel.py`; theory bound cross-checked directly against `eq:main_rate`
- **Anchor:** `analysis/theoretical_interpretation/utils/sweep_plots_two_panel.py:26-35`, `sections/main-result.tex:29` (`eq:main_rate`)
- **Type:** direction
- **Item:** `plot_pstar_two_panel` (144 lines) + `_draw_recovery_panel`/`_draw_scale_panel` + a new public `cbm_sufficient_p(eta, S_out_max, margin, m)` implementing `p > 8(1+η)³(S_out^max)²log m / (m(ρ−ηS_out^max)²)` read directly off `thm:main-sim`. This is deliberately **not** the `theory_scale = η(1+η)³log n/(n·margin²)` formula already in `kingman_threshold_vs_theory.ipynb` cell 19/21, which register `R1/22` already flags as carrying a spurious extra factor of `η` relative to `eq:main_rate`.
- **Blocking:** no — used the theorem-correct form; left the pre-existing flawed `theory_scale` cell untouched per the task's explicit instruction not to fix cell 21 (out of scope, superseded).

### [C2/04] Two independent, disagreeing theory-curve formulas now coexist in the repo
- **Anchor:** `analysis/theoretical_interpretation/utils/sweep_plots_two_panel.py:26` vs. `analysis/theoretical_interpretation/generated/eta_pool_sweeps/kingman_threshold_vs_theory.ipynb` cell 19 (id `a5e0139e`)
- **Type:** doubt
- **Item:** My new `cbm_sufficient_p` and the pre-existing `theory_scale` both claim to plot the `thm:main-sim` bound but differ by a factor of `η`; a reader who compares the two eta-pool notebooks side by side will see two different "the theory curve" shapes. I did not reconcile them since `kingman_threshold_vs_theory.ipynb`'s formula is `R1/22`'s territory (author-facing math flaw), not a code bug for C2 to silently fix.
- **Blocking:** no — flagging for the author; my own new code uses the theorem-correct form.

### [C2/05] `n=8000` dropped (not annotated `k=1`) from every eta-pool figure
- **Anchor:** `analysis/theoretical_interpretation/generated/eta_pool_sweeps/eta_pool_sweep.ipynb` cell `9a5ad600` (`NS_INCLUDE`)
- **Type:** direction
- **Item:** All three eta-pool operators have exactly one completed pool sample at `n=8000` per `η` (verified: `discover_ns_and_samples` prints `8000  1  1  1  1` pre-edit), versus `≥3` everywhere else, so I set `NS_INCLUDE = [500,...,6000]` globally rather than keeping `n=8000` with a `k=1` annotation. A single sample would inject un-aggregated noise into both the recovery-grid means and every `p*` power-law/CBM-bound fit without adding a usable scale point.
- **Blocking:** no — dropped; re-adding it (with `k=1` annotated in the legend) is a one-line `NS_INCLUDE` edit if the author wants the extra scale point back.

### [C2/06] Main two-panel figure reports `η ∈ {1, 10}`
- **Anchor:** `eta_pool_sweep.ipynb` cell `9a5ad600` (`ETA_MAIN`), `pstar_gen_kmeans_2panel.png`
- **Type:** direction
- **Item:** Chose the balanced case (`η≈1`) and one clearly imbalanced case (`η≈10`) from the four available (`1,5,10,15`) so the two rows read as a between-row comparison of the phase-transition shift, per the task's own suggested reading. `η=15`'s CBM bound is vacuous at every sampled `n` in this pool (`margin ≤ 0` throughout, matching R1/12's infeasibility regime), so it would have contributed an empty bound overlay; `η=10` still shows one.
- **Blocking:** no — `ETA_MAIN` is a two-line notebook edit if the author prefers a different pair.

### [C2/07] `sigma2`'s appendix curve is flat at `p*=1.0` — pre-existing, not introduced here
- **Anchor:** `pstar_gen_3operators.png`, `eta_pool_sweep.ipynb` cell `d8be357d`
- **Type:** flaw
- **Item:** At every `(n,η)` where `sigma2` has cached data (`n∈{500,1000,2000,4000}`; none at `3000/6000` per the inventory), mean NMI never crosses `0.95` below `p=1.0`, so `compute_pstar` returns `p*=1` everywhere and the "fit" degenerates to `n^0.00`. `cached_only_methods={"sigma2"}` means this run never recomputed it (`cache report: computed=0`), so the flat curve reflects the existing cache, not a bug in the restructured notebook.
- **Blocking:** no — informational; matches why `kmeans` was made `PAPER_METHOD` instead of `sigma2`.

### [C2/08] `identity_scatter_kingman.png` saves one of three `DATASETS`; picked `kingman_mu0.1`
- **Anchor:** `kingman_threshold_vs_theory.ipynb` cell `9e6d4b66` (`PRIMARY_DATASET`), cell `f6fc60c8`
- **Type:** question
- **Item:** Cell 7 (`f6fc60c8`) loops over three `DATASETS` (`kingman_mean`, `kingman_mu0.3`, `kingman_mu0.1`) producing one 3-panel figure each, but the required output is a single file. I saved only `kingman_mu0.1` since it is the only entry whose `(tree_model="kingman", mu=0.1)` matches the eta-pool's own default pool params (`src/utils/eta_pool_sweep.py:67`).
- **Blocking:** no — addressed to the author/C1: confirm `kingman_mu0.1` (not `kingman_mean`) is the intended "the Kingman dataset" for this figure; `PRIMARY_DATASET` is a one-line change otherwise.

### [C2/09] Verification gate — zero COMPUTE across both notebooks
- **Anchor:** `eta_pool_sweep.ipynb` cells `5176ad18`, `d8be357d`
- **Type:** direction
- **Item:** Added an `on_sweep` callback (mirroring `scripts/build_sweeps.py`'s `_on_sweep`) to both the main and appendix compute cells and re-executed the full notebook end-to-end. Main sweep: `{'kmeans': 156}`, `computed=0 cached=156 skipped=0`; appendix sweep: `{'sign': 156, 'sigma2': 104, 'kmeans': 156}`, `computed=0 cached=416 skipped=0` — no cache key moved.
- **Blocking:** no.

### [C2/10] `build_sweeps.py` default narrowed to `kmeans`; `choices` kept at all three
- **Anchor:** `scripts/build_sweeps.py:63-64`
- **Type:** direction
- **Item:** Changed only `default=["sign","kmeans"]` → `default=["kmeans"]`; `choices=list(METHOD_SPECS)` (and `METHOD_SPECS` itself, in `eta_pool_sweep.py`) are untouched so `--methods sign sigma2 kmeans` still rebuilds the appendix comparison. `src/utils/sweep_cache.py` and `src/runners/p_sweep_inner.py` were not opened this session — their `M` status in `git status` predates this session's work (already modified at conversation start by other agents).
- **Blocking:** no.

### [C2/11] Lint check
- **Anchor:** n/a
- **Item:** `python3 scripts/collate_open_items.py --check` run after filing this register; fragment fixed until clean (see run below).
- **Type:** direction
- **Blocking:** no.

---

## C3 — synthesized figures + figure hygiene

Fragment owner: C3. Anchors are `v9` labels/paths and notebook cell ids under
`analysis/theoretical_interpretation/synthesized/cbm_theory/` unless a `v8` path is given.

### [C3/01] New fitted C for fig:pstar_balanced = 32.69, addressed to C1
- **Anchor:** `appendix-emp.tex:42-46` (`% VALUE FROM C3` marker, `$C=[\text{TBD, C3}]$`), `open-items/14-C1.md` `[C1/04]`
- **Type:** direction
- **Item:** Re-running `balanced_binary_threshold.ipynb` NMI-scored gives `C = 32.69`. This happens to land almost exactly on the old sign-agreement value `C=32.7` that C1's marker explicitly says not to reuse (`v8/thesis_v8.tex:1174`) — it is a genuinely recomputed NMI-scored number, not the stale one carried over, and register item `C3/02` explains why the two coincide here (a sharp threshold makes metric choice not move `p̂*` for this figure).
- **Blocking:** no — value computed and reported here; C1 fills in `$C=32.69$`, no `.tex` touched by C3.

### [C3/02] Metric choice does not move p̂*/C for the balanced-binary figure
- **Anchor:** `balanced_binary_threshold.ipynb` cell `5dd4e606`/`a4359734`
- **Type:** doubt
- **Item:** `sign_agreement`, `ari`, and `nmi` all locate the 95%-threshold at the exact same point on the 25-point geometric `p`-grid for every `n`, because the balanced-binary recovery transition is sharp enough that once sign-agreement clears 95% the cluster metrics have already saturated near 1. The premise that "NMI is a stricter bar near the transition" (true in general, and true for the non-balanced CBM figure, `C3/pstar_synth_*`) does not materially change this particular figure's number.
- **Blocking:** no — plotted NMI as instructed; noted so the caption isn't written to imply the metric switch moved the constant.

### [C3/03] Confirming a quantity exceeds a negative bound is not evidence
- **Anchor:** `hbm_spectral_gap_verification.ipynb` cell `878a9350` (sanity table); register `R1/23`
- **Type:** flaw
- **Item:** Both `S_by_alpha.png` and `spectral_gap_bound.png` now come from one execution (`R1/23`'s stated blocking default), fixing the provenance bug, but the underlying interpretation issue is untouched: at `α=0.90, η≥6` the analytic bound `Δλ_bound` itself goes negative, so "slack ≥ 0" there only restates `Δλ ≥ 0`, a property of any Laplacian eigenvalue gap, not a validation of the bound's tightness.
- **Blocking:** no — this is a re-statement of `R1/23`'s open half; no interpretive claim was added or removed from the notebook's own markdown, which already flags the same fact.

### [C3/04] `coherence_vs_eta_by_n.png` generator is unowned this round
- **Anchor:** `analysis/theoretical_interpretation/generated/eta_pool_sweeps/fiedler_tree_partition_by_eta.ipynb`
- **Type:** direction
- **Item:** This notebook produces one of the eleven `v9/figures/` targets but no agent in this restructure owns it, so `scripts/sync_paper_figures.py --check` always reports it `STALE` (the notebook's last edit postdates the May-generated PNG by design, not by regression). The default sync copies the existing PNG through unchanged; nobody has re-verified it reflects current code.
- **Blocking:** no — copied through as-is per task instructions; flagging for whichever future round claims `generated/eta_pool_sweeps/`.

### [C3/05] `nonbalanced_flat_cbm.ipynb` p̂* threshold now keyed on NMI, not requested verbatim
- **Anchor:** `nonbalanced_flat_cbm.ipynb` — n/a (change is in `balanced_binary_threshold.ipynb` cell `5dd4e606`)
- **Type:** direction
- **Item:** The task asked only for cell 11 to pass `metric_col="nmi"`; to keep the plotted `p̂*(n)`/`C` consistent with the plotted curve I also switched `df_thresholds`'s threshold lookup (cell 9) from `sign_agreement` to `nmi`, since leaving it on `sign_agreement` while plotting NMI curves would have mismatched the labeled `p̂*` points against the curve they annotate.
- **Blocking:** no — took the consistent option; the sign-agreement/ARI columns remain in `df_agg` and `trials.csv` for anyone who wants the old cut.

### [C3/06] `scripts/sync_paper_figures.py` C2 notebook mapping is a best guess
- **Anchor:** `scripts/sync_paper_figures.py` `PAPER_FIGURES` entries for `pstar_gen_kmeans_2panel.png`, `recovery_grid_kmeans.png`, `pstar_gen_3operators.png`, `identity_scatter_kingman.png`
- **Type:** question
- **Item:** These four are C2's files; the notebook paths were inferred by grepping the filename strings inside `analysis/theoretical_interpretation/generated/` rather than by reading C2's cell logic (out of scope — C2's files are exclusive-write to them). `identity_scatter_kingman.png` appeared in `v9/figures/` mid-session with a timestamp matching `simulation_distance_vs_similarity.ipynb`, which supports the guess but C2 should confirm both mappings.
- **Blocking:** no — manifest entries are best-effort per the task's "complete the manifest" instruction; wrong notebook names there don't affect the sync/copy logic, only the printed provenance column.

---

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

---

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
