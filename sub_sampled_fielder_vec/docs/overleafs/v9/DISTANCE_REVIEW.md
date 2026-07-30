# Referee report — the distance route (§2 construction, §4 Theorem 2, Appendix F)

Files reviewed: `sections/problem.tex` (the $\mathcal B$ construction paragraph),
`sections/main-result-dist.tex`, `sections/proof-outline.tex` (closing paragraph),
`sections/appendix-G.tex` (App. F, complete). Pre-review snapshot:
`scratchpad/distance-pre-review/`.

---

## 1. Summary of general observations

The distance route is mathematically the more interesting half of the paper: a second operator
on the same leaves, reaching the same bipartition through a different spectral mechanism, with a
split-decomposition identity (`lem:split-decomp`) that makes $\mathcal B\preceq 0$ transparent
and is proved cleanly here. That lemma and the balanced-case spectrum
(`prop:G-spectrum`) are the section's strongest material and are correct as written.

The section is not, however, in submittable condition, for four reasons of descending severity.

**(i) The displayed rate is mis-normalized, and three downstream claims contradict each other
because of it.** `eq:main-rate-dist` carries the factor $(1+\eta)^3$ — which arises only *after*
the gap bound $\Delta\lambda\ge m_{\min}\cdot\text{margin}$ has been substituted — while
simultaneously dividing by the *unsubstituted* $\Delta\lambda^{(\mathcal B)2}$ and by $m$. The
same synthesis that produces `eq:main_rate` on the Laplacian side produces
$p>C(1+\eta)d_0^2\,m\log m/\Delta\lambda^{(\mathcal B)2}$ here, with $m$ in the *numerator*. The two differ by $8m^2/(1+\eta)^2$ — six orders of magnitude at $m=10^3$.
Every inconsistency in the section traces back to this: the proof in `app:dist-main` dutifully
substitutes the mis-normalized expression, obtains $\Theta(\log m/m^3)$, and disowns the
$\Theta(\log m/m)$ claim, while `rem:gap-scope` in §4 and the closing paragraph of App. F both
assert that the two routes *share* $\Theta(\log m/m)$. With the normalization repaired, the
substitution $\Delta\lambda^{(\mathcal B)}=m(d_0-d_1)/2$ yields $p\,m/\log m=\text{const}$ —
the shared exponent is true, and the section can stop apologizing for it.

**(ii) Theorem 2 rests on a step the text itself declines to prove.** `app:dist-pernode` states
that the within-clan cancellation "is ported as an assumption the donor makes, not a proven fact
for $\mathcal B$", and correctly explains why the Laplacian argument does not transfer:
`lem:C-inc` is a row-sum identity and $\mathcal B$ has no row-sum-zero structure. The proof then
begins "Granting it". A theorem may be conditional — this one already is, on
$\Delta\lambda^{(\mathcal B)}>0$ — but the condition must be *stated*, numbered, and visible in
the theorem's hypotheses. As written, a reader who stops at §4 believes Theorem 2 is proved.

**(iii) A numbered theorem is an unverified paraphrase.** `thm:dist-partition` is attributed to
Griffing's thesis, carries the non-quantitative hypothesis "well separated from the runner-up",
and is followed by the admission "We were unable to verify the theorem's exact hypotheses
against the primary source … should be read as a paraphrase, not a transcription." Numbered
environments assert; paraphrases do not. The unconditional Rayleigh-quotient bound
$|\lambda_1(\mathcal B)|\ge 2\tau_{e^*}n_1n_2/m$ that follows it *is* proved here and should be
separated out and kept.

**(iv) The manuscript still addresses its own editorial history.** The body text refers to "the
donor", "the donor manuscript", "Unlike the donor's summary", and cites internal register items
("register item [B1/06]", "open-items [B1/03]") in eleven places. These are artifacts of the
port. A referee reads them as evidence that the section was assembled rather than written.

Beyond these, the section is over-long relative to its content: the roadmap paragraph enumerates
seven subsections in a single sentence chain, and several paragraphs argue with the source rather
than stating results. Roughly a fifth of App. F can go without losing mathematics.

---

## 2. Formal and logical issues

**F1 — `eq:main-rate-dist` mixes two normalizations.** `main-result-dist.tex:35–37`.
The shared synthesis inequality is $8\sigma^2\log m<\text{gap}^2\cdot(\min_i|u_i|)^2$. With
$\min_i|u_i|^2=1/(m\eta)$ and $\sigma^2=\tilde\rho^{(\mathcal B)}(\eta)/p$ this gives
$p>8\eta\,\tilde\rho^{(\mathcal B)}(\eta)\,m\log m/\Delta\lambda^{(\mathcal B)2}=8(1+\eta)d_0^2m\log m/\Delta\lambda^{(\mathcal B)2}$ (the $\eta$ from the margin $1/(m\eta)$ cancels the $1/\eta$ in the load). On the Laplacian
side the extra step $\Delta\lambda\ge m_{\min}(\rho-\eta S_{\rm out}^{\max})$, with
$m_{\min}^2=m^2/(1+\eta)^2$, converts $m\log m/\Delta\lambda^2$ into
$(1+\eta)^2\log m/[m\cdot\text{margin}^2]$ — this is precisely where `eq:main_rate`'s
$(1+\eta)^3$ and its $m$-in-the-denominator come from (I verified the reproduction numerically:
ratio $1.000$ across $m,\eta$). Writing the distance rate with $(1+\eta)^3$ *and*
$\Delta\lambda^{(\mathcal B)2}$ *and* $m$ in the denominator applies that conversion twice.
**Correction:** state the theorem in the abstract-gap normalization and put the substituted form
in a corollary — which is also the natural home for the $\Theta(\log m/m)$ claim and the
constant comparison.

**F2 — the proof contradicts §4 on the exponent.** `appendix-G.tex:348–354` versus
`main-result-dist.tex:86–88`. Under F1 the appendix's $\Theta(\log m/m^3)$ is an artifact of the
mis-normalized display, not a finding; it must be withdrawn, not reconciled. Note that the
passage as written is also self-undermining: it presents $\Theta(\log m/m^3)$ — a *better* rate
than the Laplacian's — as a reason for caution, when a sufficient condition that is easier to
satisfy would be a strength.

**F3 — the comparison's direction is unsupported by its own display.** `main-result-dist.tex:86–88`
claims the distance route "pays a strictly worse multiplicative constant, of order
$d_0^2/\beta_0^2$". Against the corrected rate the penalty is
$d_0^2(\rho-\eta S_{\rm out}^{\max})^2/[\beta_0^2(d_0-d_1)^2]$, not $d_0^2/\beta_0^2$ alone:
the margins differ between the two routes as well as the noise scales. At $\beta_0=0.05$,
$d_0=0.693$, $d_1=0.2$, margin $=0.6$ this is $\approx 284$, not $192$.

**F4 — the load factor $\tilde\rho^{(\mathcal B)}$ is asserted, and the natural computation
disagrees with it.** `appendix-G.tex:306–321`. The text sets
$\tilde\rho^{(\mathcal B)}(\eta)=(1+\eta)d_0^2/\eta$, i.e. the Laplacian load with
$\beta_0\to d_0$. But the Laplacian load descends from increments
$X_{ij}\propto S_{ij}(v_i-v_j)$, whereas $(E_{\mathcal D}u)_i=\sum_j(E_{\mathcal D})_{ij}u_j$
weights $u_j$ alone. Carrying that through with $c_A^2=\eta/m$, $c_B^2=1/(m\eta)$ gives, for the
binding clan $B$, $\mathrm{Var}=\tfrac{1-p}{p}d_0^2\,\eta/(1+\eta)$ — load
$\eta^2 d_0^2/(1+\eta)$ after synthesis, not $(1+\eta)d_0^2/\eta$. The asymmetry is real: the
Laplacian's differencing makes both clans contribute equally; $\mathcal B$'s does not.
**I have not changed the exponent**, because the underlying cancellation is itself unproven
(F5) and the resolution is the author's call. The theorem now carries $\tilde\rho^{(\mathcal B)}$
symbolically, so both readings are expressible, and the discrepancy is registered.

**F5 — Theorem 2's proof is conditional on an unnumbered borrowed claim.**
`appendix-G.tex:309–321` ("As asserted in the donor manuscript, not re-derived here"; "This claim
is ported as an assumption"), used at `appendix-G.tex:322` ("Granting it"). **Correction:**
promote to a numbered assumption and list it in Theorem 2's hypotheses alongside
$\Delta\lambda^{(\mathcal B)}>0$.

**F6 — `lem:C-first` is mis-described where it is applied.** `appendix-G.tex:341–344` says the
Neumann lemmas "bound $\|\hat u-u\|_\infty$ by $\|E_{\mathcal B}\|_2/\Delta\lambda^{(\mathcal B)}$
up to a coherence factor $(1+\eta)^{3/2}$". `lem:C-first` (`appendix-C.tex:40–45`) bounds
$|w_i^{(1)}|\le(|(E u)_i|+\|E\|_2|u_i|)/(\Delta\lambda-\|E\|_2)$ — a per-node quantity plus an
operator-norm term, with no coherence factor. The distance proof should mirror the Laplacian
synthesis it claims to follow: reduce to the per-node condition
$|(E_{\mathcal B}u)_i|<\Delta\lambda^{(\mathcal B)}|u_i|$ (noting, as `sec:outline` does, that
the second term is discarded) and then substitute the per-node bound.

**F7 — the coherence identity is claimed "operator-agnostic".** `appendix-G.tex:342–344`.
`prop:coherence` is proved from `lem:A-decomp`, the *Laplacian's* Fiedler coordinates under the
CBM, and $\mu(U)$ is defined for $U=[v^{(1)},v^{(2)}]$. The transfer to
$[\mathbf 1/\sqrt m,u^{(1)}]$ is legitimate *provided* $u^{(1)}$ is piecewise constant with the
same clan sizes — which `prop:G-spectrum` establishes only in the balanced case, and
`rem:G-unbalanced` explicitly leaves open at general $\eta$. As written the proof assumes at
general $\eta$ what the section elsewhere declares unknown.

**F8 — $c$ is declared and never bound.** `main-result-dist.tex:38–39`: "$C,c>0$ are universal
constants for which the failure probability is $O(m^{-c})$", while the statement's only
probabilistic qualifier is the macro `\hp`. **Correction:** "with probability at least
$1-O(m^{-c})$" inside the statement.

**F9 — entry bound omits the reweighting.** `appendix-G.tex:276–278`: "$E_{\mathcal D}$ has
entries bounded by $r(\mathcal T)$". Under inverse-probability weighting the entries are bounded
by $d_0/p$; `appendix-G.tex:308` already uses $d_0/p$. The remark is missing the factor.

**F10 — `thm:dist-partition`: numbered assertion, unverified content, unquantified hypothesis.**
`appendix-G.tex:220–238`. See §1(iii).

**F11 — a forward dependency.** `appendix-G.tex:82` uses `lem:split-decomp` ("By
\Cref{lem:split-decomp} below") to justify the sign convention two subsections before the lemma
is stated. The lemma is short and self-contained; state it first.

**F12 — duplicate and dead labels.** `main-result-dist.tex:80` carries two labels on one remark
(`rem:gap-scope`, `rem:dist-tradeoff`), both referenced, both resolving to the same number.
`appendix-G.tex:357` attaches `\label{sec:dist-compare}` to an unnumbered `\paragraph`, where a
`\Cref` would print whatever counter last advanced.

**F13 — an unreconciled empirical claim.** `appendix-G.tex:250–252` reports a median imbalance of
$34.1$ for Fiedler-on-$S$. Since `cor:infeasible` puts the certificate's reach at
$\eta=O((m/\log m)^{1/3})$, a reader will ask whether the comparison is being run in a regime
where neither theorem certifies anything. One clause suffices; silence invites the objection.

---

## 3. The distance part's introduction

The introduction is App. F's opening paragraph (`appendix-G.tex:36–45`) plus, in the main text,
the $\mathcal B$ construction paragraph in §2 and the one-sentence bridge that ends
`main-result.tex`.

*What works.* Placing the construction in §2, so that §4 can put the two theorems back to back,
is the right structure: it establishes both operators before either rate is claimed. The bridge
sentence ("the hypotheses of the two operators are not nested … a *sibling* … rather than a
corollary") does real work in two lines.

*What does not.* The App. F opening is a seven-item inventory in a single sentence chain
("…the construction in full: how $\mathcal B$ is built …, the split-decomposition identity …,
the population-level spectral landmarks …, a quoted partition-recovery result …, the sampling and
noise model …, the per-node concentration step …, and the proof …"). It tells the reader the
order of the furniture, not why the room is worth entering. An appendix opening should state the
one structural fact that makes the route work, then the one thing that is missing.

*Suggested replacement* (theme, not final wording):

> This appendix develops the distance route in full. Its structural basis is a single identity:
> for any tree-additive $\mathcal D$, the centred operator $\mathcal B=H\mathcal D H$ is a
> nonnegatively-weighted sum of rank-one negative-semidefinite terms, one per edge
> (`lem:split-decomp`). Hence $\mathcal B\preceq0$, its leading eigenvector sits at the most
> negative eigenvalue, and each edge contributes in proportion to $\tau_e n_1(e)n_2(e)$ — so the
> operator is, by construction, tuned to deep and balanced splits. Sections … build $\mathcal D$
> and $\mathcal B$ from data, compute the population spectrum exactly in the balanced flat CBM,
> and carry the sub-sampling analysis through to Theorem 2. Two population inputs remain
> unresolved at general imbalance: a lower bound on $\Delta\lambda^{(\mathcal B)}$, and the
> per-node cancellation of within-clan increments; both are isolated as explicit hypotheses
> rather than absorbed silently.

That last sentence does what the current opening lacks — it tells the referee, before they find
it themselves, exactly where the section's weight is not yet carried.

Two further additions worth making in §2 rather than the appendix: one sentence on *why* a second
operator is worth the reader's time (the distance route connects to the classical
numerical-taxonomy literature, and its sign rule targets a different edge — currently the reader
learns this only in App. F), and one clause fixing $\mathcal D$'s relation to $r(\mathcal T)$,
since the diameter reappears in the noise analysis.

---

## 4. The distance part's conclusion

The conclusion is the closing "Honest comparison with the Laplacian route" paragraph
(`appendix-G.tex:357–377`), together with `rem:gap-scope` in §4.

*What works.* The two-item structure — coherence barrier unchanged, noise scale
$\beta_0\to d_0$ — is the right frame, and the refusal to double-count $d_0^2/\beta_0^2$ against
a separate $r(\mathcal T)^2$ penalty is a genuine correction to the source that deserves to
survive.

*What does not.* First, the heading. "Honest comparison" concedes that something elsewhere was
not; a referee notices. Call it "Comparison with the Laplacian route". Second, the paragraph
argues with the source ("Unlike the donor's summary…") instead of stating a result. Third, and
most importantly, it is a comparison of constants that never states the comparison's conclusion
in one line. Fourth, it inherits F2's contradiction.

*Suggested closing* (theme):

> Where the two routes are directly comparable they differ in exactly one quantity. The coherence
> barrier is identical — $\mu(U)=(1+\eta)/2$ depends only on the piecewise-constant shape of the
> structural vector, so both routes lose feasibility at $\eta=\Omega((m/\log m)^{1/3})$. The rate
> exponent is identical: substituting the balanced-case gap into Theorem 2 gives
> $p=\Theta(\log m/m)$, as in Theorem 1. What differs is the constant, and it differs in one
> place: the cross-clan noise scale $\beta_0$ is replaced by $d_0=-\log\beta_0$, so the distance
> route pays $d_0^2(\rho-\eta S_{\rm out}^{\max})^2/[\beta_0^2(d_0-d_1)^2]$ for the same target
> — a factor of roughly $300$ at $\beta_0=0.05$. This is a statement about noise tolerance for a
> *fixed* split, and it is consistent with the separate observation that $\mathcal B$'s sign rule
> often selects a *better-balanced* split to begin with: the two comparisons answer different
> questions, and only the first is settled by the rates.

That formulation states the conclusion (identical exponent, one identified constant penalty),
quantifies it, and keeps the tension with the edge-selection finding — without arguing with a
manuscript the reader has not seen.

**Observation, outside this scope.** The paper has no Conclusion or Discussion section at all. At
Annals/JAMS level a short closing section is conventional and would carry the three things
currently scattered in remarks: what is proved versus certified (the sufficient/necessary
asymmetry), the two open population inputs (a general-$\eta$ gap bound for either operator, and
the $\ell_\infty$ criterion's slack), and the empirical gap between the $(1+\eta)^3$ inflation and
the measured factor of $\approx4.6$ across $\eta\in[1,15]$. Recorded here, not drafted.

---

## 5. Punctuation, with reasoning

Line references are to the pre-review snapshot.

1. **`appendix-G.tex:49–51`** — "Fix the same binary phylogenetic tree $\mathcal T$ on $m$ leaves
   and GTR rate matrix $Q$ as in \Cref{asm:gtr}, and the population similarity matrix $S$ of
   \Cref{sec:problem}." Two coordinated objects joined by "and", then a comma plus a third
   "and": the comma reads as closing a list that has not been opened, and "as in `asm:gtr`" is
   left straddling. → "Fix the binary phylogenetic tree $\mathcal T$ on $m$ leaves, the GTR rate
   matrix $Q$ of \Cref{asm:gtr}, and the population similarity matrix $S$ of \Cref{sec:problem}."
   (Serial list, one comma per item, modifier attached to its own noun.)

2. **`appendix-G.tex:58–65`** — a 90-word sentence with three parenthetical qualifications and a
   trailing "See open-items register entry [B1/06] for the rationale." Split after "already used
   elsewhere in this paper"; delete the register pointer (§4 of this report).

3. **`appendix-G.tex:82–84`** — "By \Cref{lem:split-decomp} below, $\Bop\preceq0$ for any
   tree-additive $\Dmat$, so its eigenvector of largest *magnitude* sits at the most negative
   eigenvalue $\lambda_1(\Bop)$; write $\uone$ for a unit-norm choice of it." An eigenvector does
   not "sit at" an eigenvalue. → "…so the eigenvalue of largest magnitude, $\lambda_1(\Bop)$, is
   the most negative one; let $\uone$ denote a corresponding unit eigenvector." The semicolon is
   correct (two independent clauses, second elaborating).

4. **`appendix-G.tex:103–105`** — "The following is what makes $\Bop\preceq0$ rigorous, and is
   kept as a short lemma because it is used twice: once to justify…, and once inside…" The comma
   before "and is kept" separates a compound predicate sharing one subject; delete it. The comma
   before "and once inside" is likewise unnecessary in a two-item list.

5. **`appendix-G.tex:116–117`** — "since $(s_e)_i^2=1$ for every $i$, $\delta_e$ has zero
   diagonal automatically, and $(\delta_e)_{ij}=1$ iff…" Three clauses, two commas, one "and":
   the middle comma makes "automatically" seem to modify the third clause. → "…for every $i$,
   $\delta_e$ automatically has zero diagonal; moreover $(\delta_e)_{ij}=1$ iff…"

6. **`appendix-G.tex:268–269`** — "Idempotence of $H$ and $\opnorm H=1$ give immediately".
   Adverb stranded after the verb, splitting it from its object. → "immediately give".

7. **`appendix-G.tex:274–281`** — `rem:dist-entryscale`: "On the distance side, $E_{\Dmat}$ has
   entries bounded by $r(\mathcal T)$, which by \Cref{app:dist-landmarks} equals $d_0$ exactly
   under \Cref{asm:clock} --- not an independent $O(\log m)$ contamination term on top of $d_0$."
   The dash-clause negates something the sentence never asserted, so it reads as a rebuttal to an
   absent interlocutor (it is one — the source manuscript). Delete the clause; state the bound
   ($d_0/p$, per F9) and its consequence.

8. **`appendix-G.tex:309–312`** — "*As asserted in the donor manuscript, not re-derived here*:
   within-clan cancellation is claimed to apply "as in the Laplacian proof" because…" Italicized
   editorial aside, colon splice, and quoted fragment from an unnamed source. Replaced wholesale
   by the numbered assumption (F5).

9. **`appendix-G.tex:329–335`** — "The scalar … is identical across all $i$; it shifts every
   coordinate of $\huone$ uniformly and so does not affect sign recovery, only the effective
   threshold in \eqref{eq:dist-entrywise}." The "only…" fragment after a comma dangles off
   "affect". → "…and so does not affect sign recovery; it changes only the effective threshold
   in \eqref{eq:dist-entrywise}."

10. **`main-result-dist.tex:46–51`** — "The two statements differ in how far their parameters can
    be resolved: \eqref{eq:main_rate} is closed-form in the model, whereas no explicit lower
    bound on $\DlamB$ valid at general imbalance $\eta$ is available (\Cref{rem:gap-scope}), so
    \eqref{eq:main-rate-dist} is left in terms of the abstract gap." Three coordinators
    (`whereas`, `so`) on one spine, with the parenthetical between the second clause and its
    consequence. Split at "available.".

11. **`main-result-dist.tex:89–94`** — "This is a statement about noise tolerance for recovering a
    *fixed* edge, and it does not contradict the separate empirical finding … that $\Bop$'s sign
    rule can select a *better-balanced* edge than $L$'s Fiedler vector before any sub-sampling
    occurs: one claim compares tolerance for a fixed target, the other which target is chosen,
    and both belong in the paper." 60 words, a colon, and an elliptical "the other which target"
    that needs a verb. → two sentences; supply "the other, which target is chosen". "Both belong
    in the paper" is a note to a co-author, not a claim; delete.

12. **`appendix-G.tex:368–372`** — "Unlike the donor's summary, this is stated as a *single*
    factor $d_0^2/\bze^2$, not a product with a separate $r(\mathcal T)^2$ term: under
    \Cref{asm:clock}, $r(\mathcal T)=d_0$ exactly (\Cref{app:dist-landmarks}), so counting … and
    … as independent penalties double-counts the same quantity." Keep the mathematics, drop the
    comparison to the source and the colon-chain: "Under \Cref{asm:clock}, $r(\mathcal T)=d_0$
    exactly, so the diameter penalty and the noise-scale penalty are the same quantity and must
    not be multiplied."

Two systematic notes. **`---` before a subordinate clause is used ten times in App. F** where a
comma or semicolon is correct; an em dash should mark an interruption, not a routine
qualification. **"Where" is used causally** ("Where the two routes are directly comparable, they
differ…") — acceptable once, mannered when repeated; prefer "When" or a plain conditional.

---

## 6. Round 2 — independent referee, and corrections to §§1–5 of this report

An independent referee (fresh context, given the pre-edit snapshot, the edited files and this
report, and asked to audit it) raised twelve findings. Ten are accepted and applied; two are
matters of presentation where I applied a variant. It also identified four errors **in this
report**, all of which are real. Recorded honestly, because two of them had propagated into the
manuscript.

### 6.1 Errors in this report (and one it caused in the manuscript)

- **F3's arithmetic mixed two values of $\beta_0$.** It read "at $\beta_0=0.05$, $d_0=0.693$" —
  but $d_0=-\log\beta_0$, so $\beta_0=0.05$ gives $d_0=2.996$; $0.693$ belongs to $\beta_0=0.5$.
  The manuscript then printed "roughly $300$ at $\bze=0.05$, $d_0-d_1\approx0.5$", and at those
  parameters $S_{\rm in}=e^{-(d_0-0.5)}=0.030$, so $\rho-\eta\Sout^{\max}=-0.070<0$:
  `asm:margin` is violated and `eq:main_rate` is unavailable there. **Fixed:** the price is now
  quoted at the paper's own synthesized-regime values ($\beta_0=0.05$, $S_{\rm in}=0.9$, hence
  $d_0=2.996$, $d_1=0.105$), where the factor is $275.0$ and the margin is $+0.80$. Also, "not
  $192$" was a slip — the pre-review text said $1.92$ and $530$, never $192$.
- **F4's explanation of the asymmetry was wrong.** It claimed "the Laplacian's differencing makes
  both clans contribute equally". `prop:C-var` gives $\sigma_A^2=\eta\,\sigma_B^2$ — they do not.
  The genuine asymmetry is that the ratio runs the *other way* for $\Bop$, which strengthens
  rather than weakens the "clan $B$ binds" step.
- **F4's disposition was wrong, not merely cautious.** Registering the load discrepancy while
  shipping the disputed value as `asm:dist-cancel` converted the wrong side of an open question
  into a hypothesis of the main theorem. See 6.2(1): the assumption was false in the one case
  where it can be checked, which made Theorem 2 and its corollary vacuous. It has been deleted.
- **F12 was reported as fixed but net-regressed**, and "$8m^2/(1+\eta)^2$" should be
  $m^2/(1+\eta)^2$ (both displays carry the same unnamed $C$; the $8$ was my synthesis constant).

### 6.2 Accepted findings, applied

1. **`asm:dist-cancel` was false, and unnecessary.** The direct computation of the cross-clan
   load gives $\eta d_0^2/(1+\eta)$, not $(1+\eta)d_0^2/\eta$ — a factor $4$ at $\eta=1$, in the
   balanced flat CBM where `prop:G-spectrum` *proves* $\uone=(e_A-e_B)/\sqrt m$, so the
   assumption was refutable, not merely unproven. More importantly no cancellation is needed:
   since $\norm{\uone}_2=1$ and $\Dmat_{ij}\le d_0$, the total per-node variance including
   within-clan terms is at most $d_0^2/p$ in one line. This is now `lem:dist-var`, the assumption
   is gone, and the resulting rate is *stronger* by a constant:
   $p>C d_0^2\eta\,m\log m/(\DlamB)^2$.
2. **Piecewise constancy of $\uone$ was used but not stated** (the margin $1/\sqrt{m\eta}$ and the
   sign rule both need it; `asm:cbm` gives it only for the Laplacian's Fiedler vector). Now
   hypothesis (ii) of Theorem 2.
3. **No operator-norm bound on $E_{\Dmat}$ existed anywhere**, though `lem:C-series`,
   `lem:C-first` and `lem:C-remainder` all require $\opnorm{E}<\Dlam$. Added `lem:dist-opnorm`
   (matrix Bernstein, mirroring `lem:B-bernstein`): $\opnorm{E_{\Bop}}\le\opnorm{E_{\Dmat}}
   \lesssim d_0\sqrt{m\log m/p}$, and the strong-signal condition is now discharged explicitly in
   Step 1 of the proof.
4. **`lem:C-first` bounds $w^{(1)}$, not $w$**, and `lem:C-remainder` was never invoked. The proof
   now cites both and discards the lower-order terms explicitly, as `sec:outline` does.
5. **"A uniform shift does not affect sign recovery" was false**, and the de-meaning remedy was
   vacuous: $\hat\Bop\1=0$ already forces $\huone$ to have mean zero, and a shift exceeding $c_B$
   flips all of clan $B$. The centering scalar is now *bounded* — $\opnorm{E_{\Dmat}}/\sqrt m
   \lesssim d_0\sqrt{\log m/p}$, the same order as the per-node term, hence a constant factor
   (`rem:dist-centering`).
6. **The Bernstein linear term used $c_B$ where the summands weight $\uone_j$** ($c_A$), an error
   of a factor $\eta$; corrected, and the crossover it implies is exactly $p\ge\eta\log m/m$,
   which is now hypothesis (iii) rather than an appeal to "the variance-dominated regime" —
   itself circular at abstract gap, since $p=\Theta(\log m/m)$ presupposes $\DlamB=\Theta(m)$.
7. **"Both routes lose feasibility at the same $\eta=\Omega((m/\log m)^{1/3})$" was a non
   sequitur**: that threshold is read off `eq:main_rate` *after* the gap substitution, which is
   unavailable for $\DlamB$; and `cor:dist-balanced` fixes $\eta=1$, so it cannot exhibit
   $\eta$-growth. The coherence *identity* transfers where hypothesis (ii) holds; the threshold
   does not, and the text now says so.
8. **The real-data figures ($34.1$ vs $1.51$) cited a cohort that does not appear in this paper**,
   and the hedge attached to them was false ($1.51$ is a balanced split). Both numbers are
   removed; `rem:dist-compare-empirical` now makes the same point structurally, from
   `lem:dist-edge`'s $\tau_e n_1(e)n_2(e)$ weighting.
9. **The flat CBM is not realizable on a binary tree** — constant within-clan distance forces each
   clan to be a star. `lem:split-decomp` is now stated for an arbitrary tree (its proof never used
   binarity), and `rem:flat-realization` records the two-star realization.
10. **Griffing's result is now a prose citation**, not a numbered proposition with a
    non-quantitative hypothesis; `lem:dist-edge` carries what this paper actually proves.

### 6.3 Applied as a variant

- The referee asked for the "two normalizations must not be mixed" warning to leave the proof
  environment: it is now `rem:dist-normalization`, after the proof. Also fixed there: $C'$ was
  undefined, and "the failure probability is the union bound" was a category error.
- Dead labels: `app:comp-split` renamed `app:dist-split` (it was in the `app:comp*` family of
  Appendices A–C) and added to the roadmap, which had also listed the subsections in the wrong
  order. Equation labels on displayed equations are kept even where currently unreferenced —
  standard practice, and they are the natural targets for future cross-references.

### 6.4 Referred to the author, not fixed here

**The identity $r(\mathcal T)=d_0$ makes `cor:tolerance` vacuous.** `cor:tolerance` (in
`bridging.tex`, proved in `appendix-A.tex`) lower-bounds $\Sin^{\min}\ge e^{-r(\mathcal T)}$ and
concludes $\eta\le\mathcal O(m^c)$. But under `asm:clock` the diameter *is* the cross-clan
distance, $r(\mathcal T)=d_0=-\log\Sout^{\max}$, so $e^{-r(\mathcal T)}=\Sout^{\max}$ and the
derived tolerance collapses to $\eta<1$. The distance appendix's identification is correct; what
it exposes is that `cor:tolerance` needs the maximum *within-clan* distance, not the diameter.
This sits outside the distance part and touches a headline claim of §3 ("under balanced
generative models the gap survives polynomial imbalance"), so it is registered
(`open-items/13-B1.md` [B1/24]) rather than edited.

---

## 7. Round 3 — second independent referee

A second referee (fresh context, given the snapshot, the edited files and §6, and asked to audit
both) confirmed this round's ten repairs by independent derivation — `lem:dist-var`'s bound,
the rate's $\eta^1m^1(\DlamB)^{-2}$ powers, `cor:dist-balanced`, `lem:dist-opnorm`'s variance
proxy, the crossover, and the $275$ (arithmetic exact, parameters `asm:margin`-feasible) — and
then found a defect that **both** previous passes missed because both audited the distance route
*against* the Laplacian route instead of auditing the step they share.

### 7.1 The remaining substantive defect: the discarded Neumann remainder (accepted)

`eq:dist-entrywise` compares an $\ell_\infty$ error against a margin of order $(m\eta)^{-1/2}$,
but the only bound available for the second-order term is through the operator norm,
$\norm{w^{(2)}}_\infty=O(\opnorm{E}^2/\Delta\lambda^2)$ (`lem:C-remainder`). At the paper's own
parameters ($\eta=1$, $\bze=0.05$, $\Sin=0.9$) and its own rate, `lem:B-bernstein` and `lem:gap`
give $\opnorm{\EL}/\Dlam\approx7.07$ — **constant in $m$, not vanishing** — so the discarded
remainder is $\approx50$ while the margin $\to0$. I verified this numerically at
$m=10^3,10^6,10^9$: the ratio is identical at every size, so it is not a small-$m$ artifact.
Requiring $\norm{w^{(2)}}_\infty=o(\min_i|\uone_i|)$ forces
$p\gtrsim d_0^2m^{3/2}\eta^{1/2}\log m/(\DlamB)^2$, i.e. $\Theta(\log m/\sqrt m)$ in the balanced
case — a factor $\sqrt m$ above the claimed rate.

This is the paper's own `R1/07`, which records the same objection on the Laplacian side. What
round 2 got wrong was to cite `lem:C-remainder` in the distance proof as though invoking it
discharged the issue: that converted an omission into an explicit invalid step. Applied:

- Theorem 2 now carries **hypothesis (iii)**, remainder control
  $\opnorm{E_{\Bop}}\le\varepsilon_m\DlamB$ with $\varepsilon_m\le c_0(m\eta)^{-1/4}$, and Step 2
  shows that (iii) is exactly what makes *both* discarded corrections lower order against the
  margin. The old hypothesis (iii) ($p\ge\eta\log m/m$) is deleted, being implied by the rate via
  $\DlamB\le md_0$ — which also discharges `lem:dist-opnorm`'s own unverified hypothesis.
- `rem:dist-remainder` states what (iii) costs ($\Theta(\log m/\sqrt m)$ if discharged through
  the operator norm), names the fix (an entry-wise leave-one-out bound,
  \cite{abbe2020entrywise}), and records that **`thm:main-sim` has the same gap** — its stated
  rate is not self-contained either.

I did not restate Theorem 1. Its exponent is the paper's central claim and the repair is a new
lemma, not an edit; that is the author's call, and it is now registered as [B1/27].

### 7.2 Also accepted and applied

- **Theorem 2 had no instance inside the paper's model class.** Hypotheses (i)–(ii) are verified
  only in the flat CBM, which `rem:flat-realization` shows is two stars — not the binary tree §2
  fixes. `cor:dist-balanced` now says so, and states that (i)–(ii) are open at general imbalance
  and on binary trees.
- **`cor:infeasible` was false as stated.** Counterexample: $m=10^6$, $\eta=100$, $\bze=0.001$,
  $\Sin=0.9$ satisfies every hypothesis and gives RHS $=1.8\times10^{-4}$. The threshold needs
  $\bze$ and the margin held away from $0$ independently of $\eta$; the invariant form is
  $\eta^3\bze^2/(\rho-\eta\Sout^{\max})^2=\Omega(m/\log m)$. Both are now stated.
- **Tree-additivity of $\Dmat$ needs `asm:clock`**, which was in force but uncited: $S_{ij}$
  depends on the MRCA depth, and $2\tau_{ij}$ is the path length only under ultrametricity.
- **$c_A,c_B$ were under-determined**: one equation, two unknowns. Orthogonality
  $\uone\perp\1$ (from $\Bop\1=0$) supplies $a c_A=b c_B$ — and it is also what gives the two
  clans opposite signs, which the sign rule needs.
- **The per-node bound was an unnumbered display cited as a result**: promoted to
  `lem:dist-pernode` with its hypotheses.
- **`rem:dist-centering` argued both sides of a non-issue** (the scalar is a term inside
  $(E_{\Bop}\uone)_i$, not a shift of $\huone$, and $\hat\Bop\1=0$ makes the shift impossible);
  the framing is gone, the bound is kept as `eq:dist-EBu`, and Step 2 no longer claims an
  identity where only the bounds agree.
- **The $275$ was presented as the ratio of the two theorems' right-hand sides.** It is the ratio
  of the two *parameter groups* at $\eta=1$; the stated constants ($8(1+\eta)^3$ versus $4C$) are
  not comparable, and by 7.1 neither is pinned down. Now said explicitly.
- **The coherence clause in the comparison was dead weight** — no quantity in
  `eq:main-rate-dist` depends on $\coh$ — and is demoted to a parenthesis.
- **The appeal to `cor:tolerance`** in `rem:infeasible-mech` now carries the caveat that the
  distance appendix's $r(\mathcal T)=d_0$ weakens that corollary (§6.4, [B1/24]).

### 7.3 Errors in §6 of this report, per the second referee — all four accepted

$S_{\rm in}=e^{-2.496}=0.0824$, not $0.030$ (the conclusion, that `asm:margin` fails at those
parameters, is unaffected); §6.2(4) claimed a fix that was not one (7.1); §6.2(6) presented a
redundant hypothesis as necessary; "six orders of magnitude" is 5.4. §6.2(9) stated the two-star
fact without drawing its consequence for Theorem 2's instantiability (now 7.2).

### 7.4 Verdict carried forward

The second referee's verdict — *not yet sound enough for submission* — stands, and I agree with
it. After this round the distance route's remaining defects are all of one kind: it is honest
about what it assumes, but two of its three hypotheses are verified only in a configuration
outside the paper's tree class, and the third (remainder control) is unverified on both routes.
The single blocking item is 7.1, and it is not specific to the distance part.
