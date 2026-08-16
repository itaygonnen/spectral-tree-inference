# The Distance Route — standalone

`distance_route.tex` is a **self-contained** document. It has its own preamble, its own
bibliography (`refs.bib`, in this directory), and shares no `\input`, no label namespace and no
`.bib` file with `../thesis_v9.tex`. Every label is prefixed `d:` and every macro is
`\providecommand`'d, so the file could later be `\input` into a larger document without a clash.

## Build

```bash
cd sub_sampled_fielder_vec/docs/overleafs/v9/distance-route
latexmk -pdf -interaction=nonstopmode distance_route.tex
grep -c "undefined" distance_route.log     # expected: 0
```

Ten pages. A truncated `.aux` from an interrupted run causes
`File ended while scanning use of \@newl@bel`; fix with `latexmk -C`, then rebuild.

## What it proves

Let `D` be the matrix of path distances between the `m` leaves of a weighted tree,
`H = I − 11ᵀ/m`, and `B = HDH`. The sign rule takes the eigenvector of `B` at the eigenvalue of
largest *magnitude* (`B ⪯ 0`, so that is the most negative one) and partitions the leaves by sign.

Two exact structure theorems carry the argument.

- **Theorem 1 (split covariance).** `B = −Σ_e w_e ê_e ê_eᵀ` with `ê_e = H·1_{A_e}` normalised and
  `w_e = 2 ℓ_e a_e b_e / m`. Proved in two lines from the path decomposition `D = Σ_e ℓ_e S_e` and a
  rank-one centring lemma. Gives the spectral gap, the trace identity `Σ_e w_e = (1/m)·1ᵀD1`, and an
  explicit preference for balanced bipartitions through the factor `a_e b_e`.
- **Theorem 2 (Gower matrix as a Laplacian pseudoinverse).** `−½HDH = (L/L_RR)†`, the Schur
  complement of the tree Laplacian (conductances `1/ℓ_e`) onto the leaves. Combined with
  Stone & Griffing 2009 Thm 3.3 this makes the sign rule **unconditionally valid**: its output is
  always the bipartition induced by a single edge of the tree, with no quantitative hypothesis.

Against these, a sampling model in which each pair is observed independently with probability `p`
and estimated by inverse-probability weighting: an operator-norm bound, a Davis–Kahan step, and a
one-step refinement whose per-leaf statistic concentrates by a scalar Bernstein inequality.
Conclusion: `p = Θ(log m / m)` suffices to reproduce exactly the partition the fully observed
operator returns.

## How it differs from the deleted `sections/appendix-G.tex`

That appendix (compiled as App F; deleted from the working tree, `\input` commented out at
`thesis_v9.tex:151`) factored `B` through the Bandelt–Dress canonical decomposition, stated its only
closed-form gap for balanced splits, and carried a per-node cancellation hypothesis that its own
referee report records as false and removed. Nothing here is ported from it.

Four substantive differences:

| | deleted App F | this document |
|---|---|---|
| structure | Bandelt–Dress canonical decomposition | Theorem 1, proved from the path decomposition; Theorem 2, proved from the tree Laplacian |
| Griffing | cited as a partition theorem | the pseudoinverse identity is **proved**; the partition theorem is cited to its actual source, Stone & Griffing 2009 Thm 3.3 |
| claim | the split, gap-conditional and balanced-only | *a* split, unconditionally (Cor. 6); the quantitative work is about surviving sub-sampling, not about identification |
| ℓ₂ → ℓ∞ | resolvent expansion, remainder discarded | Davis–Kahan + one-step refinement; no series to truncate |

The word *primary* does not appear. The document does not claim the recovered split is the root
split, and does not claim it agrees with the split a similarity Laplacian returns — §7 derives why
the two operators select different edges.

## Verification

The identities asserted in §2–§3, the landmark table in §4 and the `p ∝ κ² log m / m` scaling of
Theorem 3 were checked numerically on random weighted trees to a maximum absolute error of
`1.4e-14`. The check script is not committed; it reconstructs random binary trees with independent
uniform edge lengths and verifies, in order: the path decomposition, the rank-one centring lemma,
the split-covariance identity, the trace identity, the split Gram matrix, the Gower-pseudoinverse
identity (including the eigenvector correspondence `|cos| = 1` and `λ_top(G) = 1/λ₂(L*)`), the
validity corollary over 200 trees, and the two-stage recovery rate over `m = 256 … 2048`.

## Bibliography

`refs.bib` is deliberately independent of `../references.bib`, whose `Griffing2012` entry carries an
explicit "could not verify" comment. Every entry here has a checked DOI or resolved URL; the Griffing
thesis entry has the working NCSU repository URL. Three sources whose bibliographic data could only
be obtained second-hand (Zaretskii 1965, Buneman 1971, Critchley & Fichet 1997) are listed as
comments at the foot of the file and are not cited.
