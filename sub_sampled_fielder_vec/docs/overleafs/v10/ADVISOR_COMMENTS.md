# Advisor comments on v9 `sections/2-generative-model.tex` — disposition in v10

25 comments, all anchored in `2-generative-model.tex`. That file no longer exists: its content
was split into `sections/2-background.tex` (the objects),
`sections/3-problem-definition.tex` (the STDR algorithm, its cost, the observation model and the
recovery criterion) and `sections/4-assumptions.tex` (the hypotheses). Anchors below are v9 line
numbers.

| # | v9 anchor | v10 location | Disposition |
|---|---|---|---|
| 1 | l.7 opening ¶ | `2-background.tex:7-14` | Rewritten as a `\Cref`-per-subsection roadmap. "the clock that turns branch lengths into times" is gone; the section now says "the substitution model" and, in §4, "molecular clock". |
| 2 | l.21 STDR cite | `2-background.tex:18` | `\cite{aizenbud2023spectral}` on the tree definition replaced by `\cite{semple2003phylogenetics}`. |
| 3 | l.24 "carries" | `2-background.tex:26` | "Every edge … **is assigned** a branch length". |
| 4 | l.36 MRCA undefined | `2-background.tex:32` `def:mrca` | Defined as its own numbered definition before first use. |
| 5 | l.37 figure ref + τ clash | `2-background.tex:47`; notation throughout | Text now says "\Cref{fig:tree_topology} illustrates both conventions on a tree with $m=7$ terminal nodes". **Notation clash resolved**: branch/edge length is now $t_e$ (was $\tau_{hh'}$ and $\tau_e$); $\tau_h$ stays the node height and $\tau_{ij}$ the MRCA height. `eq:transition_matrix` reads $\exp(Q\,t_{(h,h')})$. |
| 6 | l.42 "descending from r" | `2-background.tex:49` | "The two edges **adjacent to the root**". |
| 7 | l.45 cumbersome sentence | `2-background.tex:56-62` | Replaced with the suggested wording: "We denote by $\|A\|$ and $\|B\|$ the sizes of the two sides … $A=\{x_1,\dots,x_{\|A\|}\}$ and $B=\{x_{\|A\|+1},\dots,x_m\}$." |
| 8 | l.46 $\|A\|,\|B\|$ vs $a,b$ | throughout | **Adopted in full.** $a$ and $b$ are gone from the manuscript. $\|A\|,\|B\|$ is used everywhere, including the dense algebra of §6, App A and App C: $c_A=\sqrt{\|B\|/(\|A\|m)}$, $J_{\|A\|\times\|B\|}$, $\1_{\|A\|}$, $m_{\min}=\min(\|A\|,\|B\|)$. |
| 9 | l.46 why this ordering | `2-background.tex:60-62` | Answered: the ordering is a relabelling that makes each side a contiguous index range, which is what puts the split into the $2\times2$ block form `eq:block_forms` rather than scattering it. |
| 10 | l.52 clan definition | `2-background.tex:257` `def:clan` | Replaced with the standard definition — the subtree separable by removing one edge, terminal **and internal** nodes included. |
| 11 | l.56 "share a neighbour" | `2-background.tex:261` | Spelled out: "the roots of the two subtrees are joined to a common neighbouring node", followed by a sentence saying what "the terminal nodes of a clan" means. |
| 12 | l.177 unclear ¶ | — | **Deleted.** The rooted-vs-clan paragraph is gone. The one clause worth keeping (that $S$ is a function of pairwise distances, so it can express exactly the edge cuts of the unrooted tree) survives as one sentence at the end of §2.3. |
| 13 | l.190 "data arrive" | `2-background.tex:186` | "Sequence data **are generated** by running a Markov substitution process". GTR is introduced as "the standard reversible substitution model for aligned sequence data", which states the scope positively rather than as a disclaimer. |
| 14 | l.197 "a stationary" | `2-background.tex:194` | "**there is** a stationary distribution $\pi>0$ with $\pi Q=0$". |
| 15 | l.201 generator/reversible | `2-background.tex:198-205` | Both explained: the generator's off-diagonal entries are instantaneous substitution rates and its rows sum to zero so $\exp(Qt)$ is stochastic; reversibility is written out as $\pi_a[\exp(Qt)]_{ab}=\pi_b[\exp(Qt)]_{ba}$. |
| 16 | l.233 new distance notation | `4-assumptions.tex:40-46` | **Dropped.** The $d(\cdot,\cdot)$ operator is gone; `eq:ultrametric` is now just $d_{ij}=2\tau_{ij}$, with the equal-legs argument in words. |
| 17 | l.243 opening sentence | `2-background.tex:216` | Now opens "The substitution model provides a measure of similarity between every pair of nodes of the tree, which we now define", with the relation to $\tau$ and $Q$ coming after `def:similarity`. |
| 18 | l.256 "collects" | `2-background.tex:246` | "which **consists of** the pairwise similarities over the terminal nodes". |
| 19 | l.258 define the graph | `2-background.tex:244` `def:simgraph` | The definition now leads with the graph — vertices are the terminal nodes, weight matrix is the similarity matrix — and derives the matrix from it. |
| 20 | l.263 placement | `2-background.tex:225-240` | The closed-form derivation now sits between `def:similarity` and `def:simgraph`, and carries the distance matrix `eq:distance_from_similarity` with it. |
| 21 | l.292 "due to Jaffe et al." | `2-background.tex:277` | Attribution moved into the theorem header `{\cite[Lem.~3.1]{jaffe2021spectral}}`; the prose no longer takes the authors as its subject. |
| 22 | l.292 wrong result / explain STDR | `2-background.tex:277, 296, 311` | **Both.** `lem:rank-one` keeps its structural role but loses the "this is what STDR is built on" claim, which now sits on `thm:stdr-split`. STDR itself is explained at the head of §3 with pseudocode (`alg:stdr`). |
| 23 | l.302 "under the molecular clock" | `4-assumptions.tex:22` `asm:clock` | The clock is now a named assumption in §4, cited by `\Cref{asm:clock}` wherever its consequences are used, rather than an aside in §2. |
| 24 | l.302 "more than rank-1" | `4-assumptions.tex:103-107` | **Dropped.** Replaced by the reason: "every pair separated by the primary split has the root as its most recent common ancestor, so `eq:similarity_clock` evaluates all of those pairs at the same height $\tau_r$." |
| 25 | l.305 why not just the CBM | `4-assumptions.tex:56-76` | **Answered by reordering.** §4 now opens with the constant block model of `\cite{Balakrishnan2011clustering}` as *the* model the paper works in, then states the localisation (diagonal blocks arbitrary) as its phylogenetic instance. `def:cfbm` is titled "Constant block model at a bipartition"; the "CFBM" coinage is gone. |

## Beyond the 25

- **Three-way split.** §2 `sec:background` makes no modelling assumption and holds only objects and
  the algorithm; §3 `sec:problem` states the problem; §4 `sec:assumptions` collects every
  hypothesis. §3 and §4 carry **no** subsections: each is one argument, and a heading inside it
  restarts the prose. This matches the venue — in `aizenbud2023spectral`, §1 and §2 (Problem setup,
  614 words) have no subsections either, and they appear only in the method and results sections;
  Eldridge et al. use none at all. §2 keeps three headings over 1539 words of genuinely parallel
  objects. The reading order is the outline as given: topology and the substitution process,
  similarity and distance and the graph they define, clans and `thm:stdr-split`, then the STDR
  pseudocode and the cost that motivates sampling.
- **The Fiedler-vector shape is out of the modelling section.** §4 says what is assumed and what it
  gives in linear algebra (`eq:block_forms`) and stops. The two-plateau eigenvector is
  `lem:A-decomp`(ii) in §6, where the proof uses it; `eq:two_plateau` no longer exists.
- **Figure 2 moved to §4.2**, where $\eta$ is defined — its three subcaptions are stated in $\eta$,
  which §2 has no access to. All floats are now `[t]`.
- **"Primary split" is gone.** The term implied the recovered split, which nothing guarantees. The
  bipartition induced by the root's two edges is now the **root split** — a tree object, named
  where the geometry is genuinely about the root (`prop:cfbm`). `def:cfbm` and `lem:A-decomp` are
  stated at an arbitrary bipartition $\mathcal{X}=A\cup B$; `thm:stdr-split` still promises only
  that $\Pi(S)$ is *a* clan pair.
- **§4 tightened.** The proof of `prop:cfbm` moved to `A-spectral-geometry.tex`; the plateau
  discussion deleted; $\eta$ and $\rho$ merged into one `def:imbalance`. §4 is a preamble carrying
  A1/A2 plus two subsections. C1 keeps its own counter, and a sentence now says why: A1--A2 are
  hypotheses on how the data are generated and already fix every quantity in $\rho$; C1 is an
  inequality between those quantities, checkable on a given tree, that delimits which trees the
  guarantee covers.
- **`clade` removed from §2.** The definitional aside asserting that the root split is a clade did
  no work. `\cite{wilkinson2007clades}` survives on `def:clan`. The one-clause gloss in §1 ("the
  unrooted analogue of a clade") is kept as reader orientation --- say the word and it goes too.
- **Distance route cut.** `B = HDH` is gone from the manuscript. This closes the 6 dangling labels
  v9 carried (`app:dist`, `thm:main-dist`, `cor:dist-balanced`, `lem:dist-edge`,
  `eq:main-rate-dist-balanced`, `rem:gap-scope`): v9 built with 14 undefined references, v10 builds
  with 0. `\Dmat` survives only as the distance matrix of §2.3.
