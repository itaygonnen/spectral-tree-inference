### [ABS/01] The abstract reads as machine-written, though its claims are sound
- **Anchor:** sections/0-abstract.tex:5
- **Source:** manual
- **Category:** ai-trace
- **Type:** phrasing
- **Severity:** material
- **Context:** The whole abstract: the cost problem, the sampling idea, the main rate, the
  role of imbalance, and the simulations.
- **Item:** It reads as machine-written: fundamentally bottlenecked, offers a computational
  workaround, In this work we frame, Our primary contribution is, By utilizing, Furthermore,
  and three rendering-X participles in one sentence. The claims themselves are sound, and
  the patch keeps every one.
- **Affects:** none
- **Patch:** review/patches/ABS-01.json
- **Status:** applied
- **Coverage:** the whole abstract body is replaced; the checker scored it 12.5 percent, the
  highest density in the document, yet flagged only one of its sentences
- **Verified:** yes -- every claim in the replacement checked against the body: the rate,
  the fixed-imbalance condition, the square-root-m saving, the two mechanisms and the cubic
  dependence all survive unchanged
- **Applied:** 2026-11-23
- **Decision:** approved -- approved in chat
