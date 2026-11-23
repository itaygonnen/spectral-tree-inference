### [STY/02] Six quantities written as if they act, replaced with plain verbs
- **Anchor:** sections/6-topological-scaling.tex:134
- **Source:** manual
- **Category:** register
- **Type:** phrasing
- **Severity:** editorial
- **Context:** Scattered prose across the scaling section, the proof outline and the
  empirical section, describing what mathematical quantities do.
- **Item:** The gap survives imbalance, the Davis-Kahan route forfeits a factor, scalar
  quantities sit behind the constants, and an eigenpair carries a split; one sentence also
  opens with So. The patch puts a plain verb in each place: remains positive, does not
  achieve, determine, persists under.
- **Affects:** none
- **Patch:** review/patches/STY-02.json
- **Status:** applied
- **Coverage:** 6 sites across 3 files; the other uses of carries describe a graph edge
  carrying a weight, which is standard usage and is left alone
- **Verified:** yes -- each site read in place, and line-initial so and and that are
  mid-sentence wraps were excluded by checking what precedes them
- **Applied:** 2026-11-23
- **Decision:** approved -- approved in chat
