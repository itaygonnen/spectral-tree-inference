# results/real_data/

```
runs/<timestamp>-<name>/     ONE directory per run -- every dataset that run covered
    config.json              what was asked for: datasets, gate, p-grid, reps, commit, time
    summary.json             headline numbers per dataset, and how the run ended:
                             status = completed | interrupted | failed, plus
                             trees_requested vs trees_done per dataset
    run.log                  everything the run printed
    experiment.log           timestamped, tagged detail: per-tree reference splits, eta,
                             k-means warnings, per-tree timing and ETA
    screening.csv            every tree of every dataset: eta and validity per operator
    per_tree.csv             every tree x every p: all metrics, both arms
    curves.csv               per dataset x p: median, std, quartiles and a bootstrap
                             interval for the median, across trees  <- plot from this
    recovery_curve.png       median NMI vs p, every dataset, both arms, both references

_cache/<dataset>/             machine state, not a deliverable: screen.npz and one .npz per
                             swept tree, so a killed run resumes and a second run reuses
                             what the first computed
```

Download or mail a whole run by taking its directory: a few hundred KB, no other file is
needed. `$STR_RESULTS_DIR` moves the parent of both roots (a scratch filesystem, say).

## Interrupting a run

Safe at any point. Trees are cached one at a time, so every tree that finished is kept and
the run still exports the CSVs and the plot for those trees. The tree in flight is
discarded -- nothing partial is ever written. Re-run the same command and it continues
from the cache. `summary.json` records `status: "interrupted"`, which is the only thing
distinguishing that run's CSVs from a complete run's shorter ones.

## Plotting from a run

`curves.csv` is already aggregated: one row per (dataset, p), one column per
`<metric>_<arm>_<stat>` with `stat` in `median, std, q25, q75, lo95, hi95, n`.

`q25`/`q75` bound the middle half of the **trees** — the spread of the dataset, which does
not shrink as trees are added. `lo95`/`hi95` are a bootstrap interval for the **median**
itself — how firmly this dataset pins the curve down, and that does shrink. `std` is kept
for continuity but is the weakest of the three here: NMI is bounded and often bimodal (a
tree either recovers its split or does not), so a mean ± std band leaves [0, 1].

```python
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv("curves.csv")
for dataset, g in df.groupby("dataset"):
    plt.plot(g.p, g.nmi_L_vs_fullmatrix_median, label=f"{dataset} L(S)")
    plt.plot(g.p, g.nmi_B_vs_fullmatrix_median, "--", label=f"{dataset} B=HDH")
plt.xscale("log"); plt.legend()
```

`per_tree.csv` holds the same numbers before aggregation, for per-tree spread.

## Columns

`screening.csv` — one row per tree, one column block per operator. `eta` is the partition
imbalance (larger clan / smaller clan) of that operator's split of the full matrix;
`valid` is 1 when that split is a real single-edge bipartition of the true tree.

| arm | operator | cut |
|---|---|---|
| `L` | Fiedler of `L(S) = Deg(S) − S` | k-means **or** sign |
| `Lsym` | Fiedler of `L_sym = I − Dg^-1/2 S Dg^-1/2` | k-means **or** sign |
| `B` | leading-\|λ\| eigenvector of `B = HDH` | sign (there the sign pattern *is* the partition) |

"k-means or sign" is decided per tree, not per run: whichever cuts that tree's reference
more evenly wins, and the choice is recorded in `rule_<arm>` with both candidates'
`eta_<arm>_kmeans` / `eta_<arm>_sign` beside it. k-means alone routinely isolates a single
taxon on real data (1/999) — a real pendant edge, so a validity gate passes it, but a
split no sub-sample can recover.

### Column names

Every per-p column reads `<metric>_<operator>_vs_<reference>`:

- operator: `L` = Fiedler of `L(S)`, `B` = leading eigenvector of `B = HDH`
- reference: `vs_fullmatrix` = that operator's own split of the complete matrix;
  `vs_truetree` = the true tree's top bipartition (what used to be called `gt`)
- metric: `nmi`, `ari`, `agreement` (% of taxa on the same side), `dot` (with the
  reference eigenvector), `signagreement` (L only -- on B the sign pattern IS the
  partition, so its `agreement` already is that number)

Plus, per p and operator: `eta_<op>` (imbalance of the split recovered at that p) and
`split_small_<op>` (size of its smaller side). Per tree: `rule_L` (which cut rule the L
arm used), `eta_ref_L`, `eta_ref_B`.

NMI is what the paper figure plots; the rest are stored so a follow-up question needs no
re-run.

### Linear-algebra diagnostics

The quantities the original pipeline-A `results.json` carried beside the partition
metrics, recorded at every p on both arms. "full" is the complete matrix an arm reads
(`S` on the L arm, `D` on the B arm), "sub" is its bootstrap-averaged sub-sample, and
"op" is the operator built from it (`L(S) = Deg(S) - S`, or `B = HDH`):

| column | what it is | original name |
|---|---|---|
| `sigma2_full_<op>` | σ₂ of the cross-clan block of the full matrix | `sigma2_avg_M` |
| `sigma2_sub_<op>` | the same on the averaged sub-sample | `sigma2_avg_S` |
| `opnorm_err_{mean,median,std}_<op>` | ‖sub − full‖₂ over the replicates | `*_operator_norm_error` |
| `rank_full_<op>`, `rank_sub_{mean,median,std}_<op>` | numerical rank ‖A‖²_F / ‖A‖²₂ | `empirical_rank_M`, `_S` |
| `rank_op_full_<op>`, `rank_op_sub_{mean,median,std}_<op>` | the same for the operator | `empirical_rank_L_M`, `_L_S` |

Every one is O(m²) per replicate — Lanczos for ‖·‖₂, power iteration for the error norm,
a randomised rank-2 SVD for σ₂ — never a dense factorisation. Measured cost on the L arm
at m=1000: +21% (0.48 s → 0.58 s for 4 p × 5 reps); the share falls as 1/m against the
O(m³) eigensolve, so about 4% at m=6000. On by default; `--no-extra-metrics` skips them.

A column is **additive**: it is not part of the cache key, so a tree swept before these
existed still counts as swept and its cells are simply blank. `curves.csv` carries an `n`
per metric, so a mixed export says how many trees actually contributed to each column.

**Why two references.** `vs_fullmatrix` asks whether sub-sampling kept what the complete
matrix saw; `vs_truetree` asks whether the complete matrix was seeing the right thing.
They can disagree sharply -- 0.68 against the full matrix beside 0.013 against the tree on
the m=1000 benchmark -- and the second is the harder question.

## Producing a run

```bash
cd sub_sampled_fielder_vec
python scripts/interactive_run.py          # real data -> datasets -> screening, then sweep
# or unattended, both sizes in one run:
nohup python scripts/run_sweep.py --dataset "1000 taxa,6000 taxa" \
    --stage sweep --dataset-rule both --prefix overnight &
```

Screening has to run before a sweep can be gated on it. Both are cached per tree and
resumable: interrupt either and re-run the same command. See
`sub_sampled_fielder_vec/scripts/cluster/README.md` for the ssh workflow.

The menu and the command line are two ways of filling in the same `RunSpec`
(`src/runners/experiment_run.py`) and both then call one `execute()`, so they produce
identical run directories from identical answers. `--dataset-rule` takes `both`,
`valid_L`, `valid_Lsym`, `valid_B`, `any` or `all`; `valid_S` is still accepted as the
old name for `valid_L`. `both` still means L(S) and B — the pair the figures compare —
so a selection recorded before `L_sym` was screened still means what it said.

The runner is not specific to real data: it takes a `Source` (a name, a loader, a list
of tree ids, a taxon count), and `RealLoader` / `GeneratedLoader` are two ways to
produce one. That is why `screening.csv` and `curves.csv` look the same whether the
trees came from FASTA files or a simulator.
