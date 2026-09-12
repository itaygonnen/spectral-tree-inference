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
    plt.plot(g.p, g.nmi_L_median, label=f"{dataset} L(S)")
    plt.plot(g.p, g.nmi_B_median, "--", label=f"{dataset} B=HDH")
plt.xscale("log"); plt.legend()
```

`per_tree.csv` holds the same numbers before aggregation, for per-tree spread.

## Columns

`screening.csv` — one row per tree. `eta` is the partition imbalance (larger clan /
smaller clan) of that operator's split of the full matrix; `valid` is 1 when that split is
a real single-edge bipartition of the true tree. `_L` is the Fiedler vector of `L(S)` cut
by k-means, `_B` the leading-|λ| eigenvector of `B = HDH` cut by sign.

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

**Why two references.** `vs_fullmatrix` asks whether sub-sampling kept what the complete
matrix saw; `vs_truetree` asks whether the complete matrix was seeing the right thing.
They can disagree sharply -- 0.68 against the full matrix beside 0.013 against the tree on
the m=1000 benchmark -- and the second is the harder question.

## Producing a run

```bash
cd sub_sampled_fielder_vec
python scripts/interactive_run.py          # real data -> datasets -> screening, then sweep
# or unattended, both sizes in one run:
nohup python scripts/run_real_sweep.py --dataset "1000 taxa,6000 taxa" \
    --stage sweep --dataset-rule valid_S --prefix overnight &
```

Screening has to run before a sweep can be gated on it. Both are cached per tree and
resumable: interrupt either and re-run the same command. See
`sub_sampled_fielder_vec/scripts/cluster/README.md` for the ssh workflow.
