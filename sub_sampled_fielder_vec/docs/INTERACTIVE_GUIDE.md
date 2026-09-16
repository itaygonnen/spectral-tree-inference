# Interactive launcher guide

```bash
cd sub_sampled_fielder_vec
python scripts/interactive_run.py
```

One experiment, two data sources. The first question is which:

```
Data:
  [1] real data - FASTA alignments with their true trees (data/tree_sets/)
  [2] generated data - simulated trees and sequences
```

Both branches then ask the same things and run the same code
(`src/runners/experiment_run.py`), so a simulated run and a real one produce the same run
directory, the same CSVs and the same curve — a **median over trees** with an
inter-quartile band. Before 2026-09-16 the simulated branch ran a different experiment
(one tree per size, spread from bootstrap replicates); see *The single-tree flow* below.

## What each branch asks

**Real data** — pick one or more datasets discovered under `data/tree_sets/` (any
`<name>/{fasta,newick}` pair; `$STR_DATA_DIR` overrides the location). Select several with
commas to run every size in turn.

**Generated data** — tree model(s), model properties, taxon counts, sequence length and
mutation rate, and optionally η pooling (keep only trees whose reference split lands in a
given imbalance band).

## Then, identically for both

A **screening status** panel: how many trees each source has, how many are screened, how
many each operator cuts a real tree edge on, the median η, how many clear the η cap, and
how many are already swept.

A **stage**:

- `screening` — split the full matrix with every operator, record η and whether that split
  is a genuine single-edge bipartition of the true tree. Resumable; re-running tops up
  rows that predate an operator without discarding the verdicts they already have.
- `recovery sweep` — for each selected tree, sub-sample at each `p` and re-read the
  bipartition. Needs screening first, because the gate reads its verdicts.

One **configuration**, applied to every source so sizes stay comparable. Press Enter to
accept the defaults, or edit: trees per source, the validity gate, the η cap, the p-grid,
bootstrap reps, operators, the extra diagnostics, a run name, and the output mode.

Finally a **plan** — the p-grid, the operators, and per source the chain
`N screened → N after the gate → N after the η cap`, with an estimated wall time — and a
confirmation. Everything is cached per tree, so the run is safe to interrupt and resume.

## Operators and the gate

| arm | operator | cut |
|---|---|---|
| `L` | Fiedler of `L(S) = Deg(S) − S` | k-means **or** sign, whichever splits that tree's reference more evenly |
| `Lsym` | Fiedler of `L_sym` | same |
| `B` | leading-\|λ\| eigenvector of `B = HDH` | sign (there the sign pattern *is* the partition) |

The gate decides which screened trees enter a sweep: `both` (L and B), `valid_L`,
`valid_Lsym`, `valid_B`, `any`, or `all`. The η cap weighs **the operators the gate
names** — `both` → L and B — so a lopsided `L_sym` split never removes a tree from a run
whose figure compares L and B.

## Output

```
results/real_data/runs/<timestamp>-<name>/
    config.json   what was asked for, plus the commit
    summary.json  headline numbers per source, status, the selection chain
    screening.csv every tree: η and validity, per operator
    per_tree.csv  every tree × every p, all metrics, every arm
    curves.csv    per source × p: median, std, quartiles, bootstrap CI   ← plot from this
    recovery_curve.png
    run.log       everything the run printed
    experiment.log  tagged detail, one line per (tree, arm, p)
results/real_data/_cache/<source>/   screen.npz + one .npz per swept tree (resumable)
```

`results/real_data/README.md` documents every column.

## Unattended runs

For anything longer than a session, use the non-interactive twin — same questions as
command-line flags, same code, same run directory:

```bash
python scripts/run_sweep.py --list
python scripts/run_sweep.py --dataset "6000 taxa" --stage sweep \
    --dataset-rule valid_L --max-eta 20 --dry-run     # see the selection, run nothing
nohup python scripts/run_sweep.py --dataset "6000 taxa" --stage sweep \
    --dataset-rule valid_L --limit 120 > logs/sweep.log 2>&1 &
```

`--dry-run` prints the grid, the operators and the selection chain and exits, which is the
cheap way to choose a gate before committing a cluster job to it. `scripts/run_real_sweep.py`
is kept as a shim for the name in the cluster runbook.

## The single-tree flow

The pre-2026-09 simulated flow — one tree per `(n_taxa, seq_len)`, cached matrices,
re-run-last-configuration, results under
`results/<tree_model>/<sampling>/<timestamp>-<run_name>/n{taxa}_L{seq_len}/` — is off the
menu. Its curve is one tree's, not a population median, so it was never comparable with a
real-data curve.

It is not deleted: `scripts/run_experiment.py` drives the same pipeline without prompts,
and `interactive_run.single_tree_menu()` is the old menu. It remains the only path that
produces the per-p linear-algebra diagnostics (coherence, spectral gap, dk_ratio, IPR) and
the middle-out and guardrail code.

## Troubleshooting

**"no datasets found"** — a dataset is `<root>/<name>/fasta/*.fasta` beside
`<name>/newick/*.nwk`. The launcher prints every root it searched and why each folder was
rejected; `python scripts/run_sweep.py --list` prints the same, plus `NOT ALIGNED` when a
FASTA's records differ in length.

**A gate selects nothing** — the screen carries no verdict for the operator it names
(e.g. `valid_Lsym` against a screen written before `L_sym` was added). Re-run
`--stage screen`; it tops up in place. `run_sweep.py` warns about exactly this.

**Checking it all works** — `python -m analysis.utils.selftest` walks every call site of
both branches in about a minute, in a scratch results root.
