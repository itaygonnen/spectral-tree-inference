# Benchmark launcher guide — `scripts/run_benchmark.py`

The operator-comparison benchmark: **how small can the sub-sampling fraction `p` get
before an operator stops recovering the partition the full matrix would give?** It runs
the pipeline behind the `distance_vs_similarity` notebooks — validity screen → per-operator
cohorts → bootstrap `p`-sweep → recovery curve + scale plot — without opening a notebook.

```bash
cd sub_sampled_fielder_vec
python scripts/run_benchmark.py
```

Enter takes the default at every prompt. Nothing expensive happens until you confirm.

> This is **not** `scripts/interactive_run.py` (pipeline A, one matrix at a time — see
> [INTERACTIVE_GUIDE.md](INTERACTIVE_GUIDE.md)). Different tool, different outputs.

---

## The dialogue

### Step 1 — data source

| Choice | What it means |
|---|---|
| `real` | FASTA alignments on disk (+ Newick ground truth if present) |
| `generated` | trees simulated on demand — no download |
| `synthesized` | prints a pointer to the CBM-theory notebooks; this tool does not apply |

`real` scans `../data` and `data/` and prints what it found; pick a number or type a path
to a fasta directory. Without a Newick directory the validity screen is skipped entirely
and the cohort is just the first `n_compare` trees.

### Steps 2–6 — the generated branch

| Step | Prompt | Notes |
|---|---|---|
| 2 | tree model | multi-select: `kingman`, `bd`, `lopsided`, `balanced_binary`. `lopsided`/`balanced_binary` have a fixed topology — only the alignment varies between trees |
| 3 | model properties | `pop_size` (kingman) / `birth_rate`+`death_rate` (bd) / `edge_length` (the deterministic two) |
| 4 | number of taxa | **comma-separated for several sizes** in one run. `balanced_binary` requires powers of two |
| 5 | alignment | sequence length, JC mutation rate |
| 6 | eta pooling | see below |

Every combination of (model × size × η bin) is one **cell**, and each cell gets
`trees per cell` trees. The total is *derived* from those answers and printed —
there is deliberately no "how many of the N trees" prompt on this path, because
truncating a cell-ordered list drops whole models or sizes rather than thinning them.

### Step 6 — eta pooling

A tree's **η** is the imbalance `max(n1,n2)/min(n1,n2)` of the split its reference
partition finds. Answer `n` and you get whatever imbalance the model happens to produce.
Answer `y` and trees are drawn from the shared pool under `cache/pool_sample`, keeping
only those whose η lands within ±1 of a target:

```
Eta targets — comma-separated [1,5,10,15]:
Trees per (model, size, eta) cell [10]:
Max rejection attempts per size, when the pool is short [2000]:
```

It then prints what is already cached, per cell:

```
  pool availability (mu=0.1, L=10000, pop=1.0):
    kingman          n=1000   eta01 10/1  eta05 10/1  eta10 10/1  eta15 10/1   [ready]
```

`ready` cells cost nothing. `build` cells are filled **after** you confirm, by shelling out
to `scripts/build_eta_pool.py` — the canonical pool builder, so anything built here is
visible to the notebooks that read that pool. Rejection sampling for a high η at a large
`n` is slow; that is what the attempt cap is for.

The defaults (μ=0.1, L=10000, pop_size=1.0, kingman) are the parameters **every pool on
disk was built with** — change any of them and the pool key changes, so you are building
from scratch. The availability table warns when a configuration has nothing cached.

Only models `build_eta_pool.py` accepts can be pooled; others skip straight to the
plain per-cell count.

### Step 7 — benchmark parameters

```
  [1] L + k-means
  [2] L_sym + k-means
  [3] B = HDH (Griffing-on-D)
Operators — comma-separated numbers [1,2,3]:
```

An operator you leave out is never screened *or* swept — dropping one saves a full
bootstrap sweep per tree. Then: `p`-values (default 20 log-spaced 0.01→1.0), bootstrap
replicates per `p`, the cohort cap, and the output directory.

---

## What the run does

1. **Screen** — build each tree's full matrix, derive each selected operator's partition,
   and test whether it is a real single-edge split of the ground-truth tree.
2. **Replace invalid draws** — a tree that fails the gate is a wasted draw, not a result.
   Cells short of their requested count get more candidates drawn and screened, in rounds:

   ```
   valid: L=1  G=4  / 4 screened  (per operator — cohorts are independent)
   replacing invalid draws: 3 needed across 3 cell(s) -> screening 3 more
   ...
   valid: L=4  G=9  / 9 screened
   ```

   Pooled cells draw the next unused slots in their bin; plain cells increment the tree
   index, which is unbounded. When a pooled bin runs out it says which cell and by how
   much, rather than quietly shrinking the sample.
3. **Per-operator cohorts** — each operator gets the trees *it* validates. No tree has to
   pass every gate (at high η almost none would), so the curves in one figure generally
   rest on different, differently-sized tree sets — which is why every legend entry and
   every `summary.json` block carries its own `N`.
4. **Sweep** — bootstrap `p`-sweep per tree, each tree swept only for the operators it is
   in the cohort of.
5. **Aggregate** — the two figures and the tables below.

---

## Outputs

Everything lands in the output directory:

| File | Contents |
|---|---|
| `results.json` | **every metric at every `p`**, one row per (cell, operator, `p`) |
| `n<n>_L<L>/results.json` | the same table split per size — the layout `src/utils/merge_results.py` globs for |
| `summary.json` | per operator: tree count, median `p*`, NMI at the smallest `p` |
| `recovery_curve.png` | median NMI vs `p`, ±1 std band, one curve per operator |
| `scale_plot.png` | `p*` vs `r(T)` (tree diameter), per-operator cohorts |
| `compare_sweep.npz` | every per-tree curve: `nmi_<op>`, `rT_<op>`, `trees_<op>` |
| `screen_table.csv` | per tree: per-operator validity + η. Human-readable on purpose |
| `sweeps/<id>.npz` | one file per tree — the resume cache |
| `config.json` | the settings both caches were produced under |

`results.json` columns: `p`, `operator`, `cell`, `model`, `num_taxa`, `eta`, `n_trees`,
`p_star_median`, `p_star_threshold`, and `mean`/`median`/`std` of `nmi`, `ari`,
`agreement`, `sign`, `dot`. Merge a multi-size run into one grid with:

```bash
python -m src.utils.merge_results <run_dir>      # -> results_grid_merged.json
```

The launcher also prints the headline numbers, so you do not have to open a file:

```
  operator                   trees      p*  NMI at smallest p
  B = HDH (Griffing-on-D)        7   0.050              0.995
  L + k-means                    4   0.600              0.025
```

`p*` is the smallest grid `p` whose NMI reaches 0.95, median over the operator's trees.
With a coarse `p` grid it is coarse — read it off a 20-point grid, not a 3-point one.

---

## Resuming, and what refuses to mix

Every tree is cached as it is produced, so **pointing a new run at an existing output
directory resumes it**: screened trees are not re-screened, swept trees are not re-swept.
Ctrl-C is safe.

Guards that raise rather than silently merging incomparable work:

- different `screen` or `sweep` settings (`num_gaps`, `min_split`, `p`-values, reps)
- a different `source.data_key` — for generated runs that covers the model list, sizes, μ,
  alignment length and η bins. Changing μ and resuming used to mix incomparable trees
  under the same ids.

Adding an **operator** to a resumed run is allowed and cheap: the sweep cache is merged
per tree, so a tree swept for `L,G` today and `L_sym` tomorrow keeps all three in one
file. Asking for an operator a cached file lacks is a clean miss, not a stale hit.

Sweep files written before the extra metrics existed still load — they just carry NMI
only, so `results.json` for those trees has no `ari`/`agreement`/`sign`/`dot` columns.
Re-run into a fresh directory if you want them.

---

## Cost

The sweep is `O(n³)` per tree per `p` above `p≈0.1`. Rough shape: `trees × |p| × reps`
dense eigendecompositions of an `n×n` matrix. Screening is one decomposition per tree per
operator on top.

A first run should be small — 2–4 trees per cell, `p` = `0.05,0.2,1.0`, 2 reps — then
raise the numbers against the same output directory, which reuses everything already
computed. At the defaults (20 `p`-values × 10 reps, cohort cap 100) an `n=1000` run is
hours, and stage 1 screens *every* draw before the cap applies.

---

## Where the code lives

| Piece | File |
|---|---|
| the dialogue | `scripts/run_benchmark.py` |
| generated-source prompts, replacement supplier | `src/utils/generated_prompts.py` |
| screen / cohorts / sweep / aggregate | `src/runners/benchmark_pipeline.py` |
| per-tree work units, worker pool | `src/runners/benchmark_workers.py` |
| loaders (FASTA, simulated) | `src/runners/benchmark_loaders.py` |
| eta-pool reads and top-up | `src/runners/eta_pool_bridge.py` |
| tree-id grammar | `src/utils/tree_ids.py` |
| tree + alignment generation | `src/models/generated_trees.py` |
| caches and the resume guards | `src/utils/benchmark_cache.py` |
| figures | `src/utils/benchmark_plots.py` |
| `results.json` | `src/utils/benchmark_results_json.py` |

A tree id is a recipe, not a file: `kingman_n1000_eta05_003` is model `kingman`, 1000
taxa, η bin 5, index 3. Plain ids are regenerated from the index (deterministic); ids with
an η bin are read back from the pool, because a rejection-sampled tree cannot be rebuilt
from a seed.
