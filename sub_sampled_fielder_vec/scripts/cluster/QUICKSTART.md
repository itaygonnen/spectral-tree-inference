# Running the real-data experiments — quickstart

Linux box over ssh, no scheduler needed. Everything below is copy-paste.

## 1. Get the code (once)

```bash
git clone -b sub_sampled_STDR https://github.com/itaygonnen/spectral-tree-inference.git
cd spectral-tree-inference
bash sub_sampled_fielder_vec/scripts/cluster/bootstrap.sh
```

Builds `.venv`, installs the dependencies, checks every import. Needs Python ≥ 3.10
(`module avail python`, then `PYTHON=python3.11 bash …` if the default is older).

## 2. Get the data (once per dataset)

### 2a. Copy the dataset file into the repo's `data/` folder

The dataset is one compressed file, `sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz`
(24 GB), in the shared Google Drive folder — **<FILL IN: Drive folder link>**.

Put it in `spectral-tree-inference/data/` — the same folder the unpacked data will live
in:

```bash
# from a computer that already has the file
scp sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz \
    user@cluster:~/spectral-tree-inference/data/

# or download it on the cluster from a direct link
cd ~/spectral-tree-inference/data
wget -O sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz "<FILL IN: url>"
```

It should be ~24 GB and pass this check — a transfer that stopped early looks fine until
you try to open it:

```bash
cd ~/spectral-tree-inference
ls -lh data/*.tar.gz
gzip -t data/*.tar.gz && echo "file is intact"
```

### 2b. Unpack it

From the repo, with no path needed — it takes the file sitting in `data/`:

```bash
cd ~/spectral-tree-inference
bash sub_sampled_fielder_vec/scripts/cluster/extract_archive.sh --name "6000 taxa"
```

That writes **all 3000 true trees** and, by default, **alignments 0001-0100** into
`data/cohorts/6000 taxa/` (~3 GB; unpacking everything would be ~90 GB). The compressed
file can stay where it is, or be deleted once you have what you need.
For more alignments: `--pattern 'random_tree_0[0-4]*.fasta'` gives 0001-0499. It reads the
whole compressed file once, so it takes a few minutes and it is worth choosing how many
alignments you want before running it.

**Shortcut if `rclone` is set up on the cluster**: this downloads from Drive and lays out
the folders in one step, no manual copying:

```bash
bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh "6000 taxa" --shared
```

`--help` walks through the one-time `rclone config`, including the headless case
(authorise on any machine with a browser, paste the token back).

Either way the result is `data/cohorts/<name>/{fasta,newick}/`, which is all the code
looks for. This is also the readiness check — if the cohort is listed with the right tree
count, everything is in place:

```bash
source .venv/bin/activate && cd sub_sampled_fielder_vec
python scripts/run_real_sweep.py --list
#   '6000 taxa': 100 trees, m=6000, L=5000  -> .../data/cohorts/6000 taxa
```

## 3. Run

Two stages: **screening** first (per tree: which split each operator reads off the full
matrix, and is it a real edge of the true tree), then the **recovery sweep** (re-read
that split from sub-sampled matrices, NMI vs `p`).

**Use the launcher.** It asks a handful of questions with sensible defaults, shows what
is already done, prints the estimated cost, and runs both stages:

```bash
tmux new -s str                     # keeps running if the connection drops
source ~/spectral-tree-inference/.venv/bin/activate
cd ~/spectral-tree-inference/sub_sampled_fielder_vec
python scripts/interactive_run.py
```

Pick `real data` → the cohorts (`1,2` runs both) → `screening`, and afterwards the same
launcher again for `recovery sweep`. Detach from tmux with `ctrl-b d`, come back with
`tmux attach -t str`.

The same thing without any questions, if you would rather queue it and log out — every
option above has a flag, and both share the same cache:

```bash
nohup python scripts/run_real_sweep.py --cohort "1000 taxa,6000 taxa" \
    --stage screen --workers 16 > logs/screen.log 2>&1 &

nohup python scripts/run_real_sweep.py --cohort "1000 taxa,6000 taxa" \
    --stage sweep --cohort-rule valid_S --max-eta 20 > logs/sweep.log 2>&1 &

tail -f logs/sweep.log
```

Both stages cache per tree and are **safe to interrupt**: re-run the same command and it
continues. Sizing: a screening worker needs ~1.2 GB at m=6000; the sweep is one process
and wants all the threads it can get, so do not export `OMP_NUM_THREADS=1` for it.

Cost, measured: screening ≈ 75 s per tree at m=6000; the sweep ≈ 30 min per tree on a
20-point × 10-rep grid (the `L` arm dominates). `--p-points` and `--reps` scale it
linearly, and the launcher prints an estimate before it starts.

## 4. Send the results back

Everything from one run is in a single directory:

```
results/real_data/runs/<timestamp>-<name>/
    curves.csv          median, quartiles and a bootstrap CI per cohort and p  <- plot from this
    per_tree.csv        the same numbers per tree, unaggregated
    screening.csv       every tree: eta and validity under each operator
    recovery_curve.png  median NMI vs p, every cohort, both operators
    config.json  summary.json  run.log  experiment.log
```

A few hundred KB, and tracked by git:

```bash
git add results/real_data && git commit -m "cluster run" && git push
```

(Or just `scp -r` that one directory.) `summary.json` records `status: completed |
interrupted | failed`, so a run cut short is not mistaken for a finished one.
