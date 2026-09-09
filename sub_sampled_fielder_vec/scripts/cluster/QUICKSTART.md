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

### 2a. Put the archive on the cluster

The dataset ships as one file, `sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz`
(24 GB), which lives in the shared Google Drive folder — **<FILL IN: Drive folder link>**.
Put it anywhere with room; `~/` or a scratch filesystem are both fine, and nothing below
depends on where.

```bash
# from a laptop that has the file
scp sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz user@cluster:~/

# or download it on the cluster, if rclone is configured (see get_data.sh --help)
rclone copy gdrive:"<FILL IN: path in Drive>" ~/ -P

# or, from a direct download link
wget -O sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz "<FILL IN: url>"
```

Check the file arrived whole before unpacking — a truncated transfer fails silently until
tar hits the end:

```bash
ls -lh sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz   # expect ~24 GB
gzip -t sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz && echo "archive intact"
```

### 2b. Unpack it into a cohort

Unpack it into the layout the code reads:

```bash
bash sub_sampled_fielder_vec/scripts/cluster/extract_archive.sh \
    sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz --name "6000 taxa"
```

That pulls out **all 3000 true trees** and, by default, **alignments 0001-0100** (~3 GB;
the full set is ~90 GB). More: `--pattern 'random_tree_0[0-4]*.fasta'` for 0001-0499.
One pass over the archive, a few minutes — `tar` must decompress the whole stream to find
the members, so choose the pattern before running it, not after.

**If the data is in Google Drive instead**, pull it straight onto the cluster with
`rclone`:

```bash
bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh "6000 taxa" --shared
```

`--help` walks through the one-time `rclone config`, including the headless case
(authorise on any machine with a browser, paste the token back).

Either way the result is `data/cohorts/<name>/{fasta,newick}/`, which is all the code
looks for. Confirm with:

```bash
python sub_sampled_fielder_vec/scripts/run_real_sweep.py --list
```

## 3. Check it works

From the repo root (the directory holding `setup.py`):

```bash
cd ~/spectral-tree-inference        # wherever you cloned it
source .venv/bin/activate
cd sub_sampled_fielder_vec
python -m analysis.utils.real_selftest      # ~1 min, writes only to a temp dir
```

Expect four `ok` lines and `all 4 checks passed`. It screens and sweeps two trees for
real, so it proves the data, the environment and the output writers all work together
before you start anything long.

## 4. Run

Two stages: **screening** first (per tree: which split each operator reads off the full
matrix, and is it a real edge of the true tree), then the **recovery sweep** (re-read
that split from sub-sampled matrices, NMI vs `p`).

Interactive — answer a few questions, defaults are sensible:

```bash
tmux new -s str            # so a dropped connection does not kill the run
python scripts/interactive_run.py
```

Unattended, survives logout:

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

## 5. Send the results back

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
