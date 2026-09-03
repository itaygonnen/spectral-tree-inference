# Running on a Linux box over ssh

Three commands the first time, two every session after that. No scheduler assumed — this
is for a plain ssh machine; under SLURM, wrap the same `run_real_sweep.py` line in an
`sbatch` script.

## 1. Push (from the laptop)

```bash
cd sub_sampled_fielder_vec
./scripts/cluster/push.sh user@host --data "6000 taxa"
```

`push.sh` rsyncs the repo and, optionally, one real cohort. It excludes `.git`, `.venv`,
`cache/`, `results/`, `data/` and `docs/overleafs/` — the remote copy exists to *produce*
results, not to hold the laptop's. rsync resumes, so a dropped transfer is re-run, not
restarted. The `6000 taxa` cohort is 3.4 GB (100 alignments + 3000 trees); `--data all`
sends every cohort.

To work on more than 100 alignments, copy the source archive up instead and extract it
there — `tar -xzf sim_trees_*.tar.gz '*/MSAs/random_tree_0[1-5]*.fasta'` selects a slice
without unpacking 90 GB.

## 2. Bootstrap (on the cluster, once)

```bash
ssh user@host
cd ~/spectral-tree-inference
bash sub_sampled_fielder_vec/scripts/cluster/bootstrap.sh
```

Creates `.venv`, installs `requirements-cluster.txt`, then `pip install -e . --no-deps`.
The `--no-deps` matters: `setup.py` lists `oct2py` (needs a system Octave), `toytree`,
`seaborn` and `sphinx`, none of which are on the import path of the experiments — but
`python-igraph` **is** (`spectraltree/utils.py` imports it at module level), so it is in
the requirements file. The script verifies every import and lists the cohorts it found.

Needs Python ≥ 3.10 (pydantic v2). If the login node's default is older:
`module avail python`, then `PYTHON=python3.11 bash .../bootstrap.sh`.

## 3. Run

```bash
source .venv/bin/activate && cd sub_sampled_fielder_vec
python scripts/run_real_sweep.py --list
```

**Interactive** — same launcher as on the laptop, option `d`:

```bash
tmux new -s str            # so a dropped ssh does not kill the run
python scripts/interactive_run.py
```

**Long jobs** — the non-interactive twin, which asks nothing and survives logout:

```bash
nohup python scripts/run_real_sweep.py \
    --cohort "6000 taxa" --stage screen --workers 16 > logs/real_screen.log 2>&1 &

nohup python scripts/run_real_sweep.py \
    --cohort "6000 taxa" --stage sweep --cohort-rule valid_S > logs/real_sweep.log 2>&1 &
tail -f logs/real_sweep.log
```

Both stages are cached per tree and resumable: kill either at any point and re-run the
same command to continue. The interactive menu and the script share those caches.

### Threads

The screen runs one process per tree and pins each to a single BLAS thread
(`OMP_NUM_THREADS=1`, set for you) — otherwise every worker opens a full thread pool and
a 64-core box runs *slower* than serial. Pick `--workers` from cores and memory: each
worker holds the alignment plus three m×m matrices, ~1.2 GB at m=6000.

The sweep is a single process, so it keeps all the threads it can get; do not export
`OMP_NUM_THREADS=1` in a shell you launch a sweep from.

## 4. Pull results back

```bash
rsync -avzP user@host:'~/spectral-tree-inference/sub_sampled_fielder_vec/analysis/notebooks_cache/' \
  sub_sampled_fielder_vec/analysis/notebooks_cache/
```

The notebooks then re-plot from the caches with no recompute.
