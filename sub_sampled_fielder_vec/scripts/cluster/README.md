# Running on a Linux box over ssh

**Handing this to someone else? Send them [QUICKSTART.md](QUICKSTART.md)** — the same
workflow, copy-paste only. This file is the reference behind it.

Onboarding is three commands. No scheduler is assumed — this is a plain ssh machine; under
SLURM, put the same `run_real_sweep.py` line in an `sbatch` script.

```bash
git clone -b sub_sampled_STDR https://github.com/itaygonnen/spectral-tree-inference.git
cd spectral-tree-inference
bash sub_sampled_fielder_vec/scripts/cluster/bootstrap.sh                      # 1. environment
bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh "6000 taxa" --shared  # 2. data
source .venv/bin/activate && cd sub_sampled_fielder_vec
python scripts/run_real_sweep.py --list                                        # 3. check
```

## 1. `bootstrap.sh` — environment

Creates `.venv`, installs `requirements-cluster.txt`, then `pip install -e . --no-deps`,
then verifies every import and lists the cohorts it found. Run it from the repo root.

`--no-deps` matters: `setup.py` lists `oct2py` (needs a system Octave), `toytree`,
`seaborn` and `sphinx`, none of which are on the import path of these experiments —
but `python-igraph` **is** (`spectraltree/utils.py` imports it at module level), so it is
in the requirements file.

Needs Python ≥ 3.10 (pydantic v2). If the default is older: `module avail python`, then
`PYTHON=python3.11 bash .../bootstrap.sh`.

## 2. `get_data.sh` — data, from Google Drive to the machine that computes

Pulls a cohort with `rclone` directly onto the cluster: resumable, no laptop in the
middle, and the one-time authorisation works from Windows, macOS or Linux. `--help` prints
the `rclone config` walkthrough, including the headless case (`rclone authorize "drive"`
on any machine with a browser, paste the token back). `--shared` is for a folder shared
with your account rather than sitting in My Drive.

Data lands in `<repo>/data/cohorts/<name>/`, and a cohort is exactly:

```
data/cohorts/6000 taxa/fasta/random_tree_0001.fasta   ...
data/cohorts/6000 taxa/newick/random_tree_0001.nwk    ...
```

matched by filename stem. Lay the Drive folder out that way once and every later fetch is
correct. `$STR_DATA_DIR` overrides the location if the data belongs on scratch.

`extract_archive.sh <tarball> --name "6000 taxa"` is the third route, and the usual one
when someone hands over the raw archive: it pulls every true tree plus a chosen slice of
the alignments (default 0001-0100, `--pattern` for more) out of `sim_trees_*.tar.gz`
straight into the cohort layout, in one pass.

`push.sh user@host --data "6000 taxa"` is the alternative when the data is only on a
laptop: it rsyncs the repo and the cohort over ssh, excluding caches, results and `.git`.

## 3. Check, then run

```bash
python -m analysis.utils.real_selftest   # ~1 min, writes to a scratch dir
```

Confirms the data is readable and the whole pipeline runs end to end before you start
anything long.

## 3b. Run

Interactive — the same launcher as on a laptop, option `d`:

```bash
tmux new -s str            # so a dropped ssh does not kill the run
python scripts/interactive_run.py
```

Long jobs — the non-interactive twin, which asks nothing and survives logout:

```bash
nohup python scripts/run_real_sweep.py \
    --cohort "6000 taxa" --stage screen --workers 16 > logs/real_screen.log 2>&1 &

nohup python scripts/run_real_sweep.py \
    --cohort "6000 taxa" --stage sweep --cohort-rule valid_S > logs/real_sweep.log 2>&1 &
tail -f logs/real_sweep.log
```

Both stages cache per tree and resume: kill either and re-run the same command. The menu
and the script share those caches.

### Threads

The screen runs one process per tree, each pinned to a single BLAS thread
(`OMP_NUM_THREADS=1`, set for you) — otherwise every worker opens a full thread pool and a
64-core box runs *slower* than serial. Size `--workers` by cores and memory: each worker
holds the alignment plus three m×m matrices, ~1.2 GB at m=6000.

The sweep is a single process and keeps every thread it can get, so do not export
`OMP_NUM_THREADS=1` in a shell you launch a sweep from.

## 4. Results back

```bash
rsync -avzP user@host:'~/spectral-tree-inference/sub_sampled_fielder_vec/analysis/notebooks_cache/' \
  sub_sampled_fielder_vec/analysis/notebooks_cache/
```

The notebooks then re-plot from the caches with no recompute.
