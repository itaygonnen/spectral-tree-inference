# Running the size series on a cluster

Five taxon counts — 1000, 2000, 4000, 6000, 8000 — each simulated until all four eta
buckets hold 10 trees, then screened, then swept.

Everything is configured by **`cluster/series.json`** — sizes, trees per bucket,
alignment length, the p-grid, workers. Edit that file; the scripts read it.

```bash
git clone <this repo> ~/spectral-tree-inference
cd ~/spectral-tree-inference

python -m venv ~/venvs/str && source ~/venvs/str/bin/activate
pip install -r requirements-cluster.txt      # includes msasim
pip install -e . --no-deps

# IQ-TREE is a binary, not a pip package. No compiler needed.
curl -LO https://github.com/iqtree/iqtree3/releases/download/v3.1.4/iqtree-3.1.4-Linux-intel.tar.gz
tar xzf iqtree-3.1.4-Linux-intel.tar.gz      # -Linux-arm if `uname -m` says aarch64

cd sub_sampled_fielder_vec
#  edit REPO and VENV at the top of cluster/sailfish_sizes.sbatch, and
#  anything you want changed in cluster/series.json, then:
sbatch cluster/sailfish_sizes.sbatch
```

No SLURM: `python scripts/run_series.py --config cluster/series.json`
(same work, sizes in turn). One size only: `--only 8000`. See the plan without
running it: `--dry-run`.

## What it costs

One array task per size, because they are independent and wildly unequal. Extrapolated
from measured single-core times (0.57 s per fill attempt at m=1000, 1.08 s at m=2000,
3.48 s at m=4000) and the project's measured ~30 min/tree sweep at m=6000:

| m | fill (40 trees) | sweep / tree | sweep, 50 workers |
|---|---|---|---|
| 1000 | ~5 min | <1 min | ~1 min |
| 2000 | ~12 min | 1 min | ~1 min |
| 4000 | ~30 min | 9 min | ~9 min |
| 6000 | ~55 min | 30 min | ~30 min |
| 8000 | ~1.2 h | 71 min | ~1.2 h |

**The sweep runs one tree per worker** (since 2026-11-03). With 50 cores and 40 trees per
size, every tree of a size runs at once, so a size costs about one tree's time rather than
forty — m=8000 drops from ~47 h to ~1.2 h. Memory is the limit, not cores: ~2 GB per
concurrent tree at m=8000, so 40 of them want ~80 GB, which is why the job asks for 128G.
Drop `--cpus-per-task` if your node has less.

`--time=24:00:00` is generous for these numbers. If a task is killed,
**resubmit the identical command** — the fill resumes from `manifest.json`'s `last_seed`,
the screen reuses its cached rows, the sweep skips trees already cached. Nothing is
recomputed, and partial results stay valid.

## What to check when it finishes

```bash
python scripts/check_eta_buckets.py      # every bucket at its cap, "no bucket" = 0
python scripts/run_sweep.py --list       # the five datasets and their tree counts
```

A fill that ends `INCOMPLETE` names the short bucket and says which fix it wants: a bucket
whose *window was never reached* needs a wider `--eta-tol`; one *reached but never valid*
needs more attempts. The two are different problems and the message distinguishes them.
Whether eta=1 and eta=15 fill at m=8000 is not known in advance — the eta distribution
shifts with m, and at m=1000 the measured median was 19.8, above every target.

Results land in `results/real_data/runs/<timestamp>/` (curves.csv, recovery_curve.png,
per_tree.csv, eta_buckets.png) with resumable state in `results/real_data/_cache/`.
