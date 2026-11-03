# Running the size series on a cluster

Five taxon counts — 1000, 2000, 4000, 6000, 8000 — each simulated until all four eta
buckets hold 10 trees, then screened, then swept.

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
#  edit the three paths under CONFIG in cluster/sailfish_sizes.sbatch, then:
sbatch cluster/sailfish_sizes.sbatch
```

No SLURM: `bash cluster/run_sizes.sh` (same work, sizes in turn, foreground).

## What it costs

One array task per size, because they are independent and wildly unequal. Extrapolated
from measured single-core times (0.57 s per fill attempt at m=1000, 1.08 s at m=2000,
3.48 s at m=4000) and the project's measured ~30 min/tree sweep at m=6000:

| m | fill (40 trees) | sweep / tree | sweep (40 trees) |
|---|---|---|---|
| 1000 | ~5 min | <1 min | ~5 min |
| 2000 | ~12 min | 1 min | ~40 min |
| 4000 | ~30 min | 9 min | ~6 h |
| 6000 | ~55 min | 30 min | ~20 h |
| 8000 | ~1.2 h | 71 min | ~47 h |

**~74 h of sweep in total, 90% of it in m=6000 and m=8000.** The fill parallelises across
`--cpus-per-task`; the sweep is one tree at a time, so more cores do not shorten it. The
small sizes finish within the hour and their figures are assemblable while 8000 runs.

`--time=48:00:00` covers m=8000 on these numbers with little margin. If a task is killed,
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
