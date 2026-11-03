#!/bin/bash
# One SLURM task per taxon count: simulate an eta-bucketed dataset, screen it, sweep it.
#
#   sbatch cluster/sailfish_sizes.sbatch
#
# Sizes run as independent array tasks because they are independent jobs -- m=1000
# finishes in minutes while m=8000 takes ~2 days, and nothing should wait on it.
#
# EVERY STAGE IS RESUMABLE. If a task hits the time limit, resubmit the identical
# command: the fill continues from its manifest's last_seed, the screen reuses its
# cached rows, and the sweep skips trees already in its cache. Nothing is recomputed.
#
# Edit the three paths in CONFIG, then submit. Nothing else needs changing.

#SBATCH --job-name=sailfish
#SBATCH --array=0-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/sailfish_%a.out
#SBATCH --error=logs/sailfish_%a.err

set -euo pipefail

# ---- CONFIG -----------------------------------------------------------------
REPO=$HOME/spectral-tree-inference                 # where you cloned it
VENV=$HOME/venvs/str/bin/activate                  # your own python env
IQTREE=$HOME/iqtree-3.1.4-Linux-intel/bin/iqtree3  # prebuilt binary, no compiler

SIZES=(1000 2000 4000 6000 8000)
PER_BUCKET=10          # trees per eta bucket -> 40 trees per size
SEQ_LEN=1000           # alignment length L
MAX_ATTEMPTS=20000     # fill gives up after this; it reports which bucket was short
# -----------------------------------------------------------------------------

N=${SIZES[${SLURM_ARRAY_TASK_ID:-0}]}
WORKERS=${SLURM_CPUS_PER_TASK:-4}
NAME=sailfish_yule_n${N}_L${SEQ_LEN}_g1x4_buckets

cd "$REPO/sub_sampled_fielder_vec"
mkdir -p logs
# shellcheck disable=SC1090
source "$VENV"

echo "=== $N taxa | $WORKERS workers | $(date) ==="
# 1. simulate until every eta bucket is full --------------------------------------
python scripts/generate_dataset.py \
    --n-taxa "$N" --root-seq-size "$SEQ_LEN" \
    --fill-buckets --per-bucket "$PER_BUCKET" \
    --max-attempts "$MAX_ATTEMPTS" --workers "$WORKERS" \
    --iqtree-bin "$IQTREE"

# 2. screen, then 3. sweep. "--stage both" does them in order, resumably.
python scripts/run_sweep.py --dataset "$NAME" --stage both --workers "$WORKERS"

echo "=== $N taxa done | $(date) ==="
echo "results: $REPO/results/real_data/runs/"
