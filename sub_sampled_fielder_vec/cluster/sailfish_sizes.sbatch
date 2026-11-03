#!/bin/bash
# The whole size series, configured by cluster/series.json. One array task per size.
#
#   sbatch cluster/sailfish_sizes.sbatch
#
# Sizes are independent and unequal -- m=1000 finishes in minutes, m=8000 takes hours --
# so nothing waits on the slowest. EVERY STAGE RESUMES: if a task hits the time limit,
# resubmit the identical command and nothing is recomputed.
#
# Edit the two paths below and everything else in cluster/series.json.

#SBATCH --job-name=sailfish
#SBATCH --array=0-4
#SBATCH --cpus-per-task=50
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sailfish_%a.out
#SBATCH --error=logs/sailfish_%a.err

set -euo pipefail

REPO=$HOME/spectral-tree-inference        # where you cloned it
VENV=$HOME/venvs/str/bin/activate         # your own python env

cd "$REPO/sub_sampled_fielder_vec"
mkdir -p logs
# shellcheck disable=SC1090
source "$VENV"

# The size this array task owns, read from the same JSON that configures everything else.
SIZE=$(python -c "import json,sys; print(json.load(open('cluster/series.json'))['sizes'][${SLURM_ARRAY_TASK_ID:-0}])")
WORKERS=${SLURM_CPUS_PER_TASK:-4}

echo "=== $SIZE taxa | $WORKERS workers | $(date) ==="
python scripts/run_series.py --config cluster/series.json --only "$SIZE"
echo "=== $SIZE taxa done | $(date) ==="
