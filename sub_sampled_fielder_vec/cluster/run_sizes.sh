#!/bin/bash
# The same job without SLURM: every size in turn, in the foreground.
#
#   bash cluster/run_sizes.sh              # all five sizes
#   bash cluster/run_sizes.sh 1000 2000    # just these
#
# Prefer sailfish_sizes.sbatch on a cluster -- it runs the sizes in parallel and
# survives a time limit. This exists for a single machine, or to try one size first.
# Resumable either way: re-run the same command and nothing is recomputed.

set -euo pipefail

REPO=${REPO:-$HOME/spectral-tree-inference}
IQTREE=${IQTREE:-iqtree3}                 # a name on PATH, or an absolute path
PER_BUCKET=${PER_BUCKET:-10}
SEQ_LEN=${SEQ_LEN:-1000}
WORKERS=${WORKERS:-4}
MAX_ATTEMPTS=${MAX_ATTEMPTS:-20000}

SIZES=("$@")
[ ${#SIZES[@]} -eq 0 ] && SIZES=(1000 2000 4000 6000 8000)

cd "$REPO/sub_sampled_fielder_vec"
for N in "${SIZES[@]}"; do
    NAME=sailfish_yule_n${N}_L${SEQ_LEN}_g1x4_buckets
    echo "=== $N taxa | $(date) ==="
    python scripts/generate_dataset.py \
        --n-taxa "$N" --root-seq-size "$SEQ_LEN" \
        --fill-buckets --per-bucket "$PER_BUCKET" \
        --max-attempts "$MAX_ATTEMPTS" --workers "$WORKERS" \
        --iqtree-bin "$IQTREE"
    python scripts/run_sweep.py --dataset "$NAME" --stage both --workers "$WORKERS"
done
echo "all sizes done -- results under $REPO/results/real_data/runs/"
