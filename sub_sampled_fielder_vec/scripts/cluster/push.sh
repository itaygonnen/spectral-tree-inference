#!/usr/bin/env bash
# Copy the repo (and optionally a cohort from data/cohorts/) to a Linux box over ssh.
# Prefer scripts/cluster/get_data.sh when the data is already in Google Drive -- it pulls
# straight onto the cluster instead of routing GBs through a laptop.
#
#   ./scripts/cluster/push.sh user@host                      # code only
#   ./scripts/cluster/push.sh user@host --data "6000 taxa"   # code + that cohort
#   ./scripts/cluster/push.sh user@host --dest /scratch/itay/str --data all
#
# rsync is resumable: re-run it after a dropped connection and it continues.
# Caches, results, virtualenvs and the git history are never pushed -- the point of the
# remote copy is to *produce* those.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
PKG="sub_sampled_fielder_vec"

HOST=""
DEST='~/spectral-tree-inference'
DATA=""

usage() { sed -n '2,10p' "$0"; exit "${1:-0}"; }

[[ $# -ge 1 ]] || usage 1
HOST="$1"; shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dest) DEST="$2"; shift 2 ;;
    --data) DATA="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done

echo "==> code  ->  $HOST:$DEST"
rsync -avzP --human-readable \
  --exclude '.git/' \
  --exclude '.venv/' \
  --exclude '__pycache__/' \
  --exclude '*.pyc' \
  --exclude "$PKG/cache/" \
  --exclude "$PKG/results/" \
  --exclude 'data/' \
  --exclude "$PKG/docs/overleafs/" \
  "$REPO_ROOT/" "$HOST:$DEST/"

if [[ -n "$DATA" ]]; then
  SRC_ROOT="$REPO_ROOT/data/cohorts"
  if [[ "$DATA" == "all" ]]; then
    echo "==> data (all cohorts, this is GBs)  ->  $HOST:$DEST"
    rsync -avzP --human-readable "$SRC_ROOT/" "$HOST:$DEST/data/cohorts/"
  else
    [[ -d "$SRC_ROOT/$DATA" ]] || { echo "no cohort '$DATA' under $SRC_ROOT" >&2; exit 1; }
    echo "==> data ('$DATA')  ->  $HOST:$DEST"
    ssh "$HOST" "mkdir -p '$DEST/data/cohorts'"
    rsync -avzP --human-readable "$SRC_ROOT/$DATA" "$HOST:$DEST/data/cohorts/"
  fi
fi

cat <<EOF

done. Next, on the cluster:

  ssh $HOST
  cd $DEST && bash $PKG/scripts/cluster/bootstrap.sh
  source .venv/bin/activate
  cd $PKG && python scripts/run_real_sweep.py --list
EOF
