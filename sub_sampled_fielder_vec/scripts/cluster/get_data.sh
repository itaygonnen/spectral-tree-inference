#!/usr/bin/env bash
# Fetch a cohort into <repo>/data/cohorts/ from Google Drive, straight onto the machine
# that will compute on it. No laptop in the middle, resumable, works from Linux, macOS or
# Windows (Git Bash / WSL).
#
#   bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh "6000 taxa"
#   bash .../get_data.sh "6000 taxa" --remote gdrive --shared
#
# First run only: authorise rclone (see --help output below). Everything after that is one
# command per cohort.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
DEST_ROOT="${STR_DATA_DIR:-$REPO_ROOT/data/cohorts}"

REMOTE="gdrive"
DRIVE_PATH=""
SHARED=0
COHORT=""

help() {
  cat <<'EOF'
usage: get_data.sh <cohort name> [--remote NAME] [--path "Drive/sub/dir"] [--shared]

  <cohort name>   directory to create under data/cohorts/, e.g. "6000 taxa".
                  Also the folder name looked for in Drive unless --path is given.
  --remote NAME   rclone remote (default: gdrive)
  --path P        folder inside the remote, if it is not just <cohort name>
  --shared        the folder was shared with you and is not in My Drive
                  (adds --drive-shared-with-me)

one-time setup, on the machine that will hold the data:

  1. install rclone:      curl https://rclone.org/install.sh | sudo bash
                          (no sudo? download the static binary from rclone.org/downloads
                           and put it on your PATH -- it is a single file, any OS)
  2. rclone config
       n) new remote   name: gdrive   storage: drive
       client_id / client_secret: blank      scope: 1 (or 2 for read-only)
       "Use web browser to automatically authenticate?"  -> N on a headless cluster
       it then prints a command to run on ANY machine that has a browser
       (Windows, Mac, Linux -- rclone must be installed there too):
           rclone authorize "drive"
       log in, copy the token it prints, paste it back into the cluster prompt.

the folder must be visible to the Google account you authorise with -- a folder shared
with that account is fine (use --shared).
EOF
}

[[ $# -ge 1 ]] || { help; exit 1; }
case "$1" in -h|--help) help; exit 0 ;; esac
COHORT="$1"; shift
while [[ $# -gt 0 ]]; do
  case "$1" in
    --remote) REMOTE="$2"; shift 2 ;;
    --path)   DRIVE_PATH="$2"; shift 2 ;;
    --shared) SHARED=1; shift ;;
    -h|--help) help; exit 0 ;;
    *) echo "unknown argument: $1" >&2; help; exit 1 ;;
  esac
done
: "${DRIVE_PATH:=$COHORT}"

command -v rclone >/dev/null 2>&1 || {
  echo "rclone not found -- see 'get_data.sh --help' for the one-time setup" >&2; exit 1; }
rclone listremotes | grep -qx "${REMOTE}:" || {
  echo "no rclone remote called '${REMOTE}'. Run: rclone config" >&2; exit 1; }

FLAGS=(--progress --transfers 8 --checkers 16 --drive-acknowledge-abuse)
[[ $SHARED -eq 1 ]] && FLAGS+=(--drive-shared-with-me)

DEST="$DEST_ROOT/$COHORT"
mkdir -p "$DEST"
echo "==> ${REMOTE}:${DRIVE_PATH}  ->  $DEST"
rclone copy "${REMOTE}:${DRIVE_PATH}" "$DEST" "${FLAGS[@]}"

n_fa=$(find "$DEST" -name '*.fasta' | wc -l | tr -d ' ')
n_nw=$(find "$DEST" -name '*.nwk' -o -name '*.newick' | wc -l | tr -d ' ')
echo "==> $n_fa alignments, $n_nw trees under $DEST"

if [[ ! -d "$DEST/fasta" || ! -d "$DEST/newick" ]]; then
  cat >&2 <<EOF

NOTE: a cohort must look like
    $DEST/fasta/*.fasta
    $DEST/newick/*.nwk
Move what you downloaded into those two subdirectories (or lay the Drive folder out that
way once, and every future fetch is correct). Verify with:
    python sub_sampled_fielder_vec/scripts/run_real_sweep.py --list
EOF
  exit 2
fi
echo "ready:  python sub_sampled_fielder_vec/scripts/run_real_sweep.py --list"
