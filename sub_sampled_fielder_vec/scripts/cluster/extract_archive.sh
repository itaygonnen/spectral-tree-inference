#!/usr/bin/env bash
# Unpack the dataset file (sim_trees_*.tar.gz) into a cohort this code can read.
#
#   bash scripts/cluster/extract_archive.sh --name "6000 taxa"
#   bash scripts/cluster/extract_archive.sh <path to .tar.gz> --name "6000 taxa" \
#        --pattern 'random_tree_0[0-4]*.fasta'
#
# With no path it takes the single .tar.gz sitting in the repo's data/ folder, which is
# where the dataset file is meant to be dropped -- the unpacked cohort lands beside it.
#
# The file holds every alignment (~90 GB unpacked) beside its true tree (~600 MB), laid
# out as <dataset>/MSAs/*.fasta and <dataset>/trees/*.nwk. This writes ALL the trees and
# only the alignments you ask for, into the folders the experiments expect:
#
#   data/cohorts/<name>/fasta/random_tree_0001.fasta   ...
#   data/cohorts/<name>/newick/random_tree_0001.nwk    ...
#
# Default: alignments 0001-0100, about 3 GB at m=6000. It reads the whole compressed file
# once either way -- tar has to decompress the entire stream to find the members -- so
# expect a few minutes, and run it once with the pattern you actually want.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
DEST_ROOT="${STR_DATA_DIR:-$REPO_ROOT/data/cohorts}"

usage() { sed -n '2,24p' "$0"; exit "${1:-0}"; }
case "${1:-}" in -h|--help) usage 0 ;; esac

DATASET_FILE=""
if [[ $# -ge 1 && "$1" != --* ]]; then
  DATASET_FILE="$1"; shift
fi
NAME=""
PATTERNS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --name)    NAME="$2"; shift 2 ;;
    --pattern) PATTERNS+=("$2"); shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done

DATA_DIR="$REPO_ROOT/data"
if [[ -z "$DATASET_FILE" ]]; then
  # no path given: the one .tar.gz in data/
  shopt -s nullglob
  found=("$DATA_DIR"/*.tar.gz)
  shopt -u nullglob
  case ${#found[@]} in
    0) echo "no .tar.gz in $DATA_DIR -- copy the dataset file there, or pass its path" >&2
       exit 1 ;;
    1) DATASET_FILE="${found[0]}" ;;
    *) echo "several .tar.gz files in $DATA_DIR; name the one you want:" >&2
       printf '  %s\n' "${found[@]}" >&2
       exit 1 ;;
  esac
elif [[ ! -f "$DATASET_FILE" && -f "$DATA_DIR/$DATASET_FILE" ]]; then
  DATASET_FILE="$DATA_DIR/$DATASET_FILE"      # a bare filename, relative to data/
fi
[[ -f "$DATASET_FILE" ]] || { echo "no such file: $DATASET_FILE" >&2; exit 1; }
echo "==> dataset file : $DATASET_FILE"
if [[ ${#PATTERNS[@]} -eq 0 ]]; then
  PATTERNS=('random_tree_00[0-9][0-9].fasta' 'random_tree_0100.fasta')   # 0001-0100
fi

# The archive's own top-level directory name, e.g. sim_trees_3000x6000sp_5k_JC_nohet_noindels
TOP="$(tar tzf "$DATASET_FILE" | head -1 | cut -d/ -f1)"
[[ -n "$NAME" ]] || NAME="$TOP"
DEST="$DEST_ROOT/$NAME"
STAGE="$DEST_ROOT/.staging_$$"

echo "==> dataset root : $TOP"
echo "==> cohort       : $DEST"
echo "==> alignments   : ${PATTERNS[*]}"
echo "    (reads the whole compressed file once; a few minutes)"

mkdir -p "$STAGE"
MEMBERS=("$TOP/trees")
for pat in "${PATTERNS[@]}"; do MEMBERS+=("$TOP/MSAs/$pat"); done
tar -xzf "$DATASET_FILE" -C "$STAGE" "${MEMBERS[@]}"

mkdir -p "$DEST"
[[ -d "$STAGE/$TOP/trees" ]] && mv "$STAGE/$TOP/trees" "$DEST/newick"
[[ -d "$STAGE/$TOP/MSAs" ]] && mv "$STAGE/$TOP/MSAs" "$DEST/fasta"
rmdir "$STAGE/$TOP" "$STAGE" 2>/dev/null || true

n_fa=$(find "$DEST/fasta" -name '*.fasta' 2>/dev/null | wc -l | tr -d ' ')
n_nw=$(find "$DEST/newick" -name '*.nwk' 2>/dev/null | wc -l | tr -d ' ')
echo "==> $n_fa alignments, $n_nw trees in $DEST"
echo
echo "check it is visible to the code:"
echo "  python sub_sampled_fielder_vec/scripts/run_real_sweep.py --list"
