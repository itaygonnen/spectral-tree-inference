#!/usr/bin/env bash
# Turn a sim_trees_*.tar.gz archive into a cohort this code can read.
#
#   bash scripts/cluster/extract_archive.sh sim_trees_3000x6000sp_5k_JC_nohet_noindels.tar.gz
#   bash scripts/cluster/extract_archive.sh <tarball> --name "6000 taxa" --pattern 'random_tree_0[0-4]*.fasta'
#
# The archive holds every alignment (~90 GB unpacked) beside its true tree (~600 MB), laid
# out as <archive>/MSAs/*.fasta and <archive>/trees/*.nwk. This extracts ALL the trees and
# only the alignments you ask for, into the layout the experiments expect:
#
#   data/cohorts/<name>/fasta/random_tree_0001.fasta   ...
#   data/cohorts/<name>/newick/random_tree_0001.nwk    ...
#
# Default: alignments 0001-0100, about 3 GB at m=6000. One pass over the archive either
# way -- tar has to decompress the whole stream to find the members, so expect a few
# minutes and run it once with the pattern you actually want.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
DEST_ROOT="${STR_DATA_DIR:-$REPO_ROOT/data/cohorts}"

usage() { sed -n '2,20p' "$0"; exit "${1:-0}"; }
[[ $# -ge 1 ]] || usage 1
case "$1" in -h|--help) usage 0 ;; esac

TARBALL="$1"; shift
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

[[ -f "$TARBALL" ]] || { echo "no such file: $TARBALL" >&2; exit 1; }
if [[ ${#PATTERNS[@]} -eq 0 ]]; then
  PATTERNS=('random_tree_00[0-9][0-9].fasta' 'random_tree_0100.fasta')   # 0001-0100
fi

# The archive's own top-level directory name, e.g. sim_trees_3000x6000sp_5k_JC_nohet_noindels
TOP="$(tar tzf "$TARBALL" | head -1 | cut -d/ -f1)"
[[ -n "$NAME" ]] || NAME="$TOP"
DEST="$DEST_ROOT/$NAME"
STAGE="$DEST_ROOT/.staging_$$"

echo "==> archive root : $TOP"
echo "==> cohort       : $DEST"
echo "==> alignments   : ${PATTERNS[*]}"
echo "    (one pass over the archive; a few minutes)"

mkdir -p "$STAGE"
MEMBERS=("$TOP/trees")
for pat in "${PATTERNS[@]}"; do MEMBERS+=("$TOP/MSAs/$pat"); done
tar -xzf "$TARBALL" -C "$STAGE" "${MEMBERS[@]}"

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
