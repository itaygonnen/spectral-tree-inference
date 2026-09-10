#!/usr/bin/env bash
# Build the runtime on a fresh Linux box. Run it from the REPO ROOT (the directory that
# holds setup.py), once per machine:
#
#   cd ~/spectral-tree-inference && bash sub_sampled_fielder_vec/scripts/cluster/bootstrap.sh
#
# Creates .venv, installs requirements-cluster.txt, installs spectraltree with --no-deps
# (setup.py pulls oct2py/toytree/sphinx, which are not needed and need system Octave),
# then verifies the imports the experiments actually use.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$REPO_ROOT"

PY="${PYTHON:-}"
if [[ -z "$PY" ]]; then
  for cand in python3.12 python3.11 python3.10 python3; do
    command -v "$cand" >/dev/null 2>&1 && { PY="$cand"; break; }
  done
fi
[[ -n "$PY" ]] || { echo "no python3 found; module load python, or set PYTHON=..." >&2; exit 1; }

echo "==> python: $($PY -V) at $(command -v "$PY")"
"$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' || {
  echo "need python >= 3.10 (pydantic v2, PEP 604 types). Try: module avail python" >&2
  exit 1
}

if [[ ! -d .venv ]]; then
  echo "==> creating .venv"
  "$PY" -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

echo "==> installing dependencies"
python -m pip install --upgrade pip wheel
python -m pip install -r requirements-cluster.txt
if ! python -m pip install -e . --no-deps; then
  cat >&2 <<'EOF'

Installing the package failed. If the error mentions a missing module (sphinx, say),
setup.py wants a build-time import this branch should already have made optional --
make sure the checkout is up to date (git pull), then re-run this script.
EOF
  exit 1
fi

echo "==> verifying"
python - <<'PY'
import importlib, sys
from pathlib import Path

root = Path.cwd()
sys.path.insert(0, str(root / "sub_sampled_fielder_vec"))
sys.path.insert(0, str(root))

for mod in ("numpy", "scipy", "pandas", "sklearn", "dendropy", "igraph",
            "pydantic", "matplotlib", "spectraltree"):
    importlib.import_module(mod)
    print(f"  ok  {mod}")

# the two entry points, imported for real: this is where a missing optional dependency
# (toytree, PIL) or a stale call site shows up, not in the module list above
importlib.import_module("analysis.utils.real_interactive")
print("  ok  analysis.utils.real_interactive (the launcher's real-data branch)")

from analysis.utils.real_cohorts import list_cohorts
cohorts = list_cohorts()
print(f"  ok  analysis.utils.real_cohorts -> {len(cohorts)} cohort(s)")
for c in cohorts:
    m, L = c.shape()
    print(f"      {c.name!r}: {len(c.ids())} trees, m={m}, L={L}")
if not cohorts:
    print("      (no cohorts yet -- fetch one:")
    print("       bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh '6000 taxa')")
PY

cat <<'EOF'

ready. Every session:

  source .venv/bin/activate
  cd sub_sampled_fielder_vec

Interactive (needs the ssh session to stay open -- use tmux):
  python scripts/interactive_run.py        # "d" = real data

Data (once per cohort, straight from Google Drive):
  bash sub_sampled_fielder_vec/scripts/cluster/get_data.sh "6000 taxa" --shared

Long jobs (survive logout):
  nohup python scripts/run_real_sweep.py --cohort "6000 taxa" --stage sweep \
      --cohort-rule valid_S > logs/real_sweep.log 2>&1 &
  tail -f logs/real_sweep.log
EOF
