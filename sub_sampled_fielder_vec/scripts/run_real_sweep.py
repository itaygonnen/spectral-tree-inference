#!/usr/bin/env python3
"""Deprecated name for ``scripts/run_sweep.py``.

The runner stopped being specific to real data when the screen and the sweep started
taking a loader instead of a FASTA directory. This shim stays because it is the command
in the cluster runbook and in shell histories; it forwards every argument unchanged.
"""
import runpy
import sys
from pathlib import Path

if __name__ == "__main__":
    print("note: run_real_sweep.py is now scripts/run_sweep.py (same arguments)",
          file=sys.stderr)
    sys.argv[0] = str(Path(__file__).with_name("run_sweep.py"))
    runpy.run_path(sys.argv[0], run_name="__main__")
