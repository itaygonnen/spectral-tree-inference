"""SNJ-with-subsampling sweep driver.

Maps the SNJ paper's three metrics to a sub-sampling axis:

  Panel 1: ||R - R̂||_2          (Thm 4.2)
  Panel 2: σ₂ separation         (Fig. 3 analog)
  Panel 3: Robinson-Foulds       (paper's empirical metric)

For now sampling is plain uniform with no completion (the worst-case
baseline). A future iteration can layer IALM completion on top.
"""
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import numpy as np

from src.config.presets import custom_config
from src.runners.snj_sweep import snj_sweep_for_params


SWEEP_CONFIG: Dict[str, Any] = {
    "tree_model": "balanced_binary",
    "seq_model": "JC69",
    # Multi-N sweep; each n gets the same L. At L=2000 with mu=0.1,
    # SNJ on the full R reliably recovers the tree up to n~2048.
    "taxa_values": [64, 128, 256, 512, 1024],
    "sequence_length_values": [2000],
    "mutation_rate": 0.1,
    "bootstrap_reps": 5,
    # Denser grid below the transition (p~0.8-0.95 is where success rate
    # collapses for the no-completion baseline). Logspace + targeted points.
    "p_values": [1.0, 0.97, 0.95, 0.93, 0.90, 0.85, 0.80, 0.70, 0.50, 0.30, 0.10],
    "tree_params": {"edge_length": 1.0},
    "run_name_prefix": "snj_sweep",
    # SNJ knobs
    "snj_bifurcating": False,  # match unrooted ground-truth convention
    "snj_alpha": 1.0,
}


def _build_cfg(n: int, L: int, p_values: List[float]):
    return custom_config(
        num_taxa=n,
        sequence_length=L,
        mutation_rate=SWEEP_CONFIG["mutation_rate"],
        tree_model=SWEEP_CONFIG["tree_model"],
        seq_model=SWEEP_CONFIG["seq_model"],
        p_values=p_values,
        bootstrap_reps=SWEEP_CONFIG["bootstrap_reps"],
        run_name=SWEEP_CONFIG["run_name_prefix"],
        sampling_method="uniform",
        matrix_kind="similarity",
        edge_length=SWEEP_CONFIG["tree_params"].get("edge_length", 1.0),
    )


def main():
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_dir = Path(__file__).resolve().parents[1] / "results" / "runs" / f"{timestamp}-snj_sweep"
    base_dir.mkdir(parents=True, exist_ok=True)

    p_values = SWEEP_CONFIG["p_values"]
    runs = []
    for n in SWEEP_CONFIG["taxa_values"]:
        for L in SWEEP_CONFIG["sequence_length_values"]:
            run_dir = base_dir / f"n{n}_L{L}"
            cfg = _build_cfg(n, L, p_values)
            print(f"\n[snj_sweep] n={n}, L={L}, p_values={p_values}")
            result = snj_sweep_for_params(
                cfg, n, L, str(run_dir),
                bifurcating=bool(SWEEP_CONFIG["snj_bifurcating"]),
                snj_alpha=float(SWEEP_CONFIG["snj_alpha"]),
            )
            runs.append({"n": n, "L": L, "run_dir": str(run_dir), "result": result})

    print("\n" + "=" * 80)
    print("SNJ sweep complete.")
    for entry in runs:
        n, L = entry["n"], entry["L"]
        print(f"  n={n:>4} L={L:>5} -> {entry['run_dir']}")
    print("=" * 80)
    print(f"\nPlot with: python {Path(__file__).parent / 'plot_snj_three_panel.py'} "
          f"{base_dir}")
    return runs


if __name__ == "__main__":
    main()
