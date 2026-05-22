"""Classical-NJ-with-subsampling sweep driver.

Maps the three SNJ-paper metrics to a sub-sampling axis, with classical NJ
on the JC distance matrix instead of SNJ on the JC similarity matrix:

  Panel 1: ||D - D̂||_2          (spec-norm error on the distance matrix)
  Panel 2: Q-criterion separation (Fig. 3 analog, NJ flavor)
  Panel 3: Robinson-Foulds       (empirical metric)

For now sampling is plain uniform with no completion (the worst-case
baseline). A future iteration can layer IALM completion on top.
"""
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from src.config.presets import custom_config
from src.runners.nj_sweep import nj_sweep_for_params


SWEEP_CONFIG: Dict[str, Any] = {
    "tree_model": "balanced_binary",
    "seq_model": "JC69",
    "taxa_values": [128, 512, 1024, 2048],
    "sequence_length_values": [2000],
    "mutation_rate": 0.1,
    "bootstrap_reps": 5,
    # Wider spread spanning 4 orders of magnitude. The high-p plateau is
    # boring for the no-completion baseline (cliff at p=1), so we drop the
    # 0.7–0.97 cluster and instead push deeper into low p.
    "p_values": [1.0, 0.9, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001],
    "tree_params": {"edge_length": 1.0},
    "run_name_prefix": "nj_sweep",
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
        matrix_kind="distance",
        edge_length=SWEEP_CONFIG["tree_params"].get("edge_length", 1.0),
    )


def main():
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_dir = Path(__file__).resolve().parents[1] / "results" / "runs" / f"{timestamp}-nj_sweep"
    base_dir.mkdir(parents=True, exist_ok=True)

    p_values = SWEEP_CONFIG["p_values"]
    runs = []
    for n in SWEEP_CONFIG["taxa_values"]:
        for L in SWEEP_CONFIG["sequence_length_values"]:
            run_dir = base_dir / f"n{n}_L{L}"
            cfg = _build_cfg(n, L, p_values)
            print(f"\n[nj_sweep] n={n}, L={L}, p_values={p_values}")
            result = nj_sweep_for_params(cfg, n, L, str(run_dir))
            runs.append({"n": n, "L": L, "run_dir": str(run_dir), "result": result})

    print("\n" + "=" * 80)
    print("NJ sweep complete.")
    for entry in runs:
        n, L = entry["n"], entry["L"]
        print(f"  n={n:>4} L={L:>5} -> {entry['run_dir']}")
    print("=" * 80)
    print(f"\nPlot with: python {Path(__file__).parent / 'plot_nj_three_panel.py'} "
          f"{base_dir}")
    return runs


if __name__ == "__main__":
    main()
