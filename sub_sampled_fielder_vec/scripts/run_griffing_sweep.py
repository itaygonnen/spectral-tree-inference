"""Griffing-distance-partition sweep driver.

Distance-matrix analogue of figure_3 (Fiedler recovery on synthetic CBM):
    - tree topology generates D  (instead of synthetic flat-CBM S)
    - reference v_ref = leading eigvec of J D J  (instead of closed-form CBM Fiedler)
    - per (p, rep): leading eigvec of J D̂ J, sign-agnostic agreement with v_ref

Outputs two top-level sweep directories under ``results/`` (one per tree
model). Each is consumable by ``plot_griffing_overlay.py``.
"""
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from src.config.presets import custom_config
from src.runners.griffing_sweep import griffing_sweep_for_params


SWEEP_CONFIG: Dict[str, Any] = {
    "taxa_values": [128, 256, 512, 1024, 2048],
    "sequence_length_values": [10000],
    "mutation_rate": 0.1,
    "bootstrap_reps": 5,
    "p_values": [1.0, 0.9, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001],
    # One sweep per tree model. birth_death is the n-comparable baseline;
    # kingman is the headline test (it failed for NJ because tip branches
    # shrink as O(1/n²), but Griffing only needs the dominant split).
    "tree_models": ["birth_death", "kingman"],
    "imputation": "zero",
    "run_name_prefix": "griffing_sweep",
}


def _build_cfg(n: int, L: int, p_values: List[float], tree_model: str):
    kwargs = {}
    if tree_model == "birth_death":
        kwargs.update(birth_rate=1.0, death_rate=0.0)
    elif tree_model == "kingman":
        kwargs.update(pop_size=1.0)
    return custom_config(
        num_taxa=n,
        sequence_length=L,
        mutation_rate=SWEEP_CONFIG["mutation_rate"],
        tree_model=tree_model,
        seq_model="JC69",
        p_values=p_values,
        bootstrap_reps=SWEEP_CONFIG["bootstrap_reps"],
        run_name=SWEEP_CONFIG["run_name_prefix"],
        sampling_method="uniform",
        matrix_kind="distance",
        **kwargs,
    )


def run_for_tree_model(tree_model: str) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    base_dir = (Path(__file__).resolve().parents[1] / "results" / "runs"
                / f"{timestamp}-griffing_sweep_{tree_model}_L"
                  f"{SWEEP_CONFIG['sequence_length_values'][0]}")
    base_dir.mkdir(parents=True, exist_ok=True)
    p_values = SWEEP_CONFIG["p_values"]
    for n in SWEEP_CONFIG["taxa_values"]:
        for L in SWEEP_CONFIG["sequence_length_values"]:
            run_dir = base_dir / f"n{n}_L{L}"
            cfg = _build_cfg(n, L, p_values, tree_model)
            print(f"\n[griffing_sweep] tree={tree_model} n={n} L={L}", flush=True)
            griffing_sweep_for_params(cfg, n, L, str(run_dir),
                                      imputation=SWEEP_CONFIG["imputation"])
    return base_dir


def main():
    written = []
    for tree_model in SWEEP_CONFIG["tree_models"]:
        out = run_for_tree_model(tree_model)
        written.append(out)
    print("\n" + "=" * 80)
    print("Griffing sweeps complete.")
    for d in written:
        print(f"  {d}")
    print("=" * 80)
    for d in written:
        print(f"\nPlot with: python {Path(__file__).parent / 'plot_griffing_overlay.py'} {d}")
    return written


if __name__ == "__main__":
    main()
