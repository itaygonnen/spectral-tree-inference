"""Distance-matrix sub-sampling experiment.

Mirrors ``run_experiment.py`` but sets ``matrix_kind="distance"`` so the
pipeline sub-samples the paralinear distance matrix D and post-transforms
to S = exp(-α · D̂) before the Fiedler / σ₂ step. See
``docs/papers/SNJ_Jaffe_Kluger.pdf`` and ``docs/papers/README.md`` for the
SNJ α parameter (``M^α = exp(-α D)``) that this experiment exercises.

Sub-sampling is non-commutative with the exp() kernel under IPW debiasing,
so sub-sampling D then exponentiating produces a different empirical
σ₂/partition than sub-sampling the similarity directly — the question this
experiment quantifies.
"""
import os
import sys
from typing import Dict, Any

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np

from src.runners.experiment_runner_utils import (
    extract_config_values,
    generate_run_prefix,
    setup_experiment_directory,
    run_single_experiment,
    auto_generate_plots,
)


WIDE_SWEEP_P_VALUES = list(np.logspace(-4, 0, 20))

SWEEP_CONFIG: Dict[str, Any] = {
    "tree_model": "balanced_binary",
    "taxa_values": [256],
    "sequence_length_values": [400],
    "mutation_rate": 0.05,
    "bootstrap_reps": 10,
    "num_workers": 4,
    "use_middle_out": False,
    "run_name_prefix": "distance_alpha_sweep",
    "p_values": [0.05, 0.1, 0.25, 0.5, 0.75, 1.0],
    "tree_params": {"edge_length": 1.0},
    "coherence_k": 4,
    "num_gaps": 0,
    "guardrails_enabled": False,

    # Sampling
    "sampling_method": "uniform",

    # Distance-matrix knobs (the point of this script)
    "matrix_kind": "distance",
    "distance_alpha": 1.0,
}


def main():
    config = SWEEP_CONFIG.copy()
    cfg_vals = extract_config_values(config, WIDE_SWEEP_P_VALUES)
    prefix = generate_run_prefix(config, cfg_vals["mutation_rate"], cfg_vals["taxa_values"])
    base_dir = setup_experiment_directory(
        config, prefix, SWEEP_CONFIG,
        script_dir=os.path.dirname(__file__),
    )

    multi_run_results = []
    for n_taxa in cfg_vals["taxa_values"]:
        for seq_len in cfg_vals["sequence_length_values"]:
            result = run_single_experiment(
                config, base_dir, n_taxa, seq_len,
                cfg_vals["tree_model"], cfg_vals["mutation_rate"],
                cfg_vals["p_values"], cfg_vals["bootstrap_reps"],
                cfg_vals["num_workers"], cfg_vals["use_middle_out"],
                prefix, cfg_vals["tree_kwargs"],
            )
            multi_run_results.append(result)

    print(f"\n{'='*80}")
    print(f"Distance-matrix sweep done ({cfg_vals['tree_model']}, α={config['distance_alpha']}).")
    for entry in multi_run_results:
        print(f"n={entry['num_taxa']:>4}, L={entry['sequence_length']:>5} -> {entry['run_dir']}")
    print(f"{'='*80}")

    auto_generate_plots(base_dir)
    return multi_run_results


if __name__ == "__main__":
    main()
