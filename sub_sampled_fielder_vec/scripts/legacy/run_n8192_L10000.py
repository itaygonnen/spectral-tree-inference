"""Distance sweep at n=8192, L=10000 (continuation of the n=512..4096 series).

Same SamplingConfig and p grid as run_n512_to_n4096_L10000.py so results
slot into the same family of plots. balanced_binary requires n = 2^k; user
asked for n=8000, rounded up to the next valid power-of-two.
"""
import os, sys, time, traceback
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from src.runners.experiment_runner_utils import (
    extract_config_values,
    generate_run_prefix,
    setup_experiment_directory,
    run_single_experiment,
    auto_generate_plots,
)


SWEEP = {
    "tree_model": "balanced_binary",
    "taxa_values": [8192],
    "sequence_length_values": [10000],
    "mutation_rate": 0.05,
    "bootstrap_reps": 10,
    "num_workers": 4,
    "use_middle_out": False,
    "run_name_prefix": "distance_a1_n8192_L10000",
    "p_values": [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0],
    "tree_params": {"edge_length": 1.0},
    "coherence_k": 4,
    "num_gaps": 0,
    "guardrails_enabled": False,
    "sampling_method": "uniform",
    "matrix_kind": "distance",
    "distance_alpha": 1.0,
}


def main():
    t_total = time.time()
    cfg_vals = extract_config_values(SWEEP, [])
    prefix = generate_run_prefix(SWEEP, cfg_vals["mutation_rate"], cfg_vals["taxa_values"])
    base_dir = setup_experiment_directory(
        SWEEP, prefix, SWEEP,
        script_dir=os.path.dirname(__file__),
    )
    print(f"\n>>> DISTANCE n=8192 sweep: {base_dir}\n", flush=True)

    for n_taxa in cfg_vals["taxa_values"]:
        for seq_len in cfg_vals["sequence_length_values"]:
            t0 = time.time()
            try:
                r = run_single_experiment(
                    SWEEP, base_dir, n_taxa, seq_len,
                    cfg_vals["tree_model"], cfg_vals["mutation_rate"],
                    cfg_vals["p_values"], cfg_vals["bootstrap_reps"],
                    cfg_vals["num_workers"], cfg_vals["use_middle_out"],
                    prefix, cfg_vals["tree_kwargs"],
                )
                dt = time.time() - t0
                print(f"  [distance] n={n_taxa} L={seq_len} done in {dt/60:.1f} min  -> {r['run_dir']}", flush=True)
            except Exception:
                traceback.print_exc()
                print(f"  [distance] n={n_taxa} FAILED after {(time.time()-t0)/60:.1f} min", flush=True)

    try:
        auto_generate_plots(base_dir)
    except Exception:
        traceback.print_exc()

    elapsed = (time.time() - t_total) / 60
    print("\n" + "=" * 80, flush=True)
    print(f"SWEEP DONE in {elapsed:.1f} min", flush=True)
    print(f"  distance -> {base_dir}", flush=True)
    print("=" * 80, flush=True)


if __name__ == "__main__":
    main()
