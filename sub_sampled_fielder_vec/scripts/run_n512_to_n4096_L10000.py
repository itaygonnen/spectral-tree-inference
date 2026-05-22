"""Distance vs Similarity sweep at n=512..4096, L=10000.

Drives the full ExperimentRunner for both matrix_kind values with otherwise
identical parameters so the two runs share tree + observations (seed=42 in
both) and can be compared head-to-head.

Outputs:
- results/balanced_binary/uniform/<ts>-similarity_n512_4096_L10000/...
- results/balanced_binary/distance_a1p000/uniform/<ts>-distance_a1_n512_4096_L10000/...
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


TAXA   = [512, 1024, 2048, 4096, 8192]
SEQLEN = 10000
MU     = 0.1
REPS   = 10
PVALS  = [0.005, 0.01, 0.02, 0.05, 0.07, 0.085, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50, 0.75, 1.0]  # dense near the p≈0.1 transition


def make_sweep(matrix_kind: str, alpha: float = 1.0):
    prefix_kind = "distance_a1" if matrix_kind == "distance" else "similarity"
    return {
        "tree_model": "balanced_binary",
        "taxa_values": TAXA,
        "sequence_length_values": [SEQLEN],
        "mutation_rate": MU,
        "bootstrap_reps": REPS,
        "num_workers": 4,
        "use_middle_out": False,
        "run_name_prefix": f"{prefix_kind}_n512_4096_L{SEQLEN}",
        "p_values": PVALS,
        "tree_params": {"edge_length": 1.0},
        "coherence_k": 4,
        "num_gaps": 0,
        "guardrails_enabled": False,
        "sampling_method": "uniform",
        "matrix_kind": matrix_kind,
        "distance_alpha": alpha,
    }


def run_sweep(matrix_kind: str):
    SWEEP = make_sweep(matrix_kind)
    cfg_vals = extract_config_values(SWEEP, [])
    prefix = generate_run_prefix(SWEEP, cfg_vals["mutation_rate"], cfg_vals["taxa_values"])
    base_dir = setup_experiment_directory(
        SWEEP, prefix, SWEEP,
        script_dir=os.path.dirname(__file__),
    )
    print(f"\n>>> {matrix_kind.upper()} sweep: {base_dir}\n", flush=True)

    results = []
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
                print(f"  [{matrix_kind}] n={n_taxa} L={seq_len} done in {dt/60:.1f} min  -> {r['run_dir']}", flush=True)
                results.append(r)
            except Exception:
                traceback.print_exc()
                print(f"  [{matrix_kind}] n={n_taxa} FAILED after {(time.time()-t0)/60:.1f} min", flush=True)

    try:
        auto_generate_plots(base_dir)
    except Exception:
        traceback.print_exc()

    return base_dir, results


def main():
    """Run distance sweep only. Pass --with-similarity to also run the baseline."""
    t_total = time.time()
    run_similarity = "--with-similarity" in sys.argv

    print("=" * 80, flush=True)
    print(f"Starting full sweep: taxa={TAXA}, L={SEQLEN}, mu={MU}, reps={REPS}", flush=True)
    print(f"p values: {PVALS}", flush=True)
    print(f"with_similarity: {run_similarity}", flush=True)
    print("=" * 80, flush=True)

    dist_base, _ = run_sweep("distance")
    sim_base = None
    if run_similarity:
        sim_base, _ = run_sweep("similarity")

    elapsed = (time.time() - t_total) / 60
    print("\n" + "=" * 80, flush=True)
    print(f"SWEEPS DONE in {elapsed:.1f} min", flush=True)
    print(f"  distance    -> {dist_base}", flush=True)
    if sim_base:
        print(f"  similarity  -> {sim_base}", flush=True)
    print("=" * 80, flush=True)


if __name__ == "__main__":
    main()
