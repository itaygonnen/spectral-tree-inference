"""Build a grid comparison plot from the latest distance + similarity full sweeps.

Walks ``results/balanced_binary/`` looking for the most recent run dirs whose
prefix matches the ``run_n512_to_n4096_L10000.py`` driver, loads each
``results_grid_merged.json``, and plots sign_agreement / dot_product / partition_agreement
vs p with one curve per n_taxa, colored by matrix_kind.
"""
import os, sys, json, glob
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1] / "results" / "runs" / "balanced_binary"


def latest_run(subdir: str, prefix: str) -> Path:
    """Find the most recent run dir under ROOT/<subdir>/uniform/ matching prefix."""
    candidates = sorted((ROOT / subdir / "uniform").glob(f"*-{prefix}*"))
    if not candidates:
        return None
    return candidates[-1]


def load_merged(run_dir: Path):
    path = run_dir / "results_grid_merged.json"
    with open(path) as f:
        return json.load(f)


def main():
    dist_dir = latest_run("distance_a1p000", "distance_a1_n512_4096_L10000")
    sim_dir  = latest_run("uniform", "similarity_n512_4096_L10000")
    if not dist_dir or not sim_dir:
        print(f"[!] missing run dirs: distance={dist_dir}, similarity={sim_dir}")
        return 1

    print(f"distance:   {dist_dir}")
    print(f"similarity: {sim_dir}")

    dist = load_merged(dist_dir)
    sim  = load_merged(sim_dir)

    # Group rows by n_taxa
    def by_n(payload):
        out = {}
        for row in payload["rows"]:
            out.setdefault(int(row["num_taxa"]), []).append(row)
        for k in out:
            out[k].sort(key=lambda r: r["p"])
        return out

    dist_by_n = by_n(dist)
    sim_by_n  = by_n(sim)
    ns = sorted(set(dist_by_n) | set(sim_by_n))

    metrics = [
        ("sign_agreement",        "sign agreement (%)", "Sign agreement vs reference Fiedler"),
        ("dot_product",           "dot product",        "Fiedler dot-product vs reference"),
        ("partition_agreement_M", "partition agr (%)",  "Partition agreement on M"),
    ]

    fig, axes = plt.subplots(len(ns), len(metrics), figsize=(5 * len(metrics), 3.5 * len(ns)), sharex=True)
    if len(ns) == 1:
        axes = axes[None, :]

    for r, n in enumerate(ns):
        for c, (key, ylabel, title) in enumerate(metrics):
            ax = axes[r, c]
            if n in sim_by_n:
                xs = [row["p"] for row in sim_by_n[n]]
                ys = [row.get(key) if row.get(key) is not None else float("nan") for row in sim_by_n[n]]
                ax.plot(xs, ys, "-o", label="similarity", color="C0")
            if n in dist_by_n:
                xs = [row["p"] for row in dist_by_n[n]]
                ys = [row.get(key) if row.get(key) is not None else float("nan") for row in dist_by_n[n]]
                ax.plot(xs, ys, "--s", label="distance α=1", color="C1")
            ax.set_xscale("log")
            ax.grid(alpha=0.3)
            if r == 0:
                ax.set_title(title)
            if c == 0:
                ax.set_ylabel(f"n={n}\n\n{ylabel}")
            else:
                ax.set_ylabel(ylabel)
            if r == len(ns) - 1:
                ax.set_xlabel("sampling rate p")
            if r == 0 and c == 0:
                ax.legend()

    fig.suptitle("Balanced binary  L=10000, μ=0.05, 10 reps  —  distance α=1 vs similarity", fontsize=13)
    fig.tight_layout()
    out = ROOT / "comparison_distance_vs_similarity_n512_4096.png"
    fig.savefig(out, dpi=140)
    print(f"\n[saved] {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
