"""Run the 600-tree real-dataset B-matrix clan-partition benchmark.

Usage
-----
    python scripts/run_real_data_benchmark.py \\
        --data-root /path/to/drive/data \\
        --out-dir results/runs/real_data_benchmark \\
        [--p-values 0.1 0.2 0.5 1.0] \\
        [--bootstrap-reps 5] \\
        [--max-trees 10]

Data root is expected to contain:
    <data_root>/200 taxa/newick/random_tree_1.nwk ...
    <data_root>/200 taxa/fasta/random_tree_1.fasta ...
    <data_root>/1000 taxa/newick/random_tree_1.nwk ...
    <data_root>/1000 taxa/fasta/...

Results land in <out-dir>/n200/ and <out-dir>/n1000/.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.runners.real_data_bpart import run_real_data_benchmark


_SIZE_CLASSES = ["200 taxa", "1000 taxa"]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Real-data B-matrix clan-partition benchmark"
    )
    parser.add_argument("--data-root", required=True,
                        help="Root folder with '200 taxa/' and '1000 taxa/' subdirs")
    parser.add_argument("--out-dir", default="results/runs/real_data_benchmark",
                        help="Output directory (default: results/runs/real_data_benchmark)")
    parser.add_argument("--p-values", nargs="+", type=float,
                        default=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                        help="Sub-sampling fractions to sweep")
    parser.add_argument("--bootstrap-reps", type=int, default=5,
                        help="Bootstrap replicates per p value per tree (default: 5)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-trees", type=int, default=None,
                        help="Cap number of trees (for quick validation runs)")
    parser.add_argument("--size-class", choices=["200", "1000", "both"],
                        default="both", help="Which size class to run (default: both)")
    parser.add_argument("--fasta-ext", default=".fasta",
                        help="FASTA file extension (default: .fasta)")
    args = parser.parse_args()

    data_root = Path(args.data_root)
    out_dir = Path(args.out_dir)
    size_classes = (
        _SIZE_CLASSES if args.size_class == "both"
        else [f"{args.size_class} taxa"]
    )

    for size in size_classes:
        n_str = size.split()[0]
        newick_dir = data_root / size / "newick"
        fasta_dir = data_root / size / "fasta"
        run_dir = out_dir / f"n{n_str}"

        if not newick_dir.exists():
            print(f"[SKIP] {newick_dir} not found — download the dataset first.")
            continue
        if not fasta_dir.exists():
            print(f"[SKIP] {fasta_dir} not found — check dataset structure.")
            continue

        print(f"\n=== Size class: {size} ===")
        run_real_data_benchmark(
            newick_dir=str(newick_dir),
            fasta_dir=str(fasta_dir),
            run_dir=str(run_dir),
            p_values=args.p_values,
            bootstrap_reps=args.bootstrap_reps,
            seed=args.seed,
            max_trees=args.max_trees,
            fasta_ext=args.fasta_ext,
        )
        print(f"Results → {run_dir}/aggregate_results.json")


if __name__ == "__main__":
    main()
