# classical/benchmark.py
"""
Classical ViennaRNA Benchmark Runner.
Can be run two ways:
  python -m classical.benchmark --sequence "GCGCAUACGC"
  python classical/benchmark.py --sequence "GCGCAUACGC"
Both work: this file adds the repo root to sys.path itself, so it doesn't
depend on how it's invoked.
"""

import sys
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import ViennaRNA as RNA
from classical.evaluate_energy import evaluate_structure_energy
from quantum.qubo import build_qubo_matrix


def run_benchmark_for(sequence: str) -> dict:
    # 1) Classical MFE (ViennaRNA)
    structure, mfe_energy = RNA.fold(sequence)

    # 2) Energy via our evaluator
    eval_energy = evaluate_structure_energy(sequence, structure)

    # 3) QUBO formulation
    qubo, candidate_pairs = build_qubo_matrix(sequence)

    # 4) Print summary
    print("=== Benchmark ===")
    print("Sequence       :", sequence)
    print("Length         :", len(sequence))
    print("MFE structure  :", structure)
    print("MFE energy     :", mfe_energy, "kcal/mol")
    print("Eval energy    :", eval_energy, "kcal/mol")
    print("Candidate pairs:", candidate_pairs)
    print("Num pairs      :", len(candidate_pairs))
    print("QUBO size      :", len(qubo), "non-zero terms")
    print()

    # 5) Return structured result
    return {
        "sequence": sequence,
        "length": len(sequence),
        "mfe_structure": structure,
        "mfe_energy": float(mfe_energy),
        "eval_energy": float(eval_energy),
        "candidate_pairs": candidate_pairs,
        "num_pairs": len(candidate_pairs),
        "qubo_size": len(qubo),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the classical ViennaRNA benchmark")
    parser.add_argument(
        "--sequence", type=str, default=None,
        help="Single RNA sequence to benchmark. If omitted, runs the default demo pair.",
    )
    args = parser.parse_args()

    sequences = [args.sequence.strip().upper()] if args.sequence else ["GCGCAUACGC", "AUGCAUGC"]
    results = {}

    for seq in sequences:
        results[seq] = run_benchmark_for(seq)

    # Ensure data/benchmark exists
    out_dir = ROOT / "data" / "benchmark"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / "reference_mfe.json"
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"Saved reference MFE benchmarks to {out_path}")