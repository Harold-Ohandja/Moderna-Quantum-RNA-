"""
Unified Method Comparison.

Runs every method the project implements on the same set of sequences and
writes ONE combined results table, instead of the separate/incompatible
CSVs each track was previously producing on its own:

  1. ViennaRNA MFE            (classical ground truth)
  2. Classical QUBO baseline  (brute-force if small enough, else the
                                greedy heuristic -- see quantum/hybrid_solver.py)
  3. CVaR-VQE                 (primary quantum solver)
  4. QAOA                     (comparison quantum solver)

Quantum methods are skipped above QUANTUM_QUBIT_LIMIT qubits, since local
simulation time grows quickly; the skip itself is logged as a row, since
"where the quantum approach stops being practical" is one of the things
the challenge explicitly asks to characterize.

Usage:
    python scripts/run_full_comparison.py
    python scripts/run_full_comparison.py --sequences GCGCAUACGC AUGCAUGC
"""

import sys
import csv
import time
import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import RNA
from classical.evaluate_energy import evaluate_structure_energy, calculate_energy_gap
from quantum.qubo import build_qubo_matrix
from quantum.decode import decode_bitstring
from quantum.hybrid_solver import solve_qubo_auto
from quantum.vqe_solver import solve_cvar_vqe
from quantum.qaoa_solver import solve_qaoa

DEFAULT_SEQUENCES = [
    "GCGCAUACGC",              # toy sequence: quantum solvers reach the exact MFE
    "AUGCAUGC",                # known formulation-gap case: no solver matches MFE here
    "GCAUCGUAGC",              # second small sequence, small qubit count
    "CUACGGCGCGGCGCCCUUGGCGA", # longer sequence: too many qubits for local quantum sim
]

QUANTUM_QUBIT_LIMIT = 16  # local-simulation practicality cutoff for this comparison
OUT_CSV = ROOT / "output" / "full_comparison_results.csv"

FIELDNAMES = [
    "sequence", "length", "num_qubits", "method", "structure",
    "energy_kcal", "match", "gap_kcal", "runtime_sec", "notes",
]


def compare_sequence(seq: str):
    seq = seq.strip().upper()
    qubo, pairs = build_qubo_matrix(seq)
    num_qubits = len(pairs)
    ref_struct, ref_mfe = RNA.fold(seq)

    rows = []

    # 1. ViennaRNA ground truth
    rows.append({
        "sequence": seq, "length": len(seq), "num_qubits": num_qubits,
        "method": "ViennaRNA_MFE", "structure": ref_struct,
        "energy_kcal": ref_mfe, "match": True, "gap_kcal": 0.0,
        "runtime_sec": None, "notes": "classical ground truth",
    })

    # 2. Classical QUBO baseline (brute-force or heuristic)
    t0 = time.perf_counter()
    bits, cost, cmethod = solve_qubo_auto(qubo, num_qubits, brute_force_limit=12)
    bitstring = "".join(map(str, bits[::-1]))
    structure, _ = decode_bitstring(bitstring, pairs, len(seq))
    runtime = time.perf_counter() - t0
    energy = evaluate_structure_energy(seq, structure)
    gap = calculate_energy_gap(energy, ref_mfe)
    rows.append({
        "sequence": seq, "length": len(seq), "num_qubits": num_qubits,
        "method": f"classical_{cmethod}", "structure": structure,
        "energy_kcal": energy, "match": structure == ref_struct,
        "gap_kcal": gap["absolute_gap_kcal"], "runtime_sec": round(runtime, 4),
        "notes": "not quantum, exact QUBO optimum" if cmethod == "brute_force" else "not quantum, approximate",
    })

    # 3 & 4. Quantum solvers, only within the local-simulation qubit budget
    if num_qubits <= QUANTUM_QUBIT_LIMIT:
        for label, fn, kwargs in [
            ("CVaR-VQE", solve_cvar_vqe, {"alpha": 0.1, "seed": 42}),
            ("QAOA", solve_qaoa, {"seed": 42}),
        ]:
            t0 = time.perf_counter()
            out = fn(qubo, pairs, len(seq), **kwargs)
            runtime = time.perf_counter() - t0
            energy = evaluate_structure_energy(seq, out["structure"])
            gap = calculate_energy_gap(energy, ref_mfe)
            rows.append({
                "sequence": seq, "length": len(seq), "num_qubits": num_qubits,
                "method": label, "structure": out["structure"],
                "energy_kcal": energy, "match": out["structure"] == ref_struct,
                "gap_kcal": gap["absolute_gap_kcal"], "runtime_sec": round(runtime, 4),
                "notes": "",
            })
    else:
        rows.append({
            "sequence": seq, "length": len(seq), "num_qubits": num_qubits,
            "method": "CVaR-VQE / QAOA", "structure": None,
            "energy_kcal": None, "match": None, "gap_kcal": None, "runtime_sec": None,
            "notes": f"skipped: {num_qubits} qubits exceeds local-sim limit ({QUANTUM_QUBIT_LIMIT})",
        })

    return rows


def print_table(rows):
    col_widths = {f: max(len(f), max((len(str(r[f])) for r in rows), default=0)) for f in FIELDNAMES}
    header = " | ".join(f.ljust(col_widths[f]) for f in FIELDNAMES)
    print(header)
    print("-" * len(header))
    for r in rows:
        print(" | ".join(str(r[f]).ljust(col_widths[f]) for f in FIELDNAMES))


def main():
    parser = argparse.ArgumentParser(description="Run the unified method comparison")
    parser.add_argument("--sequences", nargs="+", default=None, help="RNA sequences to compare")
    args = parser.parse_args()

    sequences = args.sequences if args.sequences else DEFAULT_SEQUENCES

    all_rows = []
    for seq in sequences:
        print(f"\n=== {seq} ===")
        rows = compare_sequence(seq)
        all_rows.extend(rows)

    print("\n\n=== Combined results ===")
    print_table(all_rows)

    OUT_CSV.parent.mkdir(exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nSaved combined results to {OUT_CSV}")


if __name__ == "__main__":
    main()
