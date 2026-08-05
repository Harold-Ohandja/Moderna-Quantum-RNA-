import csv
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from quantum.qubo import build_qubo_matrix
from quantum.decode import decode_bitstring
from classical.evaluate_energy import evaluate_structure_energy
from quantum.hybrid_solver import solve_qubo_auto


OUT_CSV = ROOT / "output" / "rna_solver_results_milestone5.csv"

SEQUENCES = [
    "GCGCAUACGC",
    "GCGCGUACGC",
    "GCAUCGUAGC",
    "CUACGGCGCGGCGCCCUUGGCGA",
    "GCGCGCGAAAUUCGCGCG",
]


def run_one(seq: str, brute_force_limit: int = 12):
    start = time.perf_counter()
    qubo, pairs = build_qubo_matrix(seq)
    bits, cost, method = solve_qubo_auto(qubo, len(pairs), brute_force_limit=brute_force_limit)
    bitstring = "".join(map(str, bits[::-1]))
    structure, active = decode_bitstring(bitstring, pairs, len(seq))
    energy = evaluate_structure_energy(seq, structure)
    runtime = time.perf_counter() - start

    return {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "sequence": seq,
        "length": len(seq),
        "method": method,
        "decoded_structure": structure,
        "active_pairs": str(active),
        "energy": energy,
        "qubo_cost": cost,
        "runtime_sec": round(runtime, 6),
        "notes": "exact" if method == "brute_force" else "heuristic",
    }


def main():
    OUT_CSV.parent.mkdir(exist_ok=True)

    fieldnames = [
        "timestamp",
        "sequence",
        "length",
        "method",
        "decoded_structure",
        "active_pairs",
        "energy",
        "qubo_cost",
        "runtime_sec",
        "notes",
    ]

    rows = []
    for seq in SEQUENCES:
        row = run_one(seq)
        rows.append(row)
        print(
            f"{seq} | {row['method']} | {row['decoded_structure']} | "
            f"energy={row['energy']} | cost={row['qubo_cost']} | runtime={row['runtime_sec']}s"
        )

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved results to {OUT_CSV}")


if __name__ == "__main__":
    main()