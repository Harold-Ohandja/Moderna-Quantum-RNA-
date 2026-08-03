import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from quantum.qubo import build_qubo_matrix
from quantum.decode import decode_bitstring
from classical.evaluate_energy import evaluate_structure_energy
from quantum.hybrid_solver import solve_qubo_auto


def main():
    parser = argparse.ArgumentParser(description="Run RNA QUBO solver")
    parser.add_argument("--seq", required=True, help="RNA sequence")
    parser.add_argument("--brute-force-limit", type=int, default=12)
    args = parser.parse_args()

    seq = args.seq.strip().upper()
    qubo, pairs = build_qubo_matrix(seq)
    bits, cost, method = solve_qubo_auto(qubo, len(pairs), brute_force_limit=args.brute_force_limit)

    bitstring = ''.join(map(str, bits[::-1]))
    structure, active = decode_bitstring(bitstring, pairs, len(seq))
    energy = evaluate_structure_energy(seq, structure)

    print(f"method: {method}")
    print(f"sequence: {seq}")
    print(f"bits: {bits}")
    print(f"decoded: {structure}")
    print(f"active: {active}")
    print(f"energy: {energy}")
    print(f"qubo_cost: {cost}")


if __name__ == "__main__":
    main()