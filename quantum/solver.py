from itertools import product
from typing import Dict, List, Tuple

def brute_force_qubo_minimum(
    qubo: Dict[Tuple[int, int], float],
    num_vars: int
) -> Tuple[List[int], float]:
    best_bits = None
    best_cost = float("inf")

    for bits in product([0, 1], repeat=num_vars):
        cost = 0.0
        for (i, j), coeff in qubo.items():
            cost += coeff * bits[i] * bits[j]
        if cost < best_cost:
            best_cost = cost
            best_bits = list(bits)

    return best_bits, best_cost

if __name__ == "__main__":
    from quantum.qubo import build_qubo_matrix

    seq = "GCGCAUACGC"
    qubo, candidate_pairs = build_qubo_matrix(seq)
    best_bits, best_cost = brute_force_qubo_minimum(qubo, len(candidate_pairs))

    print("=== Brute Force Solver Test ===")
    print("Sequence       :", seq)
    print("Num vars       :", len(candidate_pairs))
    print("Best bitstring  :", best_bits)
    print("Best QUBO cost  :", best_cost)