"""
Greedy Local-Search Heuristic QUBO Solver.

NOTE: this is a classical baseline, not a quantum or quantum-inspired
algorithm (see the same note in quantum/solver.py). It's a multi-restart
1-opt bit-flip hill climber, used when a sequence is too large for the
exact brute-force solver in solver.py. Useful as a fast classical
comparison point for larger sequences; not part of the quantum contribution.
"""


def qubo_cost(qubo, bits):
    cost = 0.0
    for (i, j), coeff in qubo.items():
        cost += coeff * bits[i] * bits[j]
    return cost


def greedy_local_search_qubo(qubo, n_vars, max_iters=1000, restarts=10):
    best_bits = None
    best_cost = float('inf')

    for r in range(restarts):
        bits = [0] * n_vars
        for i in range(n_vars):
            if (i + r) % 2 == 0:
                bits[i] = 1 if (i + r) % 3 == 0 else 0

        current_cost = qubo_cost(qubo, bits)

        for _ in range(max_iters):
            candidate_best_cost = current_cost
            candidate_flip = None

            for i in range(n_vars):
                bits[i] ^= 1
                cost = qubo_cost(qubo, bits)
                bits[i] ^= 1
                if cost < candidate_best_cost:
                    candidate_best_cost = cost
                    candidate_flip = i

            if candidate_flip is None:
                break

            bits[candidate_flip] ^= 1
            current_cost = candidate_best_cost

        if current_cost < best_cost:
            best_cost = current_cost
            best_bits = bits[:]

    return best_bits, best_cost