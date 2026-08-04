"""
Auto-Switching Classical QUBO Solver.

NOTE: this dispatches between the two classical baselines above
(exact brute-force for small instances, greedy heuristic for larger ones).
Neither is quantum or quantum-inspired; see quantum/solver.py for details.
"""

from quantum.solver import brute_force_qubo_minimum
from quantum.heuristic_solver import greedy_local_search_qubo


def solve_qubo_auto(qubo, n_vars, brute_force_limit=12):
    if n_vars <= brute_force_limit:
        method = 'brute_force'
        bits, cost = brute_force_qubo_minimum(qubo, n_vars)
    else:
        method = 'heuristic'
        bits, cost = greedy_local_search_qubo(qubo, n_vars)
    return bits, cost, method