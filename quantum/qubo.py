"""
QUBO Formulation Module for RNA Folding.
Implements base-pair candidate mapping, constraint penalties,
and pseudoknot exclusion (Alevras et al. model, arXiv:2405.20328).
"""

from typing import Dict, List, Tuple

# Allowed canonical base pairs (Watson-Crick + GU Wobble)
CANONICAL_PAIRS = {("A", "U"), ("U", "A"), ("C", "G"), ("G", "C"), ("G", "U"), ("U", "G")}


def get_possible_base_pairs(sequence: str, min_loop_length: int = 3) -> List[Tuple[int, int]]:
    """
    Identifies all physically feasible candidate base pairs.
    
    Args:
        sequence (str): RNA sequence (A, C, G, U).
        min_loop_length (int): Minimum hairpin loop length (default = 3).
        
    Returns:
        List[Tuple[int, int]]: List of valid base pair index tuples (i, j).
    """
    seq = sequence.upper().replace("T", "U")
    n = len(seq)
    pairs = []

    for i in range(n):
        for j in range(i + min_loop_length + 1, n):
            if (seq[i], seq[j]) in CANONICAL_PAIRS:
                pairs.append((i, j))

    return pairs


def build_qubo_matrix(
    sequence: str,
    penalty_weight: float = 5.0,
    stacking_bonus: float = -1.5,
    pair_energy: float = -1.0,
) -> Tuple[Dict[Tuple[int, int], float], List[Tuple[int, int]]]:
    """
    Constructs the QUBO dictionary for the given RNA sequence.
    
    Args:
        sequence (str): RNA sequence.
        penalty_weight (float): Penalty (P) applied to conflicting/crossing pairs.
        stacking_bonus (float): Free energy bonus for adjacent base pairs (stacking).
        pair_energy (float): Base free energy contribution for a single pair formation.
        
    Returns:
        Tuple[Dict, List]: QUBO dictionary {(i, j): value} and list of candidate pairs.
    """
    candidate_pairs = get_possible_base_pairs(sequence)
    num_vars = len(candidate_pairs)
    qubo = {}

    # 1. Diagonal Terms (Individual pair energy)
    for i in range(num_vars):
        qubo[(i, i)] = pair_energy

    # 2. Stacking Terms (Energy bonus for contiguous base pairs)
    for i in range(num_vars):
        p1 = candidate_pairs[i]
        for j in range(i + 1, num_vars):
            p2 = candidate_pairs[j]
            # Check if p1 (i1, j1) and p2 (i2, j2) form a stacked pair ((...))
            if p2[0] == p1[0] + 1 and p2[1] == p1[1] - 1:
                qubo[(i, j)] = stacking_bonus

    # 3. Constraints: Shared-base Incompatibility & Non-crossing (Pseudoknots Excluded)
    # If two base pairs share a nucleotide OR cross each other, apply penalty_weight.
    for i in range(num_vars):
        i1, j1 = candidate_pairs[i]
        for j in range(i + 1, num_vars):
            i2, j2 = candidate_pairs[j]

            # Constraint 1: A single base can pair at most once
            shares_base = (i1 == i2) or (i1 == j2) or (j1 == i2) or (j1 == j2)

            # Constraint 2: Non-crossing condition (Excluding pseudoknots: i1 < i2 < j1 < j2)
            crosses = (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)

            if shares_base or crosses:
                # Add quadratic penalty P * x_i * x_j
                qubo[(i, j)] = qubo.get((i, j), 0.0) + penalty_weight

    return qubo, candidate_pairs


if __name__ == "__main__":
    # Quick test on the Toy Sequence (10 nt)
    toy_seq = "GCGCAUACGC"
    qubo_dict, pairs = build_qubo_matrix(toy_seq)

    print("=== QUBO FORMULATION TEST ===")
    print(f"Sequence            : {toy_seq}")
    print(f"Number of pairs     : {len(pairs)}")
    print(f"Candidate pairs     : {pairs}")
    print(f"QUBO matrix size    : {len(qubo_dict)} non-zero terms")
    print("\nQUBO Matrix Sample:")
    for k, v in list(qubo_dict.items())[:5]:
        print(f"  x_{k[0]} * x_{k[1]} : {v}")
