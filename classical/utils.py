"""
Utility functions for sequence generation and benchmark helpers.
"""

import random

def generate_random_rna(length: int, seed: int = 42) -> str:
    """
    Generates a deterministic pseudo-random RNA sequence of a given length.

    Args:
        length (int): Desired nucleotide sequence length.
        seed (int): Random seed for reproducibility.

    Returns:
        str: RNA sequence string composed of A, C, G, U.
    """
    random.seed(seed + length)
    bases = ["A", "C", "G", "U"]
    return "".join(random.choice(bases) for _ in range(length))


if __name__ == "__main__":
    print("Test sequence (10 nt):", generate_random_rna(10))
