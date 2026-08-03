"""
Bitstring Decoding Module for RNA QUBO Solutions.
Converts a measured bitstring (from CVaR-VQE or QAOA) back into a
dot-bracket secondary structure, using the same candidate-pair list
produced by quantum.qubo.build_qubo_matrix.
"""

from typing import List, Tuple


def decode_bitstring(bitstring: str, candidate_pairs: List[Tuple[int, int]], seq_len: int) -> Tuple[str, List[Tuple[int, int]]]:
    """
    Decodes a measured bitstring into an RNA dot-bracket structure.

    Bit ordering follows Qiskit's convention (same one used when building the
    Hamiltonian in quantum.hamiltonian): the rightmost character of the
    bitstring corresponds to variable/qubit 0, i.e. variable i's value is
    bitstring[num_qubits - 1 - i].

    Args:
        bitstring (str): Measured bitstring, e.g. "0011100".
        candidate_pairs (List[Tuple[int, int]]): Candidate base pairs, in the
            same order returned by build_qubo_matrix.
        seq_len (int): Length of the RNA sequence.

    Returns:
        Tuple[str, List[Tuple[int, int]]]: Dot-bracket structure string and
        the list of active (selected) base pairs.
    """
    n = len(candidate_pairs)
    if len(bitstring) != n:
        raise ValueError(
            f"Bitstring length ({len(bitstring)}) does not match "
            f"number of candidate pairs ({n})."
        )

    active_pairs = []
    for i, pair in enumerate(candidate_pairs):
        bit = bitstring[n - 1 - i]
        if bit == "1":
            active_pairs.append(pair)

    structure = ["."] * seq_len
    for (i, j) in active_pairs:
        structure[i] = "("
        structure[j] = ")"

    return "".join(structure), active_pairs


if __name__ == "__main__":
    # Sanity check against the known toy-sequence result
    pairs = [(0, 5), (0, 7), (0, 9), (1, 8), (2, 7), (2, 9), (3, 8)]
    struct, active = decode_bitstring("0011100", pairs, 10)
    print("Structure:", struct)
    print("Active pairs:", active)
    assert struct == "(((....)))", "Decode mismatch against known toy result"
    print("OK: matches known ViennaRNA MFE structure for GCGCAUACGC")
