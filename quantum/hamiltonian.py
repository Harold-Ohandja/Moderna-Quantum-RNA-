"""
Hamiltonian Module for RNA QUBO Mapping.
Converts QUBO matrices into Qiskit SparsePauliOp Hamiltonians for VQE execution.
"""

from typing import Dict, Tuple
from qiskit.quantum_info import SparsePauliOp


def qubo_to_hamiltonian(qubo: Dict[Tuple[int, int], float], num_qubits: int) -> SparsePauliOp:
    """
    Converts a QUBO dictionary into a Qiskit SparsePauliOp Hamiltonian.
    
    Mapping from binary variable x_i in {0, 1} to Pauli Z operator:
    x_i = (I - Z_i) / 2

    Args:
        qubo (Dict[Tuple[int, int], float]): QUBO matrix with (i, j) tuples as keys.
        num_qubits (int): Total number of qubits required.

    Returns:
        SparsePauliOp: The Ising Hamiltonian operator in Qiskit format.
    """
    if num_qubits == 0:
        return SparsePauliOp.from_list([("I" * 1, 0.0)])

    # Initialize identity term operator
    pauli_list = []

    for (i, j), coeff in qubo.items():
        if i == j:
            # Linear term: c_i * x_i = c_i * (I - Z_i) / 2
            # Constant part: c_i / 2
            pauli_str_I = ["I"] * num_qubits
            pauli_list.append(("".join(pauli_str_I), coeff / 2.0))

            # Z_i part: -c_i / 2 * Z_i
            pauli_str_Z = ["I"] * num_qubits
            pauli_str_Z[num_qubits - 1 - i] = "Z"
            pauli_list.append(("".join(pauli_str_Z), -coeff / 2.0))

        else:
            # Quadratic term: c_ij * x_i * x_j = c_ij * (I - Z_i)(I - Z_j) / 4
            # = c_ij/4 * (I - Z_i - Z_j + Z_i Z_j)
            
            # Constant part
            pauli_str_I = ["I"] * num_qubits
            pauli_list.append(("".join(pauli_str_I), coeff / 4.0))

            # -Z_i part
            pauli_str_Zi = ["I"] * num_qubits
            pauli_str_Zi[num_qubits - 1 - i] = "Z"
            pauli_list.append(("".join(pauli_str_Zi), -coeff / 4.0))

            # -Z_j part
            pauli_str_Zj = ["I"] * num_qubits
            pauli_str_Zj[num_qubits - 1 - j] = "Z"
            pauli_list.append(("".join(pauli_str_Zj), -coeff / 4.0))

            # Z_i Z_j part
            pauli_str_Zij = ["I"] * num_qubits
            pauli_str_Zij[num_qubits - 1 - i] = "Z"
            pauli_str_Zij[num_qubits - 1 - j] = "Z"
            pauli_list.append(("".join(pauli_str_Zij), coeff / 4.0))

    # Build operator and simplify duplicate terms
    hamiltonian = SparsePauliOp.from_list(pauli_list).simplify()
    return hamiltonian


if __name__ == "__main__":
    # Test on a small 2-qubit QUBO matrix
    test_qubo = {(0, 0): -1.0, (1, 1): -1.0, (0, 1): 2.0}
    num_q = 2

    H = qubo_to_hamiltonian(test_qubo, num_q)
    print("=== HAMILTONIAN CONVERSION TEST ===")
    print("QUBO Matrix:", test_qubo)
    print("Mapped Hamiltonian:\n", H)
