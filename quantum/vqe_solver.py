"""
CVaR-VQE Solver Module for RNA Folding.
Primary quantum solver for the project, following the approach used by
Alevras et al. (arXiv:2405.20328): a hardware-efficient ansatz optimized
with a Conditional Value at Risk (CVaR) objective, which converges faster
than plain expectation-value VQE on rugged combinatorial landscapes.
"""

from typing import Callable, Dict, List, Optional, Tuple

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import n_local
from qiskit_algorithms import SamplingVQE
from qiskit_algorithms.optimizers import COBYLA, Optimizer
from qiskit_aer.primitives import SamplerV2

from quantum.hamiltonian import qubo_to_hamiltonian
from quantum.decode import decode_bitstring


def build_ansatz(num_qubits: int, reps: int = 2) -> QuantumCircuit:
    """
    Builds a hardware-efficient Two-Local-style ansatz (RY rotations + CZ
    entanglement, linear connectivity), matching the ansatz family used in
    the reference paper.

    Uses qiskit.circuit.library.n_local rather than the deprecated TwoLocal
    class (TwoLocal is deprecated as of Qiskit 2.1 and will be removed in
    Qiskit 3.0). n_local returns a plain QuantumCircuit, which Aer's
    SamplerV2 can execute directly without needing manual decomposition.
    """
    return n_local(
        num_qubits,
        rotation_blocks="ry",
        entanglement_blocks="cz",
        entanglement="linear",
        reps=reps,
    )


def solve_cvar_vqe(
    qubo: Dict[Tuple[int, int], float],
    candidate_pairs: List[Tuple[int, int]],
    seq_len: int,
    alpha: float = 0.1,
    reps: int = 2,
    maxiter: int = 150,
    optimizer: Optional[Optimizer] = None,
    seed: Optional[int] = None,
    callback: Optional[Callable] = None,
) -> Dict:
    """
    Solves the RNA-folding QUBO with CVaR-VQE.

    Args:
        qubo: QUBO dictionary from quantum.qubo.build_qubo_matrix.
        candidate_pairs: Candidate base pairs, same order as the QUBO.
        seq_len: Length of the RNA sequence.
        alpha: CVaR quantile (0 < alpha <= 1). alpha=0.1 optimizes the mean
            of the best 10% of sampled bitstrings each iteration, which is
            the setting used in the reference paper. alpha=1.0 recovers
            standard expectation-value VQE.
        reps: Ansatz repetition depth.
        maxiter: Maximum COBYLA iterations (ignored if a custom optimizer
            is supplied).
        optimizer: Optional custom qiskit_algorithms Optimizer. Defaults to
            COBYLA(maxiter=maxiter).
        seed: Optional seed for the Aer sampler, for reproducibility.
        callback: Optional callback(eval_count, params, value, metadata)
            forwarded to SamplingVQE, useful for plotting convergence.

    Returns:
        Dict with the raw result, decoded structure, active pairs, best
        bitstring, and the QUBO-space eigenvalue.
    """
    num_qubits = len(candidate_pairs)
    hamiltonian = qubo_to_hamiltonian(qubo, num_qubits)

    ansatz = build_ansatz(num_qubits, reps=reps)
    sampler = SamplerV2(seed=seed) if seed is not None else SamplerV2()
    opt = optimizer if optimizer is not None else COBYLA(maxiter=maxiter)

    # NOTE on reproducibility: SamplerV2's `seed` only controls measurement-shot
    # randomness. If no initial_point is given, SamplingVQE draws one from
    # numpy's *global* RNG state, which depends on how many random calls
    # happened earlier in the process -- so the same `seed` can silently give
    # different results depending on what ran before it. We fix this by
    # explicitly seeding a local RNG and drawing the initial point ourselves,
    # so a given seed reproduces the same result regardless of call history.
    if seed is not None:
        rng = np.random.default_rng(seed)
        initial_point = rng.uniform(-2 * np.pi, 2 * np.pi, size=ansatz.num_parameters)
    else:
        initial_point = None

    vqe = SamplingVQE(
        sampler=sampler,
        ansatz=ansatz,
        optimizer=opt,
        aggregation=alpha,
        initial_point=initial_point,
        callback=callback,
    )
    result = vqe.compute_minimum_eigenvalue(hamiltonian)

    bitstring = result.best_measurement["bitstring"]
    structure, active_pairs = decode_bitstring(bitstring, candidate_pairs, seq_len)

    return {
        "method": "CVaR-VQE",
        "alpha": alpha,
        "qubo_eigenvalue": result.eigenvalue,
        "bitstring": bitstring,
        "structure": structure,
        "active_pairs": active_pairs,
        "raw_result": result,
    }


if __name__ == "__main__":
    import argparse
    from quantum.qubo import build_qubo_matrix
    import RNA

    parser = argparse.ArgumentParser(description="Run the CVaR-VQE RNA folding solver")
    parser.add_argument("--sequence", type=str, default="GCGCAUACGC", help="RNA sequence")
    parser.add_argument("--alpha", type=float, default=0.1, help="CVaR quantile (0 < alpha <= 1)")
    parser.add_argument("--reps", type=int, default=2, help="Ansatz repetition depth")
    parser.add_argument("--maxiter", type=int, default=150, help="Max COBYLA iterations")
    parser.add_argument("--seed", type=int, default=42, help="Seed, for reproducibility")
    args = parser.parse_args()

    seq = args.sequence.strip().upper()
    qubo, pairs = build_qubo_matrix(seq)

    out = solve_cvar_vqe(
        qubo, pairs, len(seq),
        alpha=args.alpha, reps=args.reps, maxiter=args.maxiter, seed=args.seed,
    )
    print("CVaR-VQE structure:", out["structure"])

    ref_struct, ref_energy = RNA.fold(seq)
    fc = RNA.fold_compound(seq)
    predicted_energy = fc.eval_structure(out["structure"])

    print("Reference MFE     :", ref_struct, ref_energy)
    print("Predicted energy  :", predicted_energy)
    print("Match:", out["structure"] == ref_struct)
