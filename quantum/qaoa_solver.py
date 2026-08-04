"""
QAOA Solver Module for RNA Folding.
Secondary / comparison quantum solver, used for the "compare multiple
quantum encodings and approaches" optional advanced task. Uses standard
expectation-value QAOA (qiskit_algorithms.QAOA) so its behavior can be
directly contrasted with the CVaR-VQE solver in vqe_solver.py.
"""

from typing import Callable, Dict, List, Optional, Tuple

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np
from qiskit_aer import AerSimulator
from qiskit_aer.primitives import SamplerV2
from qiskit_algorithms import QAOA
from qiskit_algorithms.optimizers import COBYLA, Optimizer
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from quantum.hamiltonian import qubo_to_hamiltonian
from quantum.decode import decode_bitstring


def solve_qaoa(
    qubo: Dict[Tuple[int, int], float],
    candidate_pairs: List[Tuple[int, int]],
    seq_len: int,
    reps: int = 2,
    maxiter: int = 150,
    aggregation: Optional[float] = None,
    optimizer: Optional[Optimizer] = None,
    seed: Optional[int] = None,
    callback: Optional[Callable] = None,
) -> Dict:
    """
    Solves the RNA-folding QUBO with QAOA.

    Args:
        qubo: QUBO dictionary from quantum.qubo.build_qubo_matrix.
        candidate_pairs: Candidate base pairs, same order as the QUBO.
        seq_len: Length of the RNA sequence.
        reps: QAOA depth parameter p.
        maxiter: Maximum COBYLA iterations (ignored if a custom optimizer
            is supplied).
        aggregation: If None, uses the standard full expectation value
            (plain QAOA). If a float in (0, 1], uses a CVaR objective at
            that alpha instead, for a CVaR-QAOA variant.
        optimizer: Optional custom qiskit_algorithms Optimizer. Defaults to
            COBYLA(maxiter=maxiter).
        seed: Optional seed for the Aer sampler, for reproducibility.
        callback: Optional callback(eval_count, params, value, metadata).

    Returns:
        Dict with the raw result, decoded structure, active pairs, best
        bitstring, and the QUBO-space eigenvalue.

    Note:
        QAOA builds its ansatz (QAOAAnsatz) internally, which Aer's
        SamplerV2 cannot execute directly unless it is decomposed first.
        We pass a preset pass manager as the `transpiler` so this happens
        automatically; without it, Aer raises
        `AerError: unknown instruction: QAOA`.
    """
    num_qubits = len(candidate_pairs)
    hamiltonian = qubo_to_hamiltonian(qubo, num_qubits)

    backend = AerSimulator()
    pass_manager = generate_preset_pass_manager(optimization_level=1, backend=backend)

    sampler = SamplerV2(seed=seed) if seed is not None else SamplerV2()
    opt = optimizer if optimizer is not None else COBYLA(maxiter=maxiter)

    # See the matching note in vqe_solver.py: seed a local RNG for the initial
    # point explicitly, rather than relying on SamplerV2's seed alone, so a
    # given seed reproduces the same result regardless of call history.
    # QAOAAnsatz has 2 parameters per rep (one gamma, one beta).
    if seed is not None:
        rng = np.random.default_rng(seed)
        initial_point = rng.uniform(-2 * np.pi, 2 * np.pi, size=2 * reps)
    else:
        initial_point = None

    qaoa = QAOA(
        sampler=sampler,
        optimizer=opt,
        reps=reps,
        aggregation=aggregation,
        initial_point=initial_point,
        transpiler=pass_manager,
        callback=callback,
    )
    result = qaoa.compute_minimum_eigenvalue(hamiltonian)

    bitstring = result.best_measurement["bitstring"]
    structure, active_pairs = decode_bitstring(bitstring, candidate_pairs, seq_len)

    return {
        "method": "CVaR-QAOA" if aggregation else "QAOA",
        "aggregation": aggregation,
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

    parser = argparse.ArgumentParser(description="Run the QAOA RNA folding solver")
    parser.add_argument("--sequence", type=str, default="GCGCAUACGC", help="RNA sequence")
    parser.add_argument("--reps", type=int, default=2, help="QAOA depth (p)")
    parser.add_argument("--maxiter", type=int, default=150, help="Max COBYLA iterations")
    parser.add_argument(
        "--cvar-alpha", type=float, default=None,
        help="If set, uses CVaR-QAOA at this alpha instead of plain QAOA",
    )
    parser.add_argument("--seed", type=int, default=42, help="Seed, for reproducibility")
    args = parser.parse_args()

    seq = args.sequence.strip().upper()
    qubo, pairs = build_qubo_matrix(seq)

    out = solve_qaoa(
        qubo, pairs, len(seq),
        reps=args.reps, maxiter=args.maxiter, aggregation=args.cvar_alpha, seed=args.seed,
    )
    print(f"{out['method']} structure:", out["structure"])

    ref_struct, ref_energy = RNA.fold(seq)
    fc = RNA.fold_compound(seq)
    predicted_energy = fc.eval_structure(out["structure"])

    print("Reference MFE   :", ref_struct, ref_energy)
    print("Predicted energy:", predicted_energy)
    print("Match:", out["structure"] == ref_struct)
