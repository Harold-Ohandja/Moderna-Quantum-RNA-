"""
Quantum Resource and Scaling Analysis Module .
Evaluates qubit count, circuit depth, QUBO variables, and runtime scaling.
Identifies memory limits where local statevector simulation becomes impractical.
"""

import json
import os
import sys
import time
from typing import Dict, List
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# Internal module imports
from classical.utils import generate_random_rna
from quantum.qubo import build_qubo_matrix
from quantum.hamiltonian import qubo_to_hamiltonian


def run_scaling_analysis(
    sequence_lengths: List[int] = [6, 8, 10, 12, 14, 16, 18, 20, 22, 24],
    sim_threshold_qubits: int = 16
) -> List[Dict]:
    """
    Analyzes quantum resource growth across multiple RNA sequence lengths.
    """
    # Ensure output directories exist
    os.makedirs("data/results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    scaling_results = []

    print("=================================================================")
    print("       MILESTONE 6: QUANTUM RESOURCE & SCALING ANALYSIS         ")
    print("=================================================================")

    for length in sequence_lengths:
        sequence = generate_random_rna(length)
        
        start_time = time.time()
        
        # 1. Candidate Base Pairs & QUBO Matrix
        qubo_dict, candidate_pairs = build_qubo_matrix(sequence)
        num_variables = len(candidate_pairs)
        num_qubits = num_variables
        
        # 2. Convert to Ising Hamiltonian
        
        hamiltonian = qubo_to_hamiltonian(qubo_dict, num_qubits)
        
        build_time = time.time() - start_time
        
        # 3. Circuit Depth Estimation (n_local ansatz, RY+CZ, reps=2 -- see quantum/vqe_solver.py)
        # Depth formula for n_local(ry, cz, linear, reps=2): ~ 2 * 2 + 1 = 5 layers
        estimated_depth = 5 if num_qubits > 0 else 0

        # Determine Simulation Feasibility
        if num_qubits > sim_threshold_qubits:
            status = f"Simulation Impractical (> {sim_threshold_qubits} Qubits)"
        else:
            status = "Feasible for Local Simulation"

        record = {
            "sequence_length": length,
            "sequence": sequence,
            "num_candidate_pairs": num_variables,
            "num_qubits": num_qubits,
            "estimated_circuit_depth": estimated_depth,
            "qubo_terms_count": len(qubo_dict),
            "formulation_time_sec": round(build_time, 5),
            "status": status
        }
        
        scaling_results.append(record)
        
        print(f"Seq Length: {length:2d} nt | Qubits: {num_qubits:2d} | "
              f"QUBO Terms: {len(qubo_dict):3d} | Status: {status}")

    # Save metrics to JSON
    json_path = "data/results/scaling_data.json"
    with open(json_path, "w") as f:
        json.dump(scaling_results, f, indent=2)
        
    print("=================================================================")
    print(f"[OK] Scaling data saved to {json_path}")
    
    return scaling_results


def generate_scaling_plots(json_path: str = "data/results/scaling_data.json"):
    """
    Generates high-resolution figures for the report.
    """
    if not os.path.exists(json_path):
        print(f"Error: {json_path} not found. Run analysis first.")
        return

    with open(json_path, "r") as f:
        data = json.load(f)

    lengths = [d["sequence_length"] for d in data]
    qubits = [d["num_qubits"] for d in data]
    qubo_terms = [d["qubo_terms_count"] for d in data]

    # Plot 1: Sequence Length vs Qubit Count (# Variables)
    plt.figure(figsize=(8, 5))
    plt.plot(lengths, qubits, marker="o", color="#0A285F", linewidth=2.5, label="Qubits (Candidate Pairs)")
    plt.axhline(y=16, color="red", linestyle="--", linewidth=1.5, label="Local Sim Limit (16 Qubits)")
    plt.xlabel("RNA Sequence Length (nucleotides)", fontsize=11, fontweight="bold")
    plt.ylabel("Required Qubits / Binary Variables", fontsize=11, fontweight="bold")
    plt.title("Quantum Variable Scaling vs. RNA Length", fontsize=13, fontweight="bold", pad=12)
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(fontsize=10)
    plt.tight_layout()
    
    plot1_path = "figures/scaling_qubits_depth.png"
    plt.savefig(plot1_path, dpi=300)
    plt.close()

    # Plot 2: Sequence Length vs QUBO Complexity (Non-zero Terms)
    plt.figure(figsize=(8, 5))
    plt.plot(lengths, qubo_terms, marker="s", color="#008080", linewidth=2.5, label="QUBO Matrix Terms")
    plt.xlabel("RNA Sequence Length (nucleotides)", fontsize=11, fontweight="bold")
    plt.ylabel("Number of Non-Zero QUBO Terms", fontsize=11, fontweight="bold")
    plt.title("QUBO Matrix Interactions Complexity Scaling", fontsize=13, fontweight="bold", pad=12)
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(fontsize=10)
    plt.tight_layout()

    plot2_path = "figures/scaling_qubo_terms.png"
    plt.savefig(plot2_path, dpi=300)
    plt.close()

    print(f"[figures] Successfully saved to:\n   - {plot1_path}\n   - {plot2_path}")


if __name__ == "__main__":
    run_scaling_analysis()
    generate_scaling_plots()
