# Milestone 6 Report — Quantum Solvers Implementation, Resource Scaling & Comparative Benchmarking

## Overview
Milestone 6 marked the transition from classical baseline validation to full-scale quantum simulation. The primary objective was to integrate, execute, and validate the Variational Quantum Eigensolver (VQE) and the Quantum Alternating Operator Ansatz (QAOA) modules (`quantum/vqe_solver.py` and `quantum/qaoa_solver.py`) using Qiskit. 

Additionally, a comprehensive resource scaling analysis was performed via `quantum/resource_analysis.py` to evaluate qubit counts, circuit depth, QUBO variables, interaction complexity, and runtime constraints, identifying the exact boundary where local classical simulation becomes impractical.

---

## Objectives
The key milestones achieved during Milestone 6 were:
- **VQE & CVaR-VQE Integration:** Implement and execute VQE using a hardware-efficient `n_local` ansatz and the COBYLA classical optimizer.
- **QAOA Solver Implementation:** Deploy QAOA with alternating problem and mixer Hamiltonians to evaluate quantum approximate optimization performance.
- **Quantum-to-Classical Alignment:** Validate that the quantum ground state bitstrings map to valid RNA secondary structures via `quantum/decode.py`.
- **Resource & Scaling Benchmark:** Evaluate qubit count, circuit depth, QUBO matrix variables, and runtimes across sequence lengths ($L \in [6, 24]\text{ nt}$).
- **Simulation Limit Logging:** Log the hardware memory wall where local statevector simulation becomes infeasible ($N > 16$ qubits).

---

## Quantum Solvers Architecture

### 1. Variational Quantum Eigensolver (VQE)
The VQE module (`quantum/vqe_solver.py`) computes the ground state of the RNA Ising Hamiltonian $H_P$ by minimizing the expectation value:

$$E(\theta) = \langle \psi(\theta) | H_P | \psi(\theta) \rangle$$

- **Ansatz:** `qiskit.circuit.library.n_local` (RY rotations + linear CZ entanglement). Note: the older `TwoLocal` class was deliberately avoided, it is deprecated as of Qiskit 2.1 and, more importantly, its wrapped-instruction form fails outright under Aer's `SamplerV2` (`AerError: unknown instruction: TwoLocal`) unless manually decomposed first. `n_local` returns a plain `QuantumCircuit` and has no such issue.
- **Optimizer:** COBYLA.
- **CVaR Extension:** Conditional Value-at-Risk objective (alpha=0.1 by default), which optimizes the mean energy of only the best-sampled bitstrings each iteration rather than the full expectation value, converging faster on this kind of rugged combinatorial landscape. Following Alevras et al., *mRNA secondary structure prediction using utility-scale quantum computers* (arXiv:2405.20328), the Moderna + IBM Quantum paper this challenge is based on.

### 2. Quantum Alternating Operator Ansatz (QAOA)
The QAOA module (`quantum/qaoa_solver.py`) applies $p$ alternating layers of problem ($H_P$) and mixer ($H_M = \sum X_i$) Hamiltonians:

$$|\gamma, \beta\rangle = \prod_{k=1}^p e^{-i \beta_k H_M} e^{-i \gamma_k H_P} |+\rangle^{\otimes N}$$

---

## 📊 Quantum Resource Scaling & Hardware Limits Analysis

Using `quantum/resource_analysis.py`, we evaluated how computational requirements scale with RNA sequence length.

### 1. Qubit Requirements Scaling vs. Local Simulation Limits

![Quantum Variable Scaling vs RNA Length](figures/scaling_qubits_depth.png)

*Figure 1: Scaling of required qubits (candidate base pairs) relative to RNA sequence length, highlighting the 16-qubit local statevector simulation threshold.*

#### Observations & Logs:
* **Variable Scaling:** The number of qubits corresponds to the candidate base pairs, scaling from **2 qubits (6 nt)** up to **78 qubits (24 nt)**.
* **The 16-Qubit Wall:** Local statevector simulation requires storing $2^N$ complex amplitudes. Above 16 qubits, memory consumption exceeds standard workstation limits.
* **Terminal Status Logs:**
  * **16 nt:** 30 qubits | `Status: Simulation Impractical (> 16 Qubits)`
  * **18 nt:** 37 qubits | `Status: Simulation Impractical (> 16 Qubits)`
  * **22 nt:** 48 qubits | `Status: Simulation Impractical (> 16 Qubits)`
  * **24 nt:** 78 qubits | `Status: Simulation Impractical (> 16 Qubits)`

---

### 2. QUBO Matrix Interaction Complexity Scaling

![QUBO Matrix Interactions Complexity Scaling](figures/scaling_qubo_terms.png)

*Figure 2: Scaling of non-zero interactions in the QUBO matrix as a function of RNA sequence length.*

#### Observations & Logs:
* **Interaction Density:** Penalty constraints (overlapping base pairs and pseudoknots) cause quadratic growth in non-zero QUBO matrix terms.
* **Term Count:** Non-zero terms grow from **3 (6 nt)** and **26 (10 nt)** to **1,934 non-zero terms (24 nt)**.
* **Circuit Impact:** Each non-zero term requires a 2-qubit $ZZ$ interaction layer. For 24 nt, 1,934 $ZZ$ gates create significant circuit depth and CNOT gate overhead.

---

### 3. Scaling Data Summary Table

| RNA Length (nt) | Required Qubits ($N$) | Non-Zero QUBO Terms | Formulation Time (s) | Local Simulation Status |
| :---: | :---: | :---: | :---: | :--- |
| **6** | 2 | 3 | < 0.001 | ✅ Feasible |
| **8** | 4 | 9 | < 0.001 | ✅ Feasible |
| **10** | 9 | 26 | ~ 0.002 | ✅ Feasible |
| **12** | 9 | 26 | ~ 0.002 | ✅ Feasible |
| **14** | 19 | 148 | ~ 0.005 | ❌ Impractical ($>16$ Qubits) |
| **16** | 30 | 382 | ~ 0.012 | ❌ Impractical ($>16$ Qubits) |
| **18** | 37 | 475 | ~ 0.018 | ❌ Impractical ($>16$ Qubits) |
| **20** | 36 | 464 | ~ 0.016 | ❌ Impractical ($>16$ Qubits) |
| **22** | 48 | 768 | ~ 0.035 | ❌ Impractical ($>16$ Qubits) |
| **24** | 78 | 1,934 | ~ 0.080 | ❌ Impractical ($>16$ Qubits) |

---

## Key Experimental Results & Comparative Benchmarks

*(All figures below were generated by actually running each solver against the true ViennaRNA MFE reference, not estimated.)*

| Sequence | Solver Method | Qubits ($N$) | Bitstring | Decoded Structure | Predicted Energy (kcal/mol) | True ViennaRNA MFE (kcal/mol) | QUBO Cost | Runtime (s) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `GCGCAUACGC` | Exact (Brute-Force) | 7 | `0011100` | `(((....)))` | -1.30 | -1.30 (`(((....)))`) | -6.0 | 0.0002 |
| `GCGCAUACGC` | CVaR-VQE (COBYLA) | 7 | `0011100` | `(((....)))` | -1.30 | -1.30 (`(((....)))`) | -6.0 | 0.79 |
| `GCGCAUACGC` | QAOA (p=2) | 7 | `0011100` | `(((....)))` | -1.30 | -1.30 (`(((....)))`) | 10.9 | 3.56 |
| `GCAUCGUAGC` | Exact (Brute-Force) | 10 | `0001010100` | `((.(...)))` | 4.20 | 0.00 (`..........`) | -4.5 | 0.003 |
| `GCAUCGUAGC` | CVaR-VQE (COBYLA) | 10 | `0001010100` | `((.(...)))` | 4.20 | 0.00 (`..........`) | -4.5 | 1.47 |
| `GCAUCGUAGC` | QAOA (p=2) | 10 | `0001010100` | `((.(...)))` | 4.20 | 0.00 (`..........`) | 17.7 | 0.81 |
| `GCGCGCGAAAUUCGCGCG` | Heuristic (1-opt) | 37 | -- | `.((((((....).))))` | -4.30 | -10.70 (`.(((((((...)))))))`) | -12.0 | 0.13 |

**Reading this honestly:** on `GCGCAUACGC`, every method (exact classical, CVaR-VQE, QAOA) independently reaches the exact same correct answer. On `GCAUCGUAGC`, every method (again including exact brute-force) converges on the *same wrong* structure, which is strong evidence this is a QUBO formulation gap rather than any single solver failing, since the exact solver can't do better than the objective it was given. On the longer sequence, the heuristic solver (used because 37 qubits is outside this project's local-simulation budget) finds a structure that's directionally reasonable but far from the true MFE, both because of the same formulation gap and because 1-opt local search isn't guaranteed to find the QUBO's own global optimum the way brute-force does. See `notebooks/04_Full_Comparison.ipynb` for the fully reproducible version of this comparison, including additional sequences.

---

## Code Additions & Modifications (Milestone 6)
- **`quantum/vqe_solver.py`:** Modular CVaR-VQE implementation using an `n_local` ansatz and COBYLA.
- **`quantum/qaoa_solver.py`:** QAOA solver supporting parameterized problem/mixer layers.
- **`quantum/resource_analysis.py`:** Resource profiling script for qubit count, circuit depth, QUBO term density, and simulation feasibility logging.
- **`data/results/scaling_data.json`:** Raw JSON dataset containing resource scaling logs.

---

## Milestone 6 Summary
Milestone 6 successfully demonstrated end-to-end quantum execution for RNA secondary structure prediction. Both VQE and QAOA matched the exact classical QUBO ground state on simulator backends. Furthermore, the resource analysis clearly mapped the exponential growth of quantum variables and interaction terms, identifying $N = 16$ qubits as the hard threshold where local statevector simulation becomes impractical.
