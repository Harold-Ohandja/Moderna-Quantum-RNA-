
# Moderna Challenge Report: Day 1–2 Progress

**Project:** Optimization of mRNA Secondary Structure Prediction Using Quantum Computing  
**Challenge:** WISER Global Quantum+AI Program 2026 — Moderna Challenge  
**Team:** Pushkar Kumar and Harold   
**Status:** Classical baseline and QUBO formulation validated on toy sequences

## 1. Objective

This project targets the Moderna WISER challenge of predicting RNA secondary structure with a quantum or quantum-inspired optimization method and benchmarking results against ViennaRNA minimum free energy (MFE) folding. The goal is not to outperform classical methods, but to reproduce known structures on short RNA sequences and analyze how the computational cost scales with sequence length.

The work completed so far covers the first two stages of the plan:
- Day 1: classical benchmarking with ViennaRNA.
- Day 2: QUBO formulation with candidate base-pair variables, stacking rewards, and structural constraints.

## 2. Challenge Framing

The challenge asks participants to:
- Formulate RNA secondary structure prediction as an optimization problem.
- Use a quantum or quantum-inspired method to identify low-energy candidate structures.
- Benchmark against ViennaRNA MFE structures.
- Analyze scaling in terms of qubit count, circuit depth, number of variables, and runtime.

The internal team plan identified the IBM–Moderna reference work as the main technical anchor:
- Alevras et al., *mRNA secondary structure prediction using utility-scale quantum computers* (arXiv:2405.20328).
- The selected formulation is a QUBO model with binary variables representing candidate base pairs.
- The planned solver approach is CVaR-VQE, though solver work is not yet implemented.

## 3. Repository Status

The repository is still being built toward the ideal structure described in the README. At the moment, the following practical status applies:

### Present and working
- `classical/generate_mfe.py`
- `classical/evaluate_energy.py`
- `quantum/qubo.py`

### Intended but not yet implemented
- `classical/benchmark.py`
- `classical/utils.py`
- `quantum/hamiltonian.py`
- `quantum/qaoa_solver.py`
- `quantum/vqe_solver.py`
- `quantum/resource_analysis.py`
- `notebooks/`
- `figures/`
- `report/`
- `presentation/`

This means the README is currently more ambitious than the codebase. That is fine at this stage, as long as the implementation continues to close the gap step by step.

## 4. Day 1 Work: Classical Benchmarking

The first milestone was to verify that ViennaRNA works correctly in the local environment and can generate MFE structures and energies.

### 4.1 Environment check
The local Python virtual environment was confirmed to be working. The required libraries were installed successfully, and the project modules imported without error.

### 4.2 Toy sequence benchmark
The 10-nt toy sequence used for the initial classical test was:

- Sequence: `GCGCAUACGC`

ViennaRNA generated the following MFE result:

- Structure: `(((....)))`
- Energy: `-1.30 kcal/mol`

This result was also independently checked using the custom energy-evaluation module, which returned `-1.3` kcal/mol for the same sequence and structure. That confirms consistency between the ViennaRNA folding path and the energy evaluation path.

### 4.3 Moderna official example sequence
The 44-nt example sequence from the challenge materials was also tested:

- Sequence: `GGAGCAAAACUUGUCGAUUGAGAACAAAAUACAGAAUUUGCUUG`

The current script prints the following reference output:

- Structure: `.(((((((..((((...(((....)))...))))..))))))).`
- Energy: `-7.90 kcal/mol`

This gives the project a direct classical reference point for the challenge’s official example sequence.

### 4.4 Day 1 conclusion
Day 1 is effectively complete. ViennaRNA is functioning correctly, classical MFE references are available, and the energy evaluator agrees with the toy benchmark output.

## 5. Day 2 Work: QUBO Formulation

The second milestone was to define the optimization model for RNA secondary structure as a QUBO problem.

### 5.1 Candidate base-pair variables
The `quantum/qubo.py` module identifies candidate base pairs that satisfy:
- Watson–Crick pairing: A–U, U–A, C–G, G–C
- Wobble pairing: G–U, U–G
- Minimum loop-length constraint

Each candidate pair becomes a binary variable:
- `x_k = 1` if candidate pair `k` forms.
- `x_k = 0` otherwise.

This is the core variable mapping needed for a QUBO or Ising-style optimizer.

### 5.2 Energy terms
The implemented QUBO includes:
- A negative diagonal term for each candidate pair, encouraging pair formation.
- A stacking bonus for adjacent nested pairs.
- Penalty terms for invalid combinations.

The current simplified energy model is intentionally lightweight. It is sufficient for proof-of-concept testing on toy sequences and can be refined later toward a closer nearest-neighbor thermodynamic model.

### 5.3 Structural constraints
Two key constraints are already encoded:
- A nucleotide can participate in at most one base pair.
- Crossing pairs are penalized, which excludes pseudoknots.

This is aligned with the team plan, where pseudoknot exclusion was explicitly listed as a valid simplification for the initial version.

### 5.4 Toy QUBO test: `GCGCAUACGC`
For the 10-nt toy sequence:

- Number of candidate pairs: `7`
- Candidate pairs: `[(0, 5), (0, 7), (0, 9), (1, 8), (2, 7), (2, 9), (3, 8)]`
- QUBO non-zero terms: `26`

This shows that the formulation is non-trivial even on a small sequence and that the QUBO matrix is being populated correctly.

### 5.5 Second QUBO test: `AUGCAUGC`
For the second short sequence:

- Sequence: `AUGCAUGC`
- Number of candidate pairs: `3`
- Candidate pairs: `[(0, 5), (1, 6), (2, 7)]`
- QUBO non-zero terms: `6`

This confirms that the QUBO formulation responds naturally to sequence length and composition.

### 5.6 Day 2 conclusion
Day 2 is largely complete at the formulation level. The QUBO construction exists, works on more than one toy sequence, and captures the intended structural constraints.

## 6. Classical Energy Evaluation

The `classical/evaluate_energy.py` module provides two useful functions:
- `evaluate_structure_energy(sequence, structure)`
- `calculate_energy_gap(predicted_energy, reference_mfe_energy)`

The first validates sequence/structure length and uses ViennaRNA’s thermodynamic model to compute free energy for a given dot-bracket structure. The second computes the absolute energy gap and relative error percentage compared to a reference MFE value.

### Verified example
For `GCGCAUACGC` and `(((....)))`, the energy evaluator returned:
- `-1.3`

This exactly matches the toy MFE output from ViennaRNA after rounding.

This is important because it confirms the project has a reliable path for comparing future quantum-predicted structures against classical reference energies.

## 7. Current Technical Interpretation

At this stage, the project has a clear three-part foundation:

1. **Classical baseline**: ViennaRNA folding works.
2. **Energy scoring**: custom evaluation and gap metrics work.
3. **Optimization formulation**: QUBO construction works on toy sequences.

This is enough to say that the core scientific framing of the project is validated. What remains is solver integration and scaling analysis.

## 8. Gaps Remaining

The following are still not implemented or not yet connected end-to-end:
- A general benchmark harness that ties together sequence input, MFE folding, candidate scoring, and result logging.
- QUBO-to-solver integration.
- CVaR-VQE or QAOA implementation.
- Resource-scaling analysis.
- Structured plots and final report artifacts.

These missing pieces are expected for the next days of work.

## 9. Recommended Next Steps

The next phase should proceed slowly and in order:

1. Create a small classical benchmark harness.
2. Refactor the current scripts into reusable functions if needed.
3. Add a simple decoding path from a candidate bitstring to a dot-bracket structure.
4. Move to solver implementation.
5. Start measuring scaling behavior on short sequences.

## 10. Summary

The first two days of the project are in good shape:
- Classical MFE benchmarking is working.
- Energy evaluation is working.
- QUBO formulation is working.
- Toy sequence tests match expectations.

The repo still needs to grow toward the ideal structure listed in the README, but the foundational pieces are now verified and ready for the next stage.