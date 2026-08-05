#  Optimization of mRNA Secondary Structure Prediction Using Quantum Computing 

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Qiskit](https://img.shields.io/badge/Qiskit-2.5%2B-purple)
![ViennaRNA](https://img.shields.io/badge/ViennaRNA-2.7%2B-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

> **WISER Global Quantum+AI Program 2026** | Challenge Moderna  
> **Repository:** [`Moderna-Quantum-RNA-`](https://github.com/Harold-Ohandja/Moderna-Quantum-RNA-)

---
## Team

- **Pushkar Kumar** — quantum/ML implementation lead: QUBO formulation, CVaR-VQE and QAOA solvers, notebooks
- **Harold Ohandja** — classical solver track: brute-force/heuristic baselines, benchmark automation, resource analysis

---

## Objectives

The main objectives of this project are:

- Understand the biological principles of RNA secondary structure.
- Generate classical benchmark structures using ViennaRNA.
- Formulate RNA folding as an optimization problem.
- Encode the optimization problem into a QUBO/Ising model.
- Implement a quantum (or quantum-inspired) optimization algorithm.
- Compare quantum predictions against classical MFE structures.
- Analyze the scalability and quantum resource requirements.

---

##  Executive Summary

Predicting mRNA secondary structure is a crucial computational bottleneck in the design of mRNA-based therapeutics and vaccines. Traditional tools, such as **ViennaRNA**, rely on dynamic programming algorithms to identify Minimum Free Energy (MFE) structures under standard thermodynamic models. 

This project explores a quantum optimization approach for mRNA folding by mapping base-pairing decisions to a Quadratic Unconstrained Binary Optimization (**QUBO**) model. Inspired directly by research from IBM Quantum and Moderna [[1](#1-alevras-et-al-2024)], we implement a **CVaR-VQE (Conditional Value at Risk Variational Quantum Eigensolver)** algorithm using a hardware-efficient ansatz. We benchmark the quantum predictions against classical MFE structures, evaluate energy gaps ($\Delta E$), and conduct a resource scaling analysis.

---

##  Theoretical Background & QUBO Formulation

### 1. Binary Variable Mapping
Instead of mapping individual nucleotides to qubits, variables represent **candidate base pairs** $(i, j)$ that satisfy physical constraints ($\vert{}j - i\vert{} \ge 4$ and Watson-Crick / Wobble base pairing rules $\text{A-U, G-C, G-U}$).

$$x_k = \begin{cases} 1 & \text{if candidate base pair } k = (i, j) \text{ forms} \\ 0 & \text{otherwise} \end{cases}$$

### 2. Objective Function & Hamiltonian
The total Hamiltonian combines thermodynamic binding energies with penalty constraints:

$$H_{\text{total}} = H_{\text{energy}} + H_{\text{overlap}} + H_{\text{pseudoknots}}$$

* **Thermodynamic Energy ($H_{\text{energy}}$):** Favors stable pairings ($\text{G-C} < \text{A-U} < \text{G-U}$) and accounts for stacking bonuses [[1](#1-alevras-et-al-2024)].
* **Overlap Constraint ($H_{\text{overlap}}$):** Penalizes configurations where a single base forms multiple bonds ($P_{\text{overlap}} \cdot x_{k_1} x_{k_2}$) [[1](#1-alevras-et-al-2024)].
* **Pseudoknot Exclusion ($H_{\text{pseudoknots}}$):** Penalizes crossing pairs $(i_1 < i_2 < j_1 < j_2)$ to maintain nested secondary structures [[1](#1-alevras-et-al-2024)].

### 3. CVaR-VQE Algorithm
Standard VQE minimizes expected energy over all samples, which can struggle in rugged combinatorial landscapes [[1](#1-alevras-et-al-2024)]. We employ **CVaR-VQE**, optimizing only the expectation value of the best $\alpha$-quantile ($\alpha = 0.1$) of sampled bitstrings, drastically accelerating convergence toward the ground state [[1](#1-alevras-et-al-2024)].

---

##  Repository Structure

```text
Moderna-Quantum-RNA-/
├── README.md                     # Project summary and documentation
├── LICENSE                       # MIT License
├── requirements.txt              # Python dependencies (pip)
├── environment.yml               # Conda environment configuration
│
├── data/
│   └── benchmark/
│       └── reference_mfe.json    # ViennaRNA reference outputs (auto-generated)
│
├── classical/
│   ├── generate_mfe.py           # ViennaRNA MFE structure generation
│   ├── evaluate_energy.py        # Energy gap (ΔE) calculator
│   ├── benchmark.py              # Classical benchmark runner (CLI: --sequence)
│   └── utils.py                  # Sequence generators & helpers
│
├── quantum/
│   ├── qubo.py                   # QUBO matrix formulation & penalties
│   ├── hamiltonian.py            # Ising Hamiltonian mapping
│   ├── decode.py                 # Bitstring -> dot-bracket structure decoder
│   ├── vqe_solver.py             # CVaR-VQE solver, primary quantum method (CLI)
│   ├── qaoa_solver.py            # QAOA / CVaR-QAOA solver (CLI)
│   ├── resource_analysis.py      # Qubit/circuit-depth scaling analysis
│   ├── solver.py                 # Classical baseline: exact brute-force QUBO solver
│   ├── heuristic_solver.py       # Classical baseline: greedy local-search solver
│   └── hybrid_solver.py          # Auto-switches between the two classical baselines
│
├── notebooks/
│   ├── 01_Classical_Benchmark.ipynb  # ViennaRNA baseline + QUBO formulation walkthrough
│   ├── 02_CVaR_VQE_Solver.ipynb      # Primary quantum solver, with results & limitations
│   ├── 03_QAOA_Comparison.ipynb      # QAOA vs. CVaR-VQE, encoding/tradeoff discussion
│   └── 04_Full_Comparison.ipynb      # All methods side by side — start here
│
├── scripts/
│   ├── run_solver.py             # Run the classical baseline on one sequence (CLI)
│   ├── run_milestone5_benchmarks.py  # Batch classical baseline runs -> CSV
│   └── run_full_comparison.py    # All methods, all sequences, one combined CSV
│
├── output/                       # CSV/markdown results logs (generated by the scripts above)
└── reports/                      # Milestone-by-milestone progress reports
```

Note: `quantum/solver.py`, `heuristic_solver.py`, and `hybrid_solver.py` are classical
baselines (exact brute-force and greedy heuristic search), not quantum or
quantum-inspired algorithms. They live in `quantum/` for now alongside the actual
quantum solvers for convenience, and exist to give a second, independent classical
comparison point on the exact same QUBO the quantum solvers use.
## 🛠️ Getting Started

### Option A: No terminal needed (recommended)

Every notebook in `notebooks/` runs standalone on **Google Colab**:

1. Open the notebook file on GitHub (e.g. `notebooks/04_Full_Comparison.ipynb`).
2. Click the "Open in Colab" badge at the top (or open [colab.research.google.com](https://colab.research.google.com), then File → Open notebook → GitHub, and paste this repo's URL).
3. Run the first cell (labeled "Setup"). It installs everything needed and clones the repo automatically.
4. Then: Runtime → Run all.

That's the whole process, no installation, no command line. Results (tables, charts)
appear inline in the notebook.

**If you only run one notebook, run `04_Full_Comparison.ipynb`** — it puts every
method side by side on the same sequences.

### Option B: Running locally (for anyone comfortable with a terminal)

```bash
git clone https://github.com/Harold-Ohandja/Moderna-Quantum-RNA-.git
cd Moderna-Quantum-RNA-
pip install -r requirements.txt
```

Or with Conda:
```bash
conda env create -f environment.yml
conda activate moderna-quantum-rna
```

Then either open the notebooks locally with `jupyter notebook`, or run any script
directly from the repo root, for example:

```bash
# Classical ViennaRNA benchmark on one sequence
python classical/benchmark.py --sequence "GCGCAUACGC"

# Primary quantum solver (CVaR-VQE)
python quantum/vqe_solver.py --sequence "GCGCAUACGC" --alpha 0.1

# Comparison quantum solver (QAOA)
python quantum/qaoa_solver.py --sequence "GCGCAUACGC"

# Everything at once (classical baseline + both quantum solvers), one combined table
python scripts/run_full_comparison.py
```

**What to expect:** each command prints the decoded structure, its energy, and how
it compares to ViennaRNA's reference MFE. `run_full_comparison.py` is the most
informative single command, it takes under 15 seconds and writes a combined CSV to
`output/full_comparison_results.csv` alongside the printed table.

##  License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

##  References & Bibliography

### 1. Alevras et al. (2024)
Alevras, A., et al. *mRNA secondary structure prediction using utility-scale quantum computers*.  
[arXiv:2405.20328 [quant-ph]](https://arxiv.org/abs/2405.20328)

### 2. ViennaRNA Package
Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011).  
*ViennaRNA Package 2.0*. Algorithms for Molecular Biology, 6(1), 1-14.  
[DOI: 10.1186/1748-7188-6-26](https://doi.org/10.1186/1748-7188-6-26)

### 3. WISER Global Quantum+AI Program
WISER Program 2026 - Challenge Moderna: *Quantum Optimization for Biological Sequences*.


##  Acknowledgments

We express our gratitude to the **WISER Global Quantum+AI Program 2026**, **Moderna**, and **IBM Quantum** for providing the challenge framework, reference models, and quantum computing resources.
