#  Optimization of mRNA Secondary Structure Prediction Using Quantum Computing 

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Qiskit](https://img.shields.io/badge/Qiskit-1.0%2B-purple)
![ViennaRNA](https://img.shields.io/badge/ViennaRNA-2.5%2B-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

> **WISER Global Quantum+AI Program 2026** | Challenge Moderna  
> **Repository:** [`Moderna-Quantum-RNA`](https://github.com/Harold-Ohandja/Moderna-Quantum-RNA)

---
## Team

- **Team Member 1:** ...
- **Team Member 2:** ...
- **Team Member 3:** ...

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
Moderna-Quantum-RNA/
├── README.md                 # Project summary and documentation
├── LICENSE                   # MIT License
├── requirements.txt          # Python dependencies
├── environment.yml           # Conda environment configuration
│
├── data/
│   ├── sequences/            # Input RNA sequences (.fasta / .json)
│   ├── benchmark/            # ViennaRNA reference outputs
│   └── results/              # Optimization and scaling experiment logs
│
├── classical/
│   ├── generate_mfe.py       # ViennaRNA MFE structure generation
│   ├── evaluate_energy.py    # Energy gap (ΔE) calculator
│   ├── benchmark.py          # Full classical benchmark suite
│   └── utils.py              # Sequence generators & helpers
│
├── quantum/
│   ├── qubo.py               # QUBO matrix formulation & penalties
│   ├── hamiltonian.py        # Ising Hamiltonian mapping
│   ├── qaoa_solver.py        # QAOA solver implementation
│   ├── vqe_solver.py         # CVaR-VQE solver (Two-Local Ansatz)
│   └── resource_analysis.py # Scaling analysis (Qubits, depth, runtime)
│
├── notebooks/
│   ├── Benchmark.ipynb       # Interactive classical benchmark demo
│   └── QAOA.ipynb            # Quantum optimization & convergence demo
│
├── figures/                  # Generated plots (Scaling, Energy Gap, Circuits)
├── report/                   # LaTeX report files (report.tex, report.pdf)
└── presentation/             # Slide deck files (slides.tex, slides.pdf)
```
 
## 🛠️ Installation & Setup

### Option 1: Using Pip
```bash
# Clone the repository
git clone https://github.com.git
cd Moderna-Quantum-RNA

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Using Conda
```bash
conda env create -f environment.yml
conda activate moderna-quantum-rna
```

##  Usage

### 1. Run Classical ViennaRNA Benchmark
```bash
python classical/benchmark.py --sequence "AUGCAU..."
```

### 2. Run Quantum CVaR-VQE Optimization
```bash
python quantum/vqe_solver.py --sequence "AUGCAU..." --alpha 0.1
```

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
