
# Moderna Challenge Report: Milestone 2 – Classical Benchmark Harness

**Project:** Optimization of mRNA Secondary Structure Prediction Using Quantum Computing  
**Challenge:** WISER Global Quantum+AI Program 2026 — Moderna Challenge  
**Team:** Pushkar Kumar and Harold  

This report explains, in non-technical language where possible, what was achieved in the project’s "Milestone 2": building a small but complete benchmark tool that connects classical RNA structure prediction, a custom energy scorer, and a mathematical model used for quantum optimization.

---

## 1. Why We Need a Benchmark Harness

The Moderna challenge is about predicting how a strand of RNA folds into a secondary structure (which bases pair, which stay unpaired) and then exploring whether quantum or quantum-inspired methods can find good folds efficiently.

To judge any quantum approach fairly, we need a **trusted classical baseline**:
- A tool that, given an RNA sequence, can produce the best-known structure and its energy using a mature classical library.
- A way to score other candidate structures and measure how far they are from this classical reference.
- A way to connect those structures to the optimization model we use on the quantum side.

Milestone 2 is about creating that baseline tool — a **benchmark harness**.

The harness is a small program that:
1. Takes a sequence.
2. Computes the classical best (Minimum Free Energy, or MFE) structure and energy.
3. Evaluates the energy of that structure with a custom scoring function.
4. Builds a mathematical model (QUBO) capturing possible base-pair choices.
5. Saves all these results in a structured file for later analysis.

---

## 2. What the Harness Actually Does

The harness lives in `classical/benchmark.py` and is run from the project root with:

```bash
python -m classical.benchmark
```

When it runs, it processes two short RNA sequences (for now) and performs three main tasks:

### 2.1 Classical folding with ViennaRNA

For each sequence, the harness calls the **ViennaRNA** library, a standard tool in computational biology, to compute:
- The predicted secondary structure in **dot-bracket notation** (characters like `(`, `)`, and `.` that describe which bases pair).
- The associated **free energy** in kilocalories per mole (kcal/mol).

For example, for the 10-nucleotide sequence `GCGCAUACGC`:
- ViennaRNA predicts the structure `(((....)))`.
- The energy is about `-1.30 kcal/mol` — a negative value indicating a relatively stable folded structure.

For the 8-nucleotide sequence `AUGCAUGC`:
- ViennaRNA predicts the structure `........` (all bases unpaired).
- The energy is `0.0 kcal/mol` — indicating that, under this model, the best configuration is essentially "no folding." [github_mcp_direct:2]

### 2.2 Energy evaluation and consistency check

The project includes a custom **energy evaluation module** that:
- Checks that the sequence and structure lengths match.
- Uses ViennaRNA’s thermodynamic model internally to compute the energy of a given structure.
- Provides a function to compute the energy gap between a candidate structure and the reference MFE.

In the harness, we use this module to recompute the energy of the MFE structure for each sequence. This serves as a consistency check:
- For `GCGCAUACGC` and `(((....)))`, the custom evaluator returns `-1.3 kcal/mol`, matching the ViennaRNA result after rounding.
- For `AUGCAUGC` and `........`, the evaluator returns `0.0 kcal/mol`, again matching ViennaRNA.

This tells us that our scoring path is reliable and can later be used to evaluate structures produced by quantum or classical optimizers.

### 2.3 QUBO formulation for optimization

The harness then builds a **QUBO (Quadratic Unconstrained Binary Optimization)** model for each sequence using the `quantum/qubo.py` module.

Conceptually:
- Each potential base pair that is physically allowed (e.g., A–U, C–G, and G–U, respecting a minimum loop length) becomes a **binary variable**: 1 means "this pair is present"; 0 means "this pair is absent".
- The QUBO model assigns energy contributions to single pairs, bonuses for stacked pairs (pairs that sit on top of each other in the structure), and penalties for invalid combinations such as:
  - A base attempting to pair with multiple partners.
  - Crossing pairs (pseudoknots) that break the nested-structure assumption.

For `GCGCAUACGC`:
- There are 7 candidate base pairs.
- The QUBO has 26 non-zero terms.

For `AUGCAUGC`:
- There are 3 candidate base pairs.
- The QUBO has 6 non-zero terms.

This shows how the optimization model scales with different sequences: more possible pairings yield more variables and a richer QUBO.

---

## 3. What the Harness Prints

When `python -m classical.benchmark` runs, it prints a summary for each sequence. For example:

```text
=== Benchmark ===
Sequence       : GCGCAUACGC
Length         : 10
MFE structure  : (((....)))
MFE energy     : -1.2999999523162842 kcal/mol
Eval energy    : -1.3 kcal/mol
Candidate pairs: [(0, 5), (0, 7), (0, 9), (1, 8), (2, 7), (2, 9), (3, 8)]
Num pairs      : 7
QUBO size      : 26 non-zero terms

=== Benchmark ===
Sequence       : AUGCAUGC
Length         : 8
MFE structure  : ........
MFE energy     : 0.0 kcal/mol
Eval energy    : 0.0 kcal/mol
Candidate pairs: [(0, 5), (1, 6), (2, 7)]
Num pairs      : 3
QUBO size      : 6 non-zero terms
```

In natural language, for each sequence the harness answers:
- What the classical model thinks is the best fold.
- How stable that fold is (energy).
- How many valid base-pair options exist.
- How complex the optimization model is (size of the QUBO).

---

## 4. What the Harness Saves to Disk

Beyond printing, the harness now writes a structured JSON file at:

```text
data/benchmark/reference_mfe.json
```

The file looks like this (simplified):

```json
{
  "GCGCAUACGC": {
    "sequence": "GCGCAUACGC",
    "length": 10,
    "mfe_structure": "(((....)))",
    "mfe_energy": -1.2999999523162842,
    "eval_energy": -1.3,
    "candidate_pairs": [[0, 5], [0, 7], [0, 9], [1, 8], [2, 7], [2, 9], [3, 8]],
    "num_pairs": 7,
    "qubo_size": 26
  },
  "AUGCAUGC": {
    "sequence": "AUGCAUGC",
    "length": 8,
    "mfe_structure": "........",
    "mfe_energy": 0.0,
    "eval_energy": 0.0,
    "candidate_pairs": [[0, 5], [1, 6], [2, 7]],
    "num_pairs": 3,
    "qubo_size": 6
  }
}
```

This JSON file is important because it:
- Captures classical reference data in a reusable format.
- Can be read by future scripts to compare quantum or optimization results against known values.
- Provides a simple starting point for plotting and scaling analysis.

In practical terms, this is now a "ground truth" file for short sequences.

---

## 5. How This Fits the Plan

The internal project plan defined Milestone 2 as:

> "Build classical benchmark harness: RNA.fold() for MFE + eval_structure() for scoring any candidate. Validate on 2–3 test sequences."  
> *(Moderna_Challenge_Team_Plan.pdf)*

The work completed matches this intent:
- The harness uses `RNA.fold()` (ViennaRNA) to obtain MFE structures and energies.
- It uses `evaluate_structure_energy()` to score these structures and confirm consistency.
- It runs and validates on two short sequences, with room for adding the 44-nt official sequence next.
- It saves results to a `data/benchmark/` directory as structured JSON.

This means "Milestone 2" — having a usable classical benchmark tool — is effectively complete.

---

## 6. Why This Matters for Non-Technical Stakeholders

For someone less familiar with the underlying math or quantum computing, the key points are:

- We now have a **reliable way to check our work**. Any future quantum or optimization method can be tested against these classical references.
- The benchmark harness provides both:
  - A human-readable summary (printed in the console), and
  - A machine-readable record (JSON file) that can feed dashboards, reports, or visualizations.
- The harness is small, transparent, and easily audited; it uses a well-known library (ViennaRNA) and simple data structures.

This foundation reduces risk: instead of experimenting blindly, we have a clear yardstick for correctness.

---

## 7. Next Steps After Milestone 2

With the classical benchmark harness in place, the next stages of the project will focus on:

1. **Adding more sequences** — including the 44-nt official example from the challenge — to the benchmark JSON.
2. **Implementing solver logic** — starting with small classical or quantum-inspired solvers that operate on the QUBO model.
3. **Decoding solver outputs** — translating solver bitstrings into base-pair selections, building dot-bracket structures, and scoring them with the energy evaluator.
4. **Comparing against reference MFE** — using the JSON data and gap metrics to quantify performance.
5. **Scaling analysis** — using this harness as a basis for measuring how resource usage grows with sequence length.

These steps will build directly on the harness described here, without reinventing the classical baseline.

---

## 8. Summary

Milestone 2’s work has produced a small but complete classical benchmark harness for RNA secondary structure. It:
- Uses ViennaRNA to predict optimal folds.
- Uses a custom module to evaluate energies and prepare gap metrics.
- Uses a QUBO formulation to represent candidate base-pair choices.
- Logs all results into a structured JSON file.

This harness now serves as the "reference layer" for the rest of the project. Any future quantum or quantum-inspired method can be judged against the data it produces, which is exactly what is needed for a rigorous and credible submission to the Moderna WISER challenge.
