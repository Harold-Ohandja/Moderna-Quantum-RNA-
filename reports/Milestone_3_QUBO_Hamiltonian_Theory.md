
# Moderna Challenge Report: Milestone 3 – Optimization Theory & Mapping

**Project:** Optimization of mRNA Secondary Structure Prediction Using Quantum Computing  
**Challenge:** WISER Global Quantum+AI Program 2026 — Moderna Challenge  
**Team:** Pushkar Kumar and Harold  

This report explains the theoretical foundation for the "Milestone 3" stage of the project: how the RNA folding problem is framed as a QUBO optimization model, how that QUBO is mapped to an Ising Hamiltonian for quantum algorithms such as CVaR-VQE, and how solutions (bitstrings) are decoded back into RNA secondary structures.

The goal is to present these ideas in a way that is understandable to non-specialists while still being precise enough for technical reviewers.

---

## 1. From RNA Folding to Binary Decisions

### 1.1 RNA secondary structure in simple terms

An RNA molecule is a chain of bases (A, U, C, G). In solution, parts of the chain can pair with each other:
- A pairs with U.
- C pairs with G.
- G can pair with U ("wobble" pairing).

These pairings cause the strand to fold into a secondary structure. Some bases are paired, some are unpaired. The folding pattern affects stability and function.

Standard tools like ViennaRNA predict a **Minimum Free Energy (MFE)** structure: the fold that makes the molecule most stable under a given thermodynamic model.

### 1.2 Turning the problem into binary choices

For optimization, instead of trying to directly choose a full structure, we break the problem into many small yes/no decisions:
- For each **physically allowed base pair** `(i, j)` (base at position `i` pairing with base at `j`), define a binary variable:

  - `x_k = 1` if that pair forms in the structure.
  - `x_k = 0` if it does not.

Here, `k` indexes the candidate pairs. The set of allowed pairs is determined by:
- base identity (A/U/C/G) and pairing rules.
- minimum loop length (to avoid unrealistically tight hairpins).

The job of the optimizer is to choose values for all `x_k` such that:
- the resulting structure is physically valid (no base pairs twice, no forbidden crossings), and
- the total energy is as low as possible.

---

## 2. QUBO Formulation (Quadratic Unconstrained Binary Optimization)

### 2.1 What is a QUBO?

A QUBO model uses binary variables `x_k ∈ {0, 1}` and defines a cost function:

\[
C_{\text{QUBO}}(x) = \sum_k Q_{kk} x_k + \sum_{k < l} Q_{kl} x_k x_l
\]

where:
- `Q_{kk}` are coefficients for individual variables (linear terms).
- `Q_{kl}` are coefficients for pairs of variables (quadratic terms).
- The goal is to **minimize** this cost.

In our context:
- Each variable `x_k` represents a candidate base pair.
- The matrix `Q` encodes:
  - energy contributions for forming pairs.
  - bonuses for well-formed stacking patterns.
  - penalties for invalid combinations (shared bases and pseudoknots).

### 2.2 Energy contributions

In the current simplified model:
- Diagonal terms `Q_{kk}` (individual pairs) are set to a negative value (e.g. `-1.0`).
  - This means forming a base pair (
`x_k = 1`) **lowers the energy**, which is thermodynamically favorable.

- Off-diagonal terms `Q_{kl}` can be:
  - A negative **stacking bonus** when two pairs `(i_1, j_1)` and `(i_2, j_2)` are adjacent in a nested way (`i_2 = i_1 + 1` and `j_2 = j_1 - 1`).
  - A positive **penalty weight** when two pairs share a base or form a crossing pattern (pseudoknot).

This structure captures the main trade-offs:
- Encourage stable nested pairings and stacking.
- Discourage invalid or biologically unrealistic configurations.

### 2.3 Constraints encoded as penalties

Physical constraints include:
- A nucleotide can pair at most once.
- Pairs should not cross (pseudoknot exclusion for the initial model).

Instead of enforcing these as “hard rules,” the QUBO model adds penalty terms:
- If two candidate pairs share a base or cross, `Q_{kl}` gets a positive penalty value.
- If both variables `x_k` and `x_l` are 1, the cost increases.

An optimizer that is minimizing `C_{QUBO}` will naturally prefer solutions that avoid penalties, thereby respecting the constraints.

This is consistent with the approach described in the IBM–Moderna reference work and the team’s Milestone 1b plan.

---

## 3. From QUBO to Ising Hamiltonian

Quantum optimization algorithms like VQE and CVaR-VQE work with **spin models**, often written in terms of Pauli Z operators (`Z`) on qubits. To use these methods, we need to convert the QUBO formulation into an **Ising Hamiltonian**:

\[
H = \text{constant} + \sum_k h_k Z_k + \sum_{k < l} J_{kl} Z_k Z_l
\]

where:
- `Z_k` acts on qubit `k` and has eigenvalues `+1` or `-1`.
- `h_k` are local field strengths.
- `J_{kl}` are coupling strengths between qubits `k` and `l`.

### 3.1 Binary variables to spin variables

We re-encode each binary variable `x_k ∈ {0, 1}` as a spin variable `z_k ∈ {+1, -1}` using:

\[
x_k = \frac{1 - z_k}{2}
\]

This mapping means:
- `x_k = 0` (pair not formed) ↔ `z_k = +1`.
- `x_k = 1` (pair formed)    ↔ `z_k = -1`.

This is standard in QUBO→Ising conversions.

### 3.2 Effect on linear terms

A diagonal QUBO term `Q_{kk} x_k` becomes:

\[
Q_{kk} x_k = Q_{kk} \cdot \frac{1 - z_k}{2} = \frac{Q_{kk}}{2} - \frac{Q_{kk}}{2} z_k
\]

This contributes:
- A **constant** energy offset `Q_{kk}/2`.
- A **linear** term `-Q_{kk}/2` multiplying `z_k`.

In the Hamiltonian:
- `z_k` is replaced by the Pauli Z operator `Z_k`.
- So the term adds `-Q_{kk}/2 * Z_k` to the Hamiltonian, plus a constant.

### 3.3 Effect on quadratic terms

An off-diagonal term `Q_{kl} x_k x_l` becomes:

\[
x_k x_l = \left(\frac{1 - z_k}{2}\right)\left(\frac{1 - z_l}{2}\right)
= \frac{1}{4} (1 - z_k - z_l + z_k z_l)
\]

Thus:

\[
Q_{kl} x_k x_l = \frac{Q_{kl}}{4}(1 - z_k - z_l + z_k z_l)
\]

This contributes:
- A **constant** term: `Q_{kl}/4`.
- Two **linear** terms: `-Q_{kl}/4 * z_k` and `-Q_{kl}/4 * z_l`.
- A **quadratic** coupling term: `Q_{kl}/4 * z_k z_l`.

In operator form:
- `z_k` becomes `Z_k`.
- `z_l` becomes `Z_l`.
- `z_k z_l` becomes the tensor product `Z_k Z_l`.

So each QUBO coefficient `Q_{kl}` yields:
- adjustments to the constant energy.
- contributions to local fields `h_k` and `h_l`.
- one coupling `J_{kl}` for `Z_k Z_l`.

### 3.4 Collecting Ising coefficients

By applying these transformations to all QUBO entries, we build:

- `constant`: the overall energy shift.
- `h_k`: effective local field on each qubit.
- `J_{kl}`: coupling strength between qubits.

The final Hamiltonian is:

\[
H = \text{constant} + \sum_k h_k Z_k + \sum_{k < l} J_{kl} Z_k Z_l
\]

This form is ready for use in VQE, CVaR-VQE, or related quantum optimization algorithms.

In code, this mapping can be encapsulated in a helper such as `qubo_to_ising(qubo_dict)`, which directly consumes the existing QUBO matrix produced by `quantum/qubo.py`.

---

## 4. From Bitstrings to RNA Structures

### 4.1 Solver output: bitstrings

Whether the optimizer is classical or quantum, its output will typically be one or more bitstrings:
- A bitstring is a sequence of 0s and 1s (e.g., `[1, 0, 1, 0, 0, 1, 0]`).
- The `k`-th bit corresponds to variable `x_k`.

Given:
- `candidate_pairs` as a list of physically allowed base pairs.
- `bitstring` as the solver’s suggestion.

We interpret:
- `x_k = 1` → candidate pair `k` is formed.
- `x_k = 0` → candidate pair `k` is not formed.

### 4.2 Selected base pairs

From the bitstring, we obtain:

```python
selected_pairs = [candidate_pairs[k] for k, xk in enumerate(bitstring) if xk == 1]
```

This gives a subset of `(i, j)` pairs that constitute the proposed structure.

In a perfect solution, these pairs:
- do not share bases.
- do not cross (no pseudoknots).

### 4.3 Building a base-pair map

To translate the selected pairs into a structure string, we build a map from each position to its partner:

```python
pair_map = {i: None for i in range(N)}  # N = sequence length
for (i, j) in selected_pairs:
    pair_map[i] = j
    pair_map[j] = i
```

Now:
- If `pair_map[i]` is `None`, base `i` is unpaired.
- If `pair_map[i] = j`, base `i` is paired with base `j`.

### 4.4 Dot-bracket representation

Dot-bracket notation encodes the structure as a string of length `N` using:
- `'.'` for unpaired bases.
- `'('` at the left partner of a pair.
- `')'` at the right partner.

We build the dot-bracket string as:

```python
dot_bracket_chars = []
for i in range(N):
    j = pair_map[i]
    if j is None:
        dot_bracket_chars.append('.')
    elif i < j:
        dot_bracket_chars.append('(')
    else:
        dot_bracket_chars.append(')')

dot_bracket = ''.join(dot_bracket_chars)
```

For example, a nested structure might produce a string like `(((....)))`, matching ViennaRNA’s output for the toy sequence.

### 4.5 Optional validation checks

Although the QUBO penalties should discourage invalid patterns, the decoder can optionally check for:
- **Shared-base conflicts**: drop or flag pairs if a base tries to pair twice.
- **Crossings (pseudoknots)**: detect index patterns like `i1 < i2 < j1 < j2` and handle them explicitly.

These checks can help diagnose solver behavior and ensure that decoded structures are meaningful.

---

## 5. Integrating with the Energy Evaluator

Once a dot-bracket structure is obtained from a solution bitstring:

1. The energy evaluator `evaluate_structure_energy(sequence, dot_bracket)` computes its free energy using ViennaRNA’s model.
2. The gap function `calculate_energy_gap(predicted_energy, reference_mfe_energy)` computes:
   - the absolute energy difference.
   - the relative error percentage.

This allows direct comparison between:
- the classical MFE structure and energy.
- the solver-proposed structure and energy.

These metrics feed into the challenge’s judging criteria regarding approximation quality and benchmarking.

---

## 6. Summary of Milestone 3 Theory

Milestone 3’s theoretical work establishes:
- A clear QUBO model for candidate base pairs, stacking, and constraints.
- A standard mapping from QUBO to an Ising Hamiltonian suitable for VQE/CVaR-VQE.
- A precise decoding pipeline from solver bitstrings back to dot-bracket RNA structures.

With these pieces defined, the project is ready to:
- Implement an initial solver (classical or CVaR-VQE) for small sequences.
- Decode and score solver outputs.
- Compare them against the classical reference benchmarks built in earlier days.

This theory lays the foundation for the next phase of the project, where actual optimization runs and scaling analyses will be performed.
