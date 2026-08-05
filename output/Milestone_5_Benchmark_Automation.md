# Milestone 5 Report — Benchmark Automation and Solver Comparison

## Summary

Milestone 5 focused on automating benchmark execution for the RNA QUBO solver pipeline. The project now supports both exact and heuristic solving paths, enabling validation on short sequences and scalable testing on longer sequences.

## Methods

The benchmark script runs a fixed set of RNA sequences through the existing pipeline:

1. Build the QUBO model.
2. Select the solver path through the hybrid solver.
3. Decode the resulting bitstring into RNA secondary structure.
4. Evaluate the decoded structure using the classical energy function.
5. Record runtime and solver metadata in CSV format.

The benchmark workflow was implemented in `scripts/run_day5_benchmarks.py`. Results are written to `output/rna_solver_results_day5.csv`.

## Benchmark Set

The following sequences are included in the Milestone 5 benchmark batch:

- `GCGCAUACGC`
- `GCGCGUACGC`
- `GCAUCGUAGC`
- `CUACGGCGCGGCGCCCUUGGCGA`
- `GCGCGCGAAAUUCGCGCG`

Short sequences are expected to use the exact brute-force solver, while longer sequences are routed through the heuristic solver.

## Results

The automated benchmark script records the following fields for each sequence:

- timestamp
- sequence
- length
- method
- decoded structure
- active pairs
- energy
- QUBO cost
- runtime in seconds
- notes

This output is stored in `output/rna_solver_results_day5.csv` for later comparison and reporting.

## Observations

The pipeline now supports repeatable batch execution without manual benchmarking. This improves consistency across experiments and makes it easier to compare exact and heuristic behavior over time.

The solver architecture remains modular:
- `quantum/solver.py` retains the exact brute-force baseline.
- `quantum/heuristic_solver.py` provides the fast approximate path.
- `quantum/hybrid_solver.py` chooses between them automatically.
- `scripts/run_solver.py` remains available for single-sequence runs.

## Limitations

The exact solver remains practical only for short sequences. For larger sequences, heuristic results are faster but are not guaranteed to match the classical minimum free energy structure.

## Next Steps

Future work should focus on:

- expanding the benchmark set,
- comparing runtime and accuracy across methods,
- adding plots or summary tables,
- and incorporating more ViennaRNA reference comparisons.

## Status

Milestone 5 introduces automated benchmarking and reproducible results tracking. The project is now set up for systematic experimentation rather than one-off manual runs.