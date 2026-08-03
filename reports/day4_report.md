# Day 4 Report — RNA QUBO Solver Validation

## Overview
Day 4 focused on validating the RNA folding pipeline end-to-end and making the project usable beyond the smallest toy sequence. The main goal was to confirm that the QUBO formulation, solver, decoder, and classical energy evaluator all work together correctly, while also preparing the codebase for larger instances that are too expensive for brute-force search.

## Objectives
The goals for the day were:
- Verify that the exact brute-force solver still works on the reference toy benchmark.
- Check that the decoded bitstring maps back to the expected RNA secondary structure.
- Confirm that the energy evaluator returns the expected classical folding energy.
- Add a faster solver path for larger sequences where brute force becomes impractical.
- Create reproducible experiment outputs for future reporting and comparison.

## Validation of the exact solver
The baseline sequence `GCGCAUACGC` was used to validate the exact brute-force solver path. The solver returned the bitstring `bits: [0, 0, 1, 1, 1, 0, 0]`, which decoded to the structure `(((....)))` with active pairs `[(0, 9), (1, 8), (2, 7)]`.

The classical energy evaluator returned `-1.3`, and the QUBO cost was `-6.0`. This confirmed that the exact pipeline was functioning correctly and that the decoder maps the QUBO output to the expected RNA structure.

## Validation of the larger-sequence path
A longer sequence, `CUACGGCGCGGCGCCCUUGGCGA`, was used to test the larger-instance workflow. A brute-force attempt on this sequence was too slow and had to be interrupted, which confirmed that exhaustive search does not scale well for this input size.

To solve this, a separate heuristic solver path was added. Using the new heuristic route, the sequence produced the structure `(((.((((...)).)).)))...` with active pairs `[(0, 19), (1, 18), (2, 17), (4, 15), (5, 14), (6, 12), (7, 11)]`. The reported energy was `-2.2`, and the QUBO cost was `-13.0`.

This result confirmed that the longer-sequence path works end-to-end and that the project can now handle both small exact benchmarks and larger approximate experiments.

## Code changes made
The following files were added during Day 4:
- `quantum/heuristic_solver.py`
- `quantum/hybrid_solver.py`
- `scripts/run_solver.py`

The existing exact solver in `quantum/solver.py` was kept as the baseline. The new solver logic was separated so the project can later report exact and heuristic results independently.

## Experiment tracking
To make future work easier to track, two output files were created in the project root `output/` folder:
- `output/rna_solver_results.csv`
- `output/rna_current_work_log.md`

The CSV file is intended to store one row per experiment, including sequence, solver method, decoded structure, active pairs, energy, QUBO cost, and notes. The markdown log summarizes the current state of the project and can be reused in the final report.

## Important debugging note
While implementing the heuristic path, one issue appeared because the QUBO was stored as a sparse dictionary rather than a dense matrix. The heuristic solver originally assumed row indexing and raised a `KeyError: 0`. This was fixed by updating the heuristic cost evaluation to iterate through `(i, j)` QUBO keys directly.

## Current project status
By the end of Day 4:
- The exact solver works on the toy benchmark.
- The decoder works and reproduces the expected dot-bracket structure.
- The energy evaluator works and returns the expected value.
- The heuristic solver works on a longer sequence.
- Result logging files are in place for future experiments and reporting.

## Next steps
The next work items are:
- Add more benchmark sequences to the CSV log.
- Compare heuristic results with ViennaRNA more systematically.
- Improve runtime and quality of the heuristic search if needed.
- Prepare summary tables and figures for the final report.
- Continue building toward the quantum or hybrid solver implementation.

## Day 4 summary
Day 4 was successful because it moved the project from a single validated toy example to a more complete experimental workflow. The solver pipeline is now modular, reproducible, and ready for further comparison work.
