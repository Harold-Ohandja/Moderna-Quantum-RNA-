# RNA QUBO Work Log

## Current status
The RNA QUBO pipeline is now split into two paths:
- Exact brute-force solver for short sequences.
- Heuristic solver for larger sequences.

The decoder and energy evaluator are working end-to-end.

## Files added
- `quantum/heuristic_solver.py`
- `quantum/hybrid_solver.py`
- `scripts/run_solver.py`

## Validation completed
### Short sequence benchmark
- Sequence: `GCGCAUACGC`
- Method: brute force
- Decoded structure: `(((....)))`
- Active pairs: `[(0, 9), (1, 8), (2, 7)]`
- Energy: `-1.3`
- QUBO cost: `-6.0`

### Longer sequence benchmark
- Sequence: `CUACGGCGCGGCGCCCUUGGCGA`
- Method: heuristic
- Decoded structure: `(((.((((...)).)).)))...`
- Active pairs: `[(0, 19), (1, 18), (2, 17), (4, 15), (5, 14), (6, 12), (7, 11)]`
- Energy: `-2.2`
- QUBO cost: `-13.0`

## Notes
- Brute force is appropriate only for small instances.
- The heuristic solver prevents long runtimes on larger sequences.
- The CSV file will be useful for tracking results across future days and experiments.

## Next likely steps
- Add more sequence rows to the CSV as experiments continue.
- Compare heuristic results against ViennaRNA benchmarks where possible.
- Prepare summary figures and a short report section later.
