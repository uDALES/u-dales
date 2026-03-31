# Tree Sparse Compare Test

This test validates the sparse vegetation/tree forcing path used by the solver.

- Entry point: `run_test.sh`
- Run mode: `1005`
- Default MPI size: `32`
- Shared case fixture: `tests/cases/525/`

The test compares the sparse vegetation forcing against the established
block/tree forcing implementation and fails if the resulting source terms differ.

The runnable test directory keeps only the harness. The base case inputs are
staged from `tests/cases/525/` into a temporary run directory before execution.
