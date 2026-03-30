# Sparse IJK Test

This case validates the IBM sparse-coordinate input reader used by the solver.

- Entry point: `run_test.sh`
- Run mode: `1004`
- Default MPI size: `4`

The test compares the generic `read_sparse_ijk()` path against the established
IBM initialization path and fails if rank-local counts, indices, or coordinates
do not match.
