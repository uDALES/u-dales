# Tree Sparse Compare Test

This test validates the sparse vegetation/tree forcing path used by the solver.

Status:

- classification: `experimental`
- role: solver-facing integration test
- expected in the experimental dispatcher path for now
- shared fixture owner: `tests/cases/525/`

- Entry point: `run_test.sh`
- Run mode: `1005`
- Default MPI size: `32`
- Shared case fixture: `tests/cases/525/`

Prerequisites:

- MPI launcher/runtime available in `PATH`
- built `u-dales` executable available at the path expected by `run_test.sh`
- project module stack loaded where required by the local environment

The test compares the sparse vegetation forcing against the established
block/tree forcing implementation and fails if the resulting source terms differ.

The runnable test directory keeps only the harness. The base case inputs are
staged from `tests/cases/525/` into a temporary run directory before execution.

Run it directly with:

```bash
cd tests/integration/tree_sparse_compare
./run_test.sh
```
