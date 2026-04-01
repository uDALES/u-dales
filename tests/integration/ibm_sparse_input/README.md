# IBM Sparse Input Test

This case validates the IBM sparse-coordinate input reader used by the solver.

Status:

- classification: `supported`
- role: solver-facing integration test
- expected in the supported dispatcher path
- shared fixture owner: `tests/cases/101/`

- Entry point: `run_test.sh`
- Run mode: `1004`
- Default MPI size: `4`
- Shared case fixture: `tests/cases/101/`

Prerequisites:

- MPI launcher/runtime available in `PATH`
- built `u-dales` executable available at the path expected by `run_test.sh`
- project module stack loaded where required by the local environment

The test compares the generic `read_sparse_ijk()` path against the established
IBM initialization path and fails if rank-local counts, indices, or coordinates
do not match.

The runnable test directory keeps only the test harness and the rank-local
golden output files (`solid_*_X_Y.txt` and `fluid_boundary_*_X_Y.txt`). The
base case inputs are staged from `tests/cases/101/` into a temporary run
directory before execution.

Run it directly with:

```bash
cd tests/integration/ibm_sparse_input
./run_test.sh
```
