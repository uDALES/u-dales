# MPI Averaging Regression

This regression compares MPI-sensitive dumped fields and selected solver log
diagnostics between two refs on compact committed cases.

It is built on the same short-run pattern as the processor-boundary parity
test. Instead of one heavier statistics case, it runs two compact cases across
four processor decompositions:

- no-tree case `100`
- tree case `526`
- `serial`, `x_split`, `y_split`, and `xy_split`
- reconstructed global `tdump` / `treedump` fields
- optional solver log diagnostics from `checksim`

The goal is to catch behavioural drift in the MPI averaging and reduction
layer without requiring a direct solver-side runmode on both compared refs.

This regression is intentionally indirect. It does not call the operators one
by one. Instead, it compares normal solver outputs that are common across both
branches.

On this branch, the direct companion coverage is:

- `runmode = 1005` (`TEST_MPI_OPERATORS`) in `src/tests.f90`
- `tests/integration/mpi_operators/run_test.sh`

Keep this regression as well, because a new runmode cannot be assumed to exist
in older reference commits.

Typical usage:

```bash
python tests/regression/mpi_averaging_regression/run_test.py HEAD^ HEAD Debug
python tests/regression/mpi_averaging_regression/run_test.py HEAD^ HEAD Debug --ci
python tests/regression/mpi_averaging_regression/run_test.py \
  HEAD^ HEAD Debug \
  --reuse-current-build \
  --cases 100,526 \
  --configs xy_split
```
