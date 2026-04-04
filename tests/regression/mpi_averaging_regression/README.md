# MPI Averaging Regression

This regression protects the recent MPI/reduction refactor against a chosen
reference branch or commit.

Current local reference for the active MPI refactor check:

- `6eddf5fa2600a5eec31c4a119495d2195d594c6b`
  `Address processor-boundary tolerance and Codex review fixes`

The automated `experimental: mpi averaging regression` suite entry uses this
commit explicitly. It does not use the generic `master` default from
`tests/run_tests.py`, because old `master` crashes on committed case `100`
during IBM initialization in `Debug`.

It is built on the same one-step execution pattern as the processor-boundary
parity test. Instead of one heavier statistics case, it runs two compact cases
across four processor decompositions:

- no-tree case `100`
- tree case `526`
- `serial`, `x_split`, `y_split`, and `xy_split`
- solver log diagnostics from `checksim`
- reconstructed global `tdump` / `treedump` fields

The goal is not broad solver acceptance. The goal is to catch behavioural drift
in the core reduction operators that were recently reorganized across
`architecture`, `modmpi`, `operators`, and `ibm`.

This regression is intentionally indirect. It does not call the operators one
by one. Instead, it compares existing solver outputs that are already common
across both branches. That is what makes the test usable against older commits:
the shared interface is the normal namelist plus the normal dumped fields.

In particular, the comparison is intended to exercise:

- `global_max`
- `global_sum`
- `reduce_xy_sum`
- `reduce_yz_sum`
- `avg_xy_fluid`
- `avg_y_fluid`
- `sum_x_fluid`
- `sum_y_fluid`

Coverage note:

- `global_max` and `global_sum` are only checked when `--compare-logs` is used
- the field operators are exercised through the existing startup, statistics,
  output, and decomposition-sensitive solver paths that produce `tdump` and
  `treedump`
- this is a branch-stable regression, not a direct unit test of each operator

Complementary direct test:

- current uDALES now also has `runmode = 1005` (`TEST_MPI_OPERATORS`) in
  `tests.f90`
- that direct solver-side test uses the committed Xie/Castro geometry from case
  `100` and checks `global_sum`, `global_max`, `reduce_xy_sum`,
  `reduce_yz_sum`, `avg_xy_fluid`, `avg_y_fluid`, `sum_x_fluid`, and `sum_y_fluid`
  against explicit masked reference calculations
- keep this regression as well, because a new runmode cannot be assumed to
  exist in older reference commits

The selected outputs cover the main staggered locations used by the fluid-only
averaging operators:

- `LOC_C`
- `LOC_U`
- `LOC_V`
- `LOC_W`
- `LOC_UV`
- `LOC_WU`
- `LOC_VW`

The script follows the same reference-build pattern as the existing
`new_vegetation_module_against_v2.2` regression, but uses the
processor-boundary-style one-step harness for the actual solver runs:

- build a reference worktree for `branch_a`
- use the workspace build for the current side unless `--ci` is requested
- run both branches on patched temporary copies of the compact cases
- reconstruct global tiled output fields
- compare reconstructed fields that directly exercise the MPI averaging/reduction layer
- optionally compare selected solver diagnostics with `--compare-logs`

For fast local iteration after the reference build already exists, the script
can reuse an existing reference build directory and skip rebuilding the current
workspace executable:

- `--reference-build-dir <path>` reuses an existing reference build directory
  containing `u-dales`
- `--reuse-current-build` reuses `build/<type>/u-dales` from the current
  workspace instead of rebuilding it
- `--cases <ids>` limits the run to a comma-separated subset such as `100,526`
- `--configs <labels>` limits the decomposition matrix to labels such as
  `xy_split` or `serial,x_split`
- `--scratch-root <path>` places worktrees and run outputs under a persistent
  directory of your choice
- run/worktree directories are kept by default so failures leave artifacts
  behind; pass `--cleanup` only once the test is stable and you no longer need
  the outputs
- when preserving artifacts, prefer a fresh scratch root to avoid quota/open
  failures from reusing an oversized debug directory

Typical usage:

```bash
python tests/regression/mpi_averaging_regression/run_test.py HEAD^ HEAD Debug
python tests/regression/mpi_averaging_regression/run_test.py \
  6eddf5fa2600a5eec31c4a119495d2195d594c6b HEAD Debug
python tests/regression/mpi_averaging_regression/run_test.py HEAD^ HEAD Debug --ci
python tests/regression/mpi_averaging_regression/run_test.py \
  HEAD^ HEAD Debug \
  --reference-build-dir /path/to/reference/build/debug \
  --reuse-current-build
python tests/regression/mpi_averaging_regression/run_test.py \
  HEAD^ HEAD Debug \
  --reference-build-dir /path/to/reference/build/debug \
  --reuse-current-build \
  --scratch-root /tmp/mpi-averaging-debug
python tests/regression/mpi_averaging_regression/run_test.py \
  HEAD^ HEAD Debug \
  --reuse-current-build \
  --cases 100,526 \
  --configs xy_split \
  --scratch-root /path/to/fresh/scratch
python tests/regression/mpi_averaging_regression/run_test.py \
  HEAD^ HEAD Debug \
  --reference-build-dir /path/to/reference/build/debug \
  --reuse-current-build \
  --compare-logs \
  --cleanup
```
