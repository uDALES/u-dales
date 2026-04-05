# Processor-Boundary Parity Test

This diagnostic integration test checks whether committed cases are invariant to
the horizontal MPI decomposition used to run the solver.

It currently covers two fixtures:

- a no-tree Xie/Castro-style urban geometry case: `100`
- a vegetation case: `526`

The no-tree comparison uses case `100`, which is the committed Xie/Castro-style
urban geometry fixture and is the stronger long-term reference for
decomposition parity. The tree-forcing comparison uses case `526` as the small
committed vegetation fixture.

It stages the committed fixture for each variant into temporary run directories
and executes the same case with:

- `1 x 1` ranks
- `2 x 1` ranks
- `1 x 2` ranks
- `2 x 2` ranks

The test file has two variants:

- a no-tree comparison on case `100` that reconstructs `tdump.*.*.100.nc` and
  compares `ut`, `vt`, and `wt` between serial and decomposed runs
- a tree-forcing comparison on case `526` based on `treedump.*.*.526.nc`

Note:

- the no-tree `tdump` comparison allows a very small residual in the global
  `wt` check for `y`/`xy` splits because `tdump` is written as single-precision
  NetCDF output (`NF90_FLOAT`) even though the solver runs in `r8`
- the vegetation `treedump` comparison remains exact in the current one-step
  diagnostic

Status:

- classification: `local developer / diagnostic`
- role: solver-facing integration test
- expected runtime target: `hpc`
- shared fixture owners:
  - `tests/cases/100/` for the no-tree Xie/Castro-style comparison
  - `tests/cases/526/` for the tree-forcing comparison

Run it directly with:

```bash
python tests/integration/processor_boundaries/test_processor_boundaries.py
```
