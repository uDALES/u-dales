# IBM Cell Lookup Test

Validates `modibm::cell_index`, which locates the grid cell containing an IBM
wall-function reconstruction point.

## What it validates

`cell_index` is written as

```fortran
cell_index = count(p >= grid)
```

This is equivalent to the backward search it replaced,

```fortran
findloc(p >= grid, .true., 1, back=.true.)
```

only when `grid` is monotonically increasing and fully initialised over its whole
declared extent. Under that condition the mask `p >= grid` is `.true.` for a
leading run of elements, so the last `.true.` position equals the number of
`.true.` values.

`findloc` is not used because the NVHPC runtime has no implementation of it for
logical arrays. It compiles without complaint and aborts at run time with

```
0: FINDLOC: unimplemented for data type
```

The test checks three things:

1. **The precondition.** Every grid array (`xh`, `xf`, `yh`, `yf`, `zh`, `zf`) is
   strictly increasing across its full extent, including halo elements. A single
   uninitialised entry would break the leading-run property that `count` relies
   on, and would not otherwise be noticed.
2. **The equivalence.** `cell_index` agrees with an explicit backward search over
   a probe sweep: below the first node, exactly on each node, one ULP either side
   of each node, midway between nodes, and above the last node. Each probe also
   asserts the independently known expected index, so an off-by-one is caught
   even though both implementations would agree with each other.
3. **The real path.** Every reconstruction point that `initibm` produced for the
   loaded case resolves to the indices it stored.

## Why it uses case 100

Step 3 only means anything on a geometry that reaches the reconstruction branch
in `initibm`, which happens when

```fortran
log(bnddst / facz0) <= 1.
```

`tests/cases/100` (`xie_castro_2008_STL.stl`) has facet sections down to
`bnddst = 0.0719` against `z0 = 0.05`, so around 2300 sections reconstruct.

Grid-aligned geometries such as `tests/cases/101` have `bnddst >= 0.25 = 5*z0`,
never reach the branch, and would give a vacuous pass. The test therefore fails
if it finds no reconstruction sections rather than reporting success, so
swapping the fixture for an aligned case is caught rather than silently
reducing coverage to nothing.

## Requirements

- A CPU debug build: `tools/build_executable.sh common debug`
- MPI (`mpiexec`). Runs on a single rank.

## Running

```bash
tests/integration/ibm_cell_lookup/run_test.sh
```

Overridable via environment: `UDALES_BUILD`, `CASE_SOURCE`, `NAMELIST_SOURCE`,
`NAMELIST`, `NPROCS`, `MPIEXEC`.

## Implementation

- Runmode `1006` (`TEST_IBM_CELL_LOOKUP` in `src/modglobal.f90`)
- `tests_ibm_cell_lookup` in `src/tests.f90`, dispatched from
  `execute_runmode_actions` in `src/program.f90`

## See also

`tests/unit/fortran_intrinsics` probes the compiler for the intrinsic gap that
motivated this test, and fails if `src/` uses an intrinsic the target compiler
cannot run.
