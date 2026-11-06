# Warm-start facet continuity (issue #269)

Regression coverage for GitHub issue #269, *"Warm start of non-neutral flow does
not work. Facet temperature needs to be updated."*

## What this validates

Before the fix, `writerestartfiles` checkpointed only the fluid fields (`initd`)
and the scalars (`inits`). No facet state was ever written, so on
`lwarmstart = .true.` the solver called `readfacetfiles` first (which
unconditionally rebuilds the wall temperature profile from `Tfacinit.inp.<expnr>`
as a linear ramp to `bldT`/`flrT`) and `readrestartfiles` never touched a facet
array. `facT` — the only prognostic state the surface energy balance integrates —
was silently reset to cold-start values while the fluid continued from the
checkpoint, so the warm start was not a continuation.

The fix adds a third restart stream, `initf<ntrun>.<expnr>`. This suite asserts
that a restarted chain reproduces the reference chain's facet temperature
trajectory.

## Fixture

`tests/cases/064`, an energy-balance case with `lEB = .true.`,
`iwalltemp = 2`, `iwallmoist = 2`, `nfaclyrs = 5`, `nfcts = 140` and
`lwriteEBfiles = .true.`, so `facT.064.nc` is written and directly comparable
between runs.

Per `tests/README.md`, the fixture is not mutated. `namoptions.064.template`
holds the test-local namelist with explicit `runmode` and decomposition; the
driver substitutes `@...@` tokens into it. It overrides the fixture with
`nprocx`/`nprocy` of 2x2 and 2x1 (instead of 4x4), fixed time stepping
(`ladaptive = .false.`, `dtmax = 0.5`) so the checkpoint and the energy-balance
updates land on exact times, and disables the statistics dumps.

## Protocol

For each of the 2x2 and 2x1 decompositions:

1. **Reference run** — cold start, `runtime = 2T`, `trestart = T`. Writes
   restart files at `t = T` and `t = 2T`, and `facT.064.nc` over `[0, 2T]`.
2. **Restart run** — fresh staged directory containing the `t = T` restart
   files, `lwarmstart = .true.`, `startfile` pointing at the `t = T` `initd...`,
   `runtime = T`.

Assertions (`compare_facT.py`):

- **The #269 regression assertion**: the first `facT` record written after the
  restart equals the reference record at the same time. On unfixed code this
  instead equals the `Tfacinit.inp.064` profile advanced by one energy-balance
  step.
- Every overlapping output time matches the reference.
- The reference trajectory actually moves over the run, so the assertions above
  cannot become vacuous if the schedule is changed.
- **Negative control**: a third run with `lreadfacT = .false.` reproduces the
  pre-fix behaviour, and the suite asserts it does *not* match the reference.
  This keeps the suite self-validating without needing a second build.
- **MPI shape check**: the first post-restart record is identical between the
  2x2 and 2x1 restarts, which guards the facet broadcasts in
  `readrestartfiles`.

A third chain checkpoints **inside** an energy-balance interval
(`T_RESTART_MID = 12` with `dtEB = 5`, so the last EB update before the
checkpoint was at `t = 10`). Between EB updates the solver accumulates the facet
heat and moisture fluxes and divides them by the elapsed interval at the next
update; the `initf` file carries that interval state (`tEB`, `tnextEB`,
`fachfi`, `facefi`). The restarted run must therefore perform its first EB
update at the reference's `t = 15`, not at `12 + 5 = 17`, and match the
reference there. Without the interval state the restarted run's first `facT`
record would sit at `t = 17`, which the comparator rejects because the reference
has no record at that time.

Two further runs check the operational edges:

- **Missing `initf...` aborts**: with `lreadfacT = .true.` the facet
  initialiser does not read `Tfacinit*.inp` on a warm start, so a restart
  directory without the facet restart file has no facet state. The run must
  abort with an error naming the file and pointing at `lreadfacT = .false.`.
- **`lreadfacT = .false.`** must not read the file even when it is present.

Tolerances are set relative to single-precision resolution at ~300 K, because
`facT.064.nc` is written as `NF90_FLOAT` (see `modstat_nc.f90`); bitwise
equality is not required. In practice the restart currently reproduces the
reference exactly.

## Prerequisites

- A built `u-dales` executable. The path comes from `UDALES_BUILD` and defaults
  to `build/debug/u-dales`.
- A working MPI launcher. The suite runs at most 4 ranks.
- A `python3` with `netCDF4` and `numpy` for `compare_facT.py`.

**The solver runtime modules must match the stack the executable was built
with.** The repo is currently inconsistent about this:
`tools/build_executable.sh icl` uses the `intel/2021a` /
`netCDF-Fortran/4.5.3-iimpi-2021a` stack, while `tools/hpc_execute.sh` and the
sibling test suites use `intel/2025a` / `netCDF-Fortran/4.6.1-iimpi-2023a`.
Running an executable against the wrong netCDF-Fortran ABI segfaults inside
`nf_def_var`. `UDALES_RUNTIME_MODULES` defaults to the 2025a stack used by the
sibling suites; override it if you built with `tools/build_executable.sh icl`.

## Running

```bash
cd tests/integration/warmstart_facets
./run_test.sh
```

On a machine where the executable was built with `tools/build_executable.sh icl`:

```bash
UDALES_BUILD=$PWD/build/debug/u-dales \
UDALES_RUNTIME_MODULES="intel/2021a netCDF/4.8.0-iimpi-2021a netCDF-Fortran/4.5.3-iimpi-2021a FFTW/3.3.9-intel-2021a" \
UDALES_PYTHON_MODULES="netcdf4-python/1.6.4-foss-2023a" \
bash tests/integration/warmstart_facets/run_test.sh
```

Other knobs: `CASE_SOURCE`, `MPIEXEC`, `MPI_LAUNCH_EXTRA_ARGS`, `TMPDIR_PARENT`,
`PYTHON`, `KEEP_RUN_DIR=1` (keep the run directory on success too), and the
schedule variables `T_RESTART`, `T_RESTART_MID`, `T_TOTAL`, `DTEB`.

### GPU build

With `UDALES_GPU_BUILD` set to a GPU executable, the suite also runs, after the
CPU chains:

- a GPU warm start from the CPU checkpoint, compared with the CPU reference
  (the facet state restored on the host must be what the device computes with;
  before this check existed the device mirrors were filled before
  `readinitfiles` and a warm start ran on a surface temperature of zero);
- a GPU reference and a GPU warm start from the GPU checkpoint, plus the
  mid-interval variant, compared with each other for continuity and with the
  CPU reference for build parity (the device-side flux accumulators must reach
  the `initf` file and be picked up again).

Cross-build comparisons use `BUILD_PARITY_ATOL` in `compare_facT.py`, which is
the continuity tolerance: the CPU and GPU builds agree exactly on this case at
float32 output resolution, and a one-facet shift of the device heat flux
measured 1e-3 to 6e-3 K. Every rank needs a GPU of its own, so the GPU decomposition follows
`UDALES_GPU_MAX_RANKS` (`>= 4`: 2x2, `2`: 2x1, `1`: 1x1) and the launch goes
through `UDALES_GPU_LAUNCH_WRAPPER` (default `tools/bind.sh`). Further knobs:
`UDALES_GPU_MPIEXEC` (default `MPIEXEC`), `UDALES_GPU_RUNTIME_MODULES` (default
none; the GPU executable finds its libraries through its rpath),
`GPU_MPI_LAUNCH_EXTRA_ARGS`, and `UDALES_GPU_SELFTEST=1`, which switches the
CUDA device self-tests on in the GPU runs (Debug builds only). One of those
self-tests, `test_facet_state_initial_upload`, asserts that the device facet
state equals the host's after initialisation, and the Debug build keeps
asserting it at every `updateDevice` and at loop start.

The `gpu-mpi` group runs this with the Debug builds on two ranks and the
self-tests on; `gpu-full` runs it with the Release builds on four ranks. See
`tests/test_suites.yml`.

## CI policy

`class: experimental`, `kind: integration`, `component: solver`,
`platform: hpc`, `cost: medium`. Case 064 is a 64^3 energy-balance case and is
heavier than the current members of the `supported` merge gate, so this suite is
deliberately **not** in that gate. It runs eight short solver invocations.

On failure the run directory is preserved and its path printed.
